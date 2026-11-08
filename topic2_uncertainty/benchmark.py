"""Audit observed-label evaluation against known finite-reference expected loss."""

from itertools import combinations

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline

from research_common.metrics import (
    binary_scores,
    binomial_interval,
    cluster_mean_interval,
    rank_agreement,
    selected_indices,
    selection_regret,
    unit_losses,
)
from research_common.models import feature_columns, fit_checked, logistic_model
from research_common.randomness import derive_seed
from research_common.splits import assert_disjoint, group_split
from research_common.validation import integer, mapping, unique_choices

from .simulation import simulate, validate_covariates, validate_scenario

MODELS = ["unweighted_logistic", "confidence_weighted_logistic", "expression_only", "random_forest"]


def fit_model(name, train, seed):
    if name not in MODELS:
        raise ValueError(f"Unsupported model: {name}")
    prefixes = ["expr"] if name == "expression_only" else ["expr", "ptm", "technical"]
    columns = feature_columns(train, prefixes)
    if name == "random_forest":
        model = Pipeline(
            [
                ("impute", SimpleImputer(add_indicator=True, keep_empty_features=True)),
                (
                    "model",
                    RandomForestClassifier(n_estimators=80, min_samples_leaf=5, random_state=seed, n_jobs=1),
                ),
            ]
        )
    else:
        model = logistic_model(seed)
    fit_args = {}
    if name == "confidence_weighted_logistic":
        # A transparent baseline, not a calibrated PU estimator or a novel architecture.
        fit_args["model__sample_weight"] = np.clip(
            train.annotation_confidence.to_numpy() * train.evidence_confidence.to_numpy(), 0.05, 1
        )
    fit_checked(model, train[columns], train.label, **fit_args)
    return model, columns


def validate_config(config):
    mapping(config)
    unknown = set(config) - {
        "seed",
        "replicates",
        "n_donors",
        "reference_donors",
        "rows_per_donor",
        "bootstrap_draws",
        "models",
        "protocols",
        "scenarios",
        "save_unit_scores",
    }
    if unknown:
        raise ValueError(f"Unknown configuration fields: {sorted(unknown)}")
    integer(config.get("seed"), "seed", maximum=2**32 - 1)
    for key, minimum in [
        ("replicates", 2),
        ("bootstrap_draws", 50),
        ("n_donors", 8),
        ("reference_donors", 8),
        ("rows_per_donor", 4),
    ]:
        integer(config.get(key), key, minimum=minimum)
    unique_choices(config.get("models"), "models", MODELS, minimum=2)
    unique_choices(config.get("protocols"), "protocols", {"random_row", "held_out_donor"})
    if not isinstance(config.get("scenarios"), list) or not config["scenarios"]:
        raise ValueError("Nonempty scenarios list required")
    for scenario in config["scenarios"]:
        validate_scenario(scenario)
        if not isinstance(scenario.get("name"), str) or not scenario["name"].strip():
            raise ValueError("Each scenario needs a nonblank name")
    names = [s["name"] for s in config["scenarios"]]
    if len(set(names)) != len(names):
        raise ValueError("Scenario names must be unique")
    if not isinstance(config.get("save_unit_scores", True), bool):
        raise ValueError("save_unit_scores must be boolean")


def summarize(scores):
    records = []
    for key, subset in scores.groupby(["scenario", "protocol", "model"], sort=True):
        n = len(subset)
        coverage = float(subset.covers_reference.mean())
        low, high = binomial_interval(int(subset.covers_reference.sum()), n)
        records.append(
            dict(zip(["scenario", "protocol", "model"], key))
            | {
                "replicates": n,
                "mean_bias": subset.bias.mean(),
                "bias_mcse": subset.bias.std(ddof=1) / np.sqrt(n),
                "mean_absolute_error": subset.absolute_error.mean(),
                "absolute_error_mcse": subset.absolute_error.std(ddof=1) / np.sqrt(n),
                "rmse": np.sqrt(np.mean(subset.bias**2)),
                "reference_containment_rate": coverage,
                "reference_containment_mcse": np.sqrt(coverage * (1 - coverage) / n),
                "reference_containment_low": low,
                "reference_containment_high": high,
            }
        )
    return pd.DataFrame(records)


def paired_comparisons(scores):
    records = []
    for (scenario, protocol), subset in scores.groupby(["scenario", "protocol"]):
        wide = subset.pivot(index="replicate", columns="model", values="absolute_error")
        for a, b in combinations(sorted(wide.columns), 2):
            delta = wide[a] - wide[b]
            records.append(
                {
                    "scenario": scenario,
                    "protocol": protocol,
                    "model_a": a,
                    "model_b": b,
                    "mean_absolute_error_difference_a_minus_b": delta.mean(),
                    "paired_mcse": delta.std(ddof=1) / np.sqrt(len(delta)),
                    "replicates": len(delta),
                }
            )
    return pd.DataFrame(records)


def run_benchmark(config, source=None, progress=None):
    validate_config(config)
    if source is not None:
        validate_covariates(source)
    rows, splits, units = [], [], []
    for scenario in config["scenarios"]:
        for rep in range(config["replicates"]):
            if progress:
                progress(f"Topic 2: {scenario['name']}, replicate {rep + 1}/{config['replicates']}")
            seed = derive_seed(config["seed"], scenario["name"], rep, "development")
            reference_seed = derive_seed(config["seed"], scenario["name"], rep, "reference")
            development = simulate(seed, scenario, config["n_donors"], config["rows_per_donor"], source)
            reference = simulate(
                reference_seed,
                scenario,
                config["reference_donors"],
                config["rows_per_donor"],
                source,
                prefix="ORACLE",
            )
            for protocol in config["protocols"]:
                donor_group = "donor_id" if source is None else "source_donor_id"
                group = "sample_id" if protocol == "random_row" else donor_group
                train, test = group_split(development, group, seed, 0.3)
                if protocol == "held_out_donor":
                    assert_disjoint(train, test, ["donor_id", donor_group])
                for partition, subset in [("train", train), ("test", test)]:
                    splits.extend(
                        {
                            "scenario": scenario["name"],
                            "replicate": rep,
                            "protocol": protocol,
                            "partition": partition,
                            **record,
                        }
                        for record in subset[
                            ["sample_id", "donor_id", "source_donor_id", "source_sample_id"]
                        ].to_dict("records")
                    )
                for name in config["models"]:
                    model, cols = fit_model(name, train, seed)
                    p_test = model.predict_proba(test[cols])[:, 1]
                    p_reference = model.predict_proba(reference[cols])[:, 1]
                    loss = (test.label.to_numpy() - p_test) ** 2
                    estimate, low, high = cluster_mean_interval(
                        loss,
                        test.donor_id,
                        seed,
                        config["bootstrap_draws"],
                        resample_groups=test[donor_group],
                    )
                    q = reference.truth_probability.to_numpy()
                    # E[(Y-p)^2 | q] = q(1-p)^2 + (1-q)p^2. No sampled outcome noise.
                    expected_loss = q * (1 - p_reference) ** 2 + (1 - q) * p_reference**2
                    reference_units = unit_losses(expected_loss, reference.donor_id)
                    test_units = unit_losses(loss, test.donor_id)
                    truth = float(reference_units.mean())
                    if config.get("save_unit_scores", True):
                        for partition, losses in [("test", test_units), ("reference", reference_units)]:
                            units.extend(
                                {
                                    "scenario": scenario["name"],
                                    "replicate": rep,
                                    "protocol": protocol,
                                    "model": name,
                                    "partition": partition,
                                    "donor_id": donor,
                                    "mean_loss": float(value),
                                }
                                for donor, value in losses.items()
                            )
                    row = {
                        "scenario": scenario["name"],
                        "replicate": rep,
                        "seed": seed,
                        "protocol": protocol,
                        "model": name,
                        "reference_seed": reference_seed,
                        "resampling_unit": donor_group,
                        "estimated_brier": estimate,
                        "n_train_donors": train.donor_id.nunique(),
                        "n_test_donors": test.donor_id.nunique(),
                        "n_test_resampling_blocks": test[donor_group].nunique(),
                        "reference_mcse_conditional": reference_units.std(ddof=1)
                        / np.sqrt(len(reference_units)),
                        "reference_expected_brier": truth,
                        "bias": estimate - truth,
                        "absolute_error": abs(estimate - truth),
                        "ci_low": low,
                        "ci_high": high,
                        "covers_reference": bool(low <= truth <= high),
                        "test_truth_brier_diagnostic": float(
                            unit_losses((test.truth_label.to_numpy() - p_test) ** 2, test.donor_id).mean()
                        ),
                        "overlapping_donors": len(set(train.donor_id) & set(test.donor_id)),
                        "overlapping_source_donors": 0
                        if source is None
                        else len(set(train.source_donor_id) & set(test.source_donor_id)),
                        "reference_source_overlap": 0
                        if source is None
                        else len(set(development.source_donor_id) & set(reference.source_donor_id)),
                    }
                    row.update({"observed_" + k: v for k, v in binary_scores(test.label, p_test).items()})
                    rows.append(row)
    scores = pd.DataFrame(rows)
    audit = []
    for (scenario, rep, protocol), subset in scores.groupby(["scenario", "replicate", "protocol"]):
        estimates, truths = subset.estimated_brier.to_numpy(), subset.reference_expected_brier.to_numpy()
        selected = selected_indices(estimates)
        audit.append(
            {
                "scenario": scenario,
                "replicate": rep,
                "protocol": protocol,
                "rank_agreement": rank_agreement(estimates, truths),
                "selection_regret": selection_regret(estimates, truths),
                "selection_tie_count": len(selected),
                "selected_best_probability": float(np.mean(truths[selected] == truths.min())),
            }
        )
    return {
        "scores": scores,
        "meta_evaluation": pd.DataFrame(audit),
        "summary": summarize(scores),
        "paired_comparisons": paired_comparisons(scores),
        "unit_scores": pd.DataFrame(
            units,
            columns=["scenario", "replicate", "protocol", "model", "partition", "donor_id", "mean_loss"],
        ),
        "split_assignments": pd.DataFrame(splits),
    }
