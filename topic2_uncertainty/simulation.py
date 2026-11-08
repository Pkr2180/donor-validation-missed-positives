"""Known-truth simulations, optionally retaining empirical donor feature blocks.

True labels and probabilities are never features. Plasmode truth refers only to
the imposed mechanism, not to genuine biological PTM interactions.
"""

import numpy as np
import pandas as pd
from scipy.special import expit

from research_common.randomness import derive_seed
from research_common.validation import frame_contract, identifiers, integer, mapping, real


def validate_covariates(frame):
    required = {"source_donor_id", "source_sample_id"}
    frame_contract(frame, required)
    identifiers(frame, required)
    if not required <= set(frame) or frame[list(required)].isna().any().any():
        raise ValueError("Plasmode input needs source_donor_id and source_sample_id")
    if frame.source_sample_id.duplicated().any() or frame.source_donor_id.nunique() < 8:
        raise ValueError("Unique source samples and at least eight source donors required")
    columns = sorted(c for c in frame if c.startswith(("expr__", "ptm__")))
    for prefix in ["expr__", "ptm__"]:
        if len([c for c in columns if c.startswith(prefix)]) < 2:
            raise ValueError(f"At least two numeric {prefix} covariates required")
    if not np.isfinite(frame[columns].to_numpy(float)).all():
        raise ValueError("Plasmode source covariates must be finite; document prior preprocessing")
    return columns


def validate_scenario(scenario):
    mapping(scenario, "scenario")
    unknown = set(scenario) - {
        "name",
        "mechanism",
        "annotation_error",
        "false_negative_rate",
        "ptm_missingness",
        "donor_effect",
    }
    if unknown:
        raise ValueError(f"Unknown scenario fields: {sorted(unknown)}")
    for key in ["annotation_error", "false_negative_rate", "ptm_missingness"]:
        real(scenario.get(key), key, maximum=1.0, upper_open=True)
    real(scenario.get("donor_effect"), "donor_effect")
    if scenario.get("mechanism") not in {"linear", "interaction", "ptm_dominant"}:
        raise ValueError("Unsupported simulation mechanism")


def simulate(seed, scenario, n_donors=24, rows_per_donor=30, source=None, prefix="DEV"):
    integer(n_donors, "n_donors", minimum=8)
    integer(rows_per_donor, "rows_per_donor", minimum=4)
    validate_scenario(scenario)
    rng = np.random.default_rng(derive_seed(seed, "latent_covariates"))
    outcome_rng = np.random.default_rng(derive_seed(seed, "outcomes"))
    annotation_rng = np.random.default_rng(derive_seed(seed, "annotation"))
    missing_rng = np.random.default_rng(derive_seed(seed, "missingness"))
    label_rng = np.random.default_rng(derive_seed(seed, "label_observation"))
    technical_rng = np.random.default_rng(derive_seed(seed, "technical_noise"))
    if source is not None:
        source_columns = validate_covariates(source)
        source = source.sort_values(["source_donor_id", "source_sample_id"]).copy()
        source[source_columns] = source[source_columns].astype(float)
        location = source[source_columns].mean()
        scale = source[source_columns].std().replace(0, 1).fillna(1)
        source_groups = dict(tuple(source.groupby("source_donor_id", sort=True)))
    rows = []
    for donor in range(n_donors):
        donor_effect = rng.normal(scale=scenario["donor_effect"])
        donor_proxy = rng.normal(size=4)
        if source is None:
            expr = rng.normal(size=(rows_per_donor, 3))
            ptm = rng.normal(size=(rows_per_donor, 2))
            source_id = "FULLY_SYNTHETIC"
            source_sample_ids = ["FULLY_SYNTHETIC"] * rows_per_donor
        else:
            source_id = rng.choice(list(source_groups))
            block = source_groups[source_id]
            sampled = block.iloc[rng.integers(0, len(block), rows_per_donor)]
            source_sample_ids = sampled.source_sample_id.to_numpy()
            z = (sampled[source_columns] - location) / scale
            expr = z[[c for c in source_columns if c.startswith("expr__")]].to_numpy()
            ptm = z[[c for c in source_columns if c.startswith("ptm__")]].to_numpy()
        eta = 0.9 * expr[:, 0] - 0.5 * expr[:, 1] + 0.8 * ptm[:, 0] + donor_effect - 0.4
        if scenario["mechanism"] == "interaction":
            eta += 1.1 * expr[:, 0] * ptm[:, 1]
        if scenario["mechanism"] == "ptm_dominant":
            eta = 1.8 * ptm[:, 0] - 0.8 * ptm[:, 1] + 0.2 * expr[:, 0] + donor_effect - 0.4
        truth_p = expit(eta)
        truth_y = (outcome_rng.random(rows_per_donor) < truth_p).astype(int)
        observed_y = truth_y * (label_rng.random(rows_per_donor) >= scenario["false_negative_rate"])
        # Assumed/oracle-assisted reliability, not an empirically calibrated score.
        error = scenario["annotation_error"]
        reliability = (
            np.ones(rows_per_donor)
            if error == 0
            else annotation_rng.beta((1 - error) * 20, error * 20, rows_per_donor)
        )
        corruption = annotation_rng.random(rows_per_donor) > reliability
        observed_expr = expr.copy()
        replacements = annotation_rng.normal(size=observed_expr.shape)
        observed_expr[corruption] = replacements[corruption]
        observed_ptm = ptm.copy()
        observed_ptm[missing_rng.random(observed_ptm.shape) < scenario["ptm_missingness"]] = np.nan
        for j in range(rows_per_donor):
            row = {
                "sample_id": f"{prefix}_{seed}_{donor}_{j}",
                "donor_id": f"{prefix}_{seed}_D{donor}",
                "source_donor_id": source_id,
                "source_sample_id": source_sample_ids[j],
                "label": int(observed_y[j]),
                "annotation_corrupted": bool(corruption[j]),
                "truth_label": int(truth_y[j]),
                "truth_probability": truth_p[j],
                "annotation_confidence": reliability[j],
                "evidence_confidence": 1.0 if observed_y[j] else 1 - scenario["false_negative_rate"],
            }
            row.update({f"expr__{k}": v for k, v in enumerate(observed_expr[j])})
            row.update({f"ptm__{k}": v for k, v in enumerate(observed_ptm[j])})
            row.update(
                {f"technical__{k}": donor_proxy[k] + technical_rng.normal(scale=0.05) for k in range(4)}
            )
            rows.append(row)
    return pd.DataFrame(rows)
