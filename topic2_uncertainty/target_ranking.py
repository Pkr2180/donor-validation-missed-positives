"""Measurement-uncertainty propagation for exploratory target network rankings."""

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from research_common.validation import frame_contract, identifiers, integer, numeric


def rank_targets(nodes, edges, draws=1000, seed=19, top_k=None):
    """Monte Carlo edge support, annotation confidence and expression measurement.

    Beta(alpha,beta) is an explicit uncertainty assumption. It is not automatically
    a calibrated posterior. No confidence interval here establishes causal targets.
    """
    nr = {"target_id", "expression_mean", "expression_sd", "annotation_alpha", "annotation_beta"}
    er = {"source", "target", "support_alpha", "support_beta", "evidence_tier"}
    frame_contract(nodes, nr)
    frame_contract(edges, er, allow_empty=True)
    identifiers(nodes, ["target_id"])
    identifiers(edges, ["source", "target", "evidence_tier"])
    integer(draws, "draws", minimum=50)
    integer(seed, "seed", maximum=2**32 - 1)
    nodes = numeric(nodes, sorted(nr - {"target_id"})).sort_values("target_id").reset_index(drop=True)
    edges = numeric(edges, ["support_alpha", "support_beta"]).copy()
    # The score uses an undirected network, so A-B and B-A are the same edge.
    pairs = [tuple(sorted((a, b))) for a, b in zip(edges.source, edges.target)]
    if len(set(pairs)) != len(pairs):
        raise ValueError("Reconcile duplicate undirected evidence edges")
    if pairs:
        edges["source"], edges["target"] = zip(*pairs)
    edges = edges.sort_values(["source", "target"]).reset_index(drop=True)
    top_k = min(10, max(1, len(nodes) // 5)) if top_k is None else top_k
    integer(top_k, "top_k", minimum=1, maximum=len(nodes))
    if not nr <= set(nodes) or not er <= set(edges):
        raise ValueError("Missing node or evidence-edge schema fields")
    if draws < 50 or nodes.empty or nodes.target_id.duplicated().any():
        raise ValueError("Unique nonempty targets and at least 50 draws required")
    if not set(edges.source).union(edges.target) <= set(nodes.target_id):
        raise ValueError("Edges reference missing targets")
    if edges.duplicated(["source", "target"]).any() or (edges.source == edges.target).any():
        raise ValueError("Reconcile duplicate evidence and remove self edges before ranking")
    if not set(edges.evidence_tier) <= {"predicted", "independent_support", "experimental"}:
        raise ValueError("Unsupported evidence tier")
    for frame, columns in [(nodes, list(nr - {"target_id"})), (edges, ["support_alpha", "support_beta"])]:
        if not np.isfinite(frame[columns].to_numpy(float)).all():
            raise ValueError("Nonfinite uncertainty inputs")
    if (nodes[["expression_mean", "expression_sd"]] < 0).any().any():
        raise ValueError("Expression inputs must be nonnegative and normalized comparably")
    if (nodes[["annotation_alpha", "annotation_beta"]] <= 0).any().any() or (
        edges[["support_alpha", "support_beta"]] <= 0
    ).any().any():
        raise ValueError("Beta distribution parameters must be positive")
    rng = np.random.default_rng(seed)
    lookup = {t: i for i, t in enumerate(nodes.target_id)}
    expression = np.maximum(
        0, rng.normal(nodes.expression_mean.to_numpy(), nodes.expression_sd.to_numpy(), (draws, len(nodes)))
    )
    annotation = rng.beta(
        nodes.annotation_alpha.to_numpy(), nodes.annotation_beta.to_numpy(), (draws, len(nodes))
    )
    degree = np.zeros((draws, len(nodes)))
    for edge in edges.itertuples():
        support = rng.beta(edge.support_alpha, edge.support_beta, draws)
        degree[:, lookup[edge.source]] += support
        degree[:, lookup[edge.target]] += support
    score = expression * annotation * (1 + degree)
    if not np.isfinite(score).all():
        raise ValueError("Score overflow; rescale expression inputs")
    ranks = rankdata(-score, axis=1, method="average")
    first = rankdata(-score, axis=1, method="min")
    last = rankdata(-score, axis=1, method="max")
    # Uniform tie-breaking at the top-k boundary, evaluated analytically.
    inclusion = np.clip((top_k - first + 1) / (last - first + 1), 0, 1)
    low, median, high = np.quantile(score, [0.025, 0.5, 0.975], axis=0)
    result = pd.DataFrame(
        {
            "target_id": nodes.target_id,
            "score_median": median,
            "score_low": low,
            "score_high": high,
            "mean_rank": ranks.mean(axis=0),
            "top_k": top_k,
            "top_k_probability": inclusion.mean(axis=0),
            "top_k_mcse": inclusion.std(axis=0, ddof=1) / np.sqrt(draws),
            "mean_rank_mcse": ranks.std(axis=0, ddof=1) / np.sqrt(draws),
            "draws": draws,
        }
    )
    result["interpretation"] = "EXPLORATORY NETWORK RANK; NOT CAUSAL OR PTM VALIDATION"
    return result.sort_values("score_median", ascending=False)
