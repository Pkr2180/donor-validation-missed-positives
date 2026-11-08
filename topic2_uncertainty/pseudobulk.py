"""Aggregate raw counts by donor and annotated cell type without cell pseudoreplication."""

import numpy as np

from research_common.validation import frame_contract, identifiers, numeric, real


def aggregate_counts(counts, annotations, minimum_confidence=0.0):
    real(minimum_confidence, "minimum_confidence", maximum=1.0)
    frame_contract(counts, ["cell_id", "gene_id", "count"])
    frame_contract(annotations, ["cell_id", "donor_id", "cell_type", "confidence"])
    identifiers(counts, ["cell_id", "gene_id"])
    identifiers(annotations, ["cell_id", "donor_id", "cell_type"])
    counts = numeric(counts[["cell_id", "gene_id", "count"]], ["count"])
    annotations = numeric(annotations[["cell_id", "donor_id", "cell_type", "confidence"]], ["confidence"])
    if not {"cell_id", "gene_id", "count"} <= set(counts):
        raise ValueError("counts needs cell_id,gene_id,count")
    needed = {"cell_id", "donor_id", "cell_type", "confidence"}
    if not needed <= set(annotations):
        raise ValueError("annotations needs cell_id,donor_id,cell_type,confidence")
    if (
        counts[["cell_id", "gene_id", "count"]].isna().any().any()
        or annotations[list(needed)].isna().any().any()
    ):
        raise ValueError("Missing counts or annotation metadata")
    if counts.duplicated(["cell_id", "gene_id"]).any() or annotations.cell_id.duplicated().any():
        raise ValueError("Duplicate cell/gene pairs or annotation cell IDs")
    values = counts["count"].to_numpy(float)
    confidence = annotations.confidence.to_numpy(float)
    if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
        raise ValueError("Use nonnegative integer raw counts, not normalized expression")
    if not np.isfinite(confidence).all() or ((confidence < 0) | (confidence > 1)).any():
        raise ValueError("Annotation confidence must be finite and in [0,1]")
    if not set(counts.cell_id) <= set(annotations.cell_id):
        raise ValueError("Some counted cells lack annotations")
    if set(counts.cell_id) != set(annotations.cell_id):
        raise ValueError("Annotated cells without counts: supply explicit zero rows for zero-library cells")
    # Accumulate Python integers to prevent silent int64 overflow.
    counts["count"] = counts["count"].map(int).astype(object)
    selected = annotations[annotations.confidence >= minimum_confidence]
    merged = counts.merge(selected, on="cell_id", validate="many_to_one")
    if merged.empty:
        raise ValueError("No cells remain after confidence filtering")
    result = merged.groupby(["donor_id", "cell_type", "gene_id"], as_index=False)["count"].sum()
    if result["count"].gt(np.iinfo(np.int64).max).any():
        raise ValueError("Aggregated counts exceed int64 range")
    result["count"] = result["count"].astype(np.int64)
    qc = (
        annotations.assign(retained=annotations.confidence >= minimum_confidence)
        .groupby(["donor_id", "cell_type"], as_index=False)
        .agg(
            n_cells=("cell_id", "count"),
            retained_cells=("retained", "sum"),
            mean_confidence=("confidence", "mean"),
        )
    )
    return result, qc
