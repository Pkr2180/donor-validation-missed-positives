"""Paper 2: compare the primary grid with the sweeps at identical settings (independent replicate sets) and report
Monte Carlo standard errors of mean rank agreement. Writes sensitivity/tables/tableS_replicate_set_consistency.csv
and tableS_rank_mcse_primary.csv. Usage (from bioinformatics_method_projects/): python modal_jobs/replicate_consistency_t2.py"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "outputs/modal_t2/uncertainty_research_sim_r1"
SENS = ROOT / "outputs/sens_t2/sens_t2"
EXTRA = ROOT / "outputs/extra_t2/extra_t2"
OUT = ROOT / "topic2_uncertainty/results/modal_research_r1/sensitivity/tables"


def meta(df):
    rows = []
    for (sc, rep, pr), g in df.groupby(["scenario", "replicate", "protocol"]):
        a, b = g.estimated_brier.to_numpy(), g.reference_expected_brier.to_numpy()
        rows.append({"scenario": sc, "replicate": rep, "protocol": pr, "rho": spearmanr(a, b).statistic if np.ptp(a) > 0 and np.ptp(b) > 0 else np.nan})
    return pd.DataFrame(rows)


def pooled(a, b):
    sa, sb = pd.read_csv(SENS / a / "scores.csv"), pd.read_csv(EXTRA / b / "scores.csv")
    return pd.concat([sa, sb.assign(replicate=sb.replicate + sa.replicate.max() + 1)])


raw = pd.read_csv(RAW / "scores.csv")
mp = meta(raw)
rows = []
for prim, (a, b), lv, label in [("clean_linear", ("fn_sweep", "fn_linear_b"), "0.0", "clean (false-negative sweep, rate 0)"),
                                 ("clean_linear", ("donor_sweep", "donor_linear_b"), "0.5", "clean (donor-effect sweep, SD 0.5)"),
                                 ("negative_label_only", ("fn_sweep", "fn_linear_b"), "0.4", "false negatives 0.4"),
                                 ("donor_only", ("donor_sweep", "donor_linear_b"), "2.0", "strong donor effect, SD 2.0")]:
    sw = pooled(a, b)
    sw = sw[sw.scenario.str.endswith("_" + lv)]
    ms = meta(sw)
    for pr in ["random_row", "held_out_donor"]:
        x = mp[(mp.scenario == prim) & (mp.protocol == pr)].rho.dropna()
        y = ms[ms.protocol == pr].rho.dropna()
        p = raw[(raw.scenario == prim) & (raw.protocol == pr)]
        s = sw[sw.protocol == pr]
        se = np.sqrt(x.var() / len(x) + y.var() / len(y))
        rows.append({"setting": label, "protocol": pr, "rank_primary": x.mean(), "rank_mcse_primary": x.std() / np.sqrt(len(x)), "rank_sweep": y.mean(),
                     "rank_mcse_sweep": y.std() / np.sqrt(len(y)), "rank_diff": x.mean() - y.mean(), "rank_diff_z": (x.mean() - y.mean()) / se,
                     "bias_primary": p.bias.mean(), "bias_sweep": s.bias.mean(), "containment_primary": p.covers_reference.mean(), "containment_sweep": s.covers_reference.mean(),
                     "replicates_primary": len(x), "replicates_sweep": len(y)})
t = pd.DataFrame(rows)
t.to_csv(OUT / "tableS_replicate_set_consistency.csv", index=False)
mc = mp.groupby(["scenario", "protocol"]).rho.agg(["mean", "std", "count"]).reset_index()
mc["mcse"] = mc["std"] / np.sqrt(mc["count"])
mc.to_csv(OUT / "tableS_rank_mcse_primary.csv", index=False)
print(t.round(4).to_string())
print(mc.round(4).to_string())
