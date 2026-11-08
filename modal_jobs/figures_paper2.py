"""Publication figures for Paper 2 (missed positives and donor-grouped validation), redrawn from the raw per-replicate run files.
Intervals are mean Â± 1.96 Monte Carlo standard errors across replicates (per-replicate values averaged over the four models).
Usage (from modal_jobs/): python figures_paper2.py     -> topic2_uncertainty/results/modal_research_r1/figures_pub/"""
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

import pubfig as pf
from pubfig import plt

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "topic2_uncertainty/results/modal_research_r1"
S = R / "sensitivity"
RAW = ROOT / "outputs/modal_t2/uncertainty_research_sim_r1"
SENS = ROOT / "outputs/sens_t2/sens_t2"
EXTRA = ROOT / "outputs/extra_t2/extra_t2"
OUT = R / "figures_pub"
PR = ["held_out_donor", "random_row"]
SC = ["clean_linear", "annotation_only", "missing_ptm_only", "negative_label_only", "donor_only", "combined_nonlinear"]
SCL = {"clean_linear": "Clean", "annotation_only": "Covariate\ncorruption", "missing_ptm_only": "PTM\nmissingness", "negative_label_only": "False\nnegatives",
       "donor_only": "Strong donor\neffect", "combined_nonlinear": "Combined"}
MODELS = ["unweighted_logistic", "confidence_weighted_logistic", "expression_only", "random_forest"]
SH = pf.design_handles(PR, pf.SPLIT_C, pf.SPLIT_L, pf.SPLIT_M)
OFF = {"held_out_donor": -0.12, "random_row": 0.12}
made = {}


def lvl(s):
    return s.str.split("_").str[-1].astype(float)


def meta(df):
    rows = []
    for (sc, rep, pr), g in df.groupby(["scenario", "replicate", "protocol"]):
        a, b = g.estimated_brier.to_numpy(), g.reference_expected_brier.to_numpy()
        from scipy.stats import spearmanr
        rho = spearmanr(a, b).statistic if np.ptp(a) > 0 and np.ptp(b) > 0 else np.nan
        sel = a == a.min()
        rows.append({"scenario": sc, "replicate": rep, "protocol": pr, "rho": rho, "hit": float(np.mean(b[sel] == b.min())), "regret": b[sel].mean() - b.min()})
    return pd.DataFrame(rows)


def per_rep(df, extra_cols=()):
    """Per replicate and split: bias and containment averaged over models, merged with rank agreement and selection outcomes."""
    cols = {"bias": ("bias", "mean"), "cont": ("covers_reference", "mean")}
    for c in extra_cols:
        cols[c] = (c, "mean")
    a = df.groupby(["scenario", "replicate", "protocol"]).agg(**cols).reset_index()
    return a.merge(meta(df), on=["scenario", "replicate", "protocol"])


def mci(v):
    v = np.asarray(v, float)
    v = v[~np.isnan(v)]
    return v.mean(), 1.96 * v.std(ddof=1) / np.sqrt(len(v))


def pooled(a, b):
    sa, sb = pd.read_csv(SENS / a / "scores.csv"), pd.read_csv(EXTRA / b / "scores.csv")
    sb = sb.assign(replicate=sb.replicate + sa.replicate.max() + 1)
    return pd.concat([sa, sb], ignore_index=True)


raw = pd.read_csv(RAW / "scores.csv")
sk = pd.read_csv(SENS / "prevalence_primary/scores_with_skill.csv")
sk["skill_bias"] = (1 - sk.estimated_brier / (sk.observed_prevalence * (1 - sk.observed_prevalence))) - (1 - sk.reference_expected_brier / (sk.reference_mean_risk * (1 - sk.reference_mean_risk)))
P = per_rep(raw)
SKR = sk.groupby(["scenario", "replicate", "protocol"]).skill_bias.mean().reset_index()


def scen_panel(ax, frame, col, ylabel, ref=None, sc=SC):
    x = np.arange(len(sc))
    for p in PR:
        mu, ci = zip(*[mci(frame[(frame.scenario == s) & (frame.protocol == p)][col]) for s in sc])
        ax.errorbar(x + OFF[p], mu, yerr=ci, fmt=pf.SPLIT_M[p], color=pf.SPLIT_C[p], ms=4, lw=1)
    if ref is not None:
        pf.refline(ax, y=ref)
    ax.set_xticks(x, [SCL[s].replace("\n", " ") for s in sc], fontsize=6, rotation=30, ha="right", rotation_mode="anchor")
    ax.set_ylabel(ylabel)


# ------------------------------------------------------------------ Fig 1 study design
def fig1():
    fig = plt.figure(figsize=(pf.DOUBLE, 3.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")
    pf.box(ax, (0.01, 0.52), 0.23, 0.44, "Data-generating mechanism\nDonor effect ~ N(0, sÂ²) on log-odds\nExpression e1â€“e3, PTM p1â€“p2\n4 donor-level technical features\nOutcome: linear, interaction\nor PTM-dominant log-odds\nLatent risk q known exactly")
    pf.box(ax, (0.01, 0.06), 0.23, 0.38, "Defects (varied one at a time)\nFalse negatives: rate 0â€“0.5\nCovariate corruption: 0â€“0.4\nPTM missingness: 0.4\nDonor-effect SD: 0.5â€“3.0", fc="#fdf0ea")
    pf.box(ax, (0.29, 0.62), 0.20, 0.34, "Development data\n40 donors Ã— 40 rows\n(20â€“160 donors in a sweep)\nobserved labels only")
    pf.box(ax, (0.29, 0.06), 0.20, 0.44, "Reference population\n300 donors Ã— 40 rows\nsame scenario, independent\nexpected Brier loss against\nthe true risk q:\nq(1âˆ’p)Â² + (1âˆ’q)pÂ²")
    pf.box(ax, (0.54, 0.62), 0.20, 0.34, "Evaluation split\nRandom-row: 30% of rows\nHeld-out-donor: 12 of 40\ndonors tested", fc="#eaf1fb")
    pf.box(ax, (0.54, 0.06), 0.20, 0.44, "Four models\nUnweighted logistic\nConfidence-weighted logistic\n(oracle-assisted)\nExpression-only logistic\nRandom forest")
    pf.box(ax, (0.79, 0.06), 0.20, 0.90, "Performance measures\n(per replicate; 100â€“200\nreplicates per setting)\n\nBias: estimate âˆ’ reference\nContainment of the reference\nby 95% donor-bootstrap\nintervals\nRank agreement: Spearman\ncorrelation of estimated\nand reference losses\nSelection regret and\nP(best model chosen)\nSkill bias vs a constant\nforecast", fc="#eef8f3")
    pf.arrow(ax, (0.24, 0.79), (0.29, 0.79))
    pf.arrow(ax, (0.24, 0.25), (0.29, 0.25))
    pf.arrow(ax, (0.125, 0.44), (0.125, 0.52))
    pf.arrow(ax, (0.49, 0.79), (0.54, 0.79))
    pf.arrow(ax, (0.64, 0.62), (0.64, 0.50))
    pf.arrow(ax, (0.49, 0.28), (0.54, 0.28))
    pf.arrow(ax, (0.74, 0.40), (0.79, 0.40))
    pf.arrow(ax, (0.74, 0.79), (0.79, 0.79))
    made["fig1"] = pf.save(fig, OUT / "Fig1_study_design")


# ------------------------------------------------------------------ Fig 2 primary grid
def fig2():
    fig, axs = plt.subplots(2, 2, figsize=(pf.DOUBLE, 4.8), gridspec_kw={"hspace": 0.75, "wspace": 0.3})
    scen_panel(axs[0, 0], P, "bias", "Bias of Brier estimate", ref=0)
    scen_panel(axs[0, 1], P, "cont", "Containment of reference", ref=0.95)
    scen_panel(axs[1, 0], P, "rho", "Rank agreement (Spearman)", ref=0)
    scen_panel(axs[1, 1], P, "hit", "P(reference-best model chosen)")
    axs[0, 1].text(5.45, 0.95, "0.95", fontsize=6, color=pf.INK2, va="bottom", ha="right")
    for a_, l_, t_ in [(axs[0, 0], "A", "Bias"), (axs[0, 1], "B", "Interval containment"), (axs[1, 0], "C", "Rank agreement"), (axs[1, 1], "D", "Selection of the best model")]:
        pf.panel(a_, l_, t_, x=-0.15)
    fig.legend(handles=SH, loc="upper center", bbox_to_anchor=(0.5, -0.01), ncol=2)
    made["fig2"] = pf.save(fig, OUT / "Fig2_primary_grid")


# ------------------------------------------------------------------ Fig 3 Brier bias versus skill bias
def fig3():
    fig, axs = plt.subplots(1, 2, figsize=(pf.DOUBLE, 2.5), gridspec_kw={"wspace": 0.3})
    scen_panel(axs[0], P, "bias", "Bias of Brier estimate", ref=0)
    scen_panel(axs[1], SKR, "skill_bias", "Bias of Brier skill", ref=0)
    pf.panel(axs[0], "A", "Brier loss", x=-0.15)
    pf.panel(axs[1], "B", "Brier skill vs constant forecast", x=-0.15)
    fig.legend(handles=SH, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2)
    made["fig3"] = pf.save(fig, OUT / "Fig3_skill_vs_brier")


def sweep_fig(df, xlabel, name, letters="ABCD", xticks=None):
    d = per_rep(df)
    d["level"] = lvl(d.scenario)
    levels = sorted(d.level.unique())
    step = (levels[1] - levels[0]) * 0.08
    fig, axs = plt.subplots(1, 4, figsize=(pf.DOUBLE, 2.25), gridspec_kw={"wspace": 0.55})
    for ax, col, lab, ref, ttl in [(axs[0], "bias", "Bias of Brier estimate", 0, "Bias"), (axs[1], "cont", "Containment", 0.95, "Containment"),
                                   (axs[2], "rho", "Rank agreement", 0, "Rank agreement"), (axs[3], "hit", "P(best model chosen)", None, "Best model chosen")]:
        for p in PR:
            mu, ci = zip(*[mci(d[(d.level == lv) & (d.protocol == p)][col]) for lv in levels])
            ax.errorbar(np.array(levels) + (-step if p == "held_out_donor" else step), mu, yerr=ci, color=pf.SPLIT_C[p], marker=pf.SPLIT_M[p], ms=3.5, lw=1, elinewidth=0.7)
        if ref is not None:
            pf.refline(ax, y=ref)
        ax.set_xticks(levels if xticks is None else xticks), ax.set_xlabel(xlabel), ax.set_ylabel(lab)
        pf.panel(ax, letters["bias cont rho hit".split().index(col)], ttl, x=-0.32)
    fig.legend(handles=SH, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2)
    made[name] = pf.save(fig, OUT / name)
    return d


def fig4():
    sweep_fig(pooled("fn_sweep", "fn_linear_b"), "False-negative rate", "Fig4_false_negative_sweep")


def fig5():
    sweep_fig(pooled("donor_sweep", "donor_linear_b"), "Donor-effect SD", "Fig5_donor_effect_sweep", xticks=[0.5, 1.0, 1.5, 2.0, 3.0])


# ------------------------------------------------------------------ Fig 6 outcome mechanisms (summary table; 100-200 replicates)
def fig6():
    mt = pd.read_csv(S / "tables/tableS12_mechanism_sweeps.csv")
    fig, axs = plt.subplots(1, 4, figsize=(pf.DOUBLE, 2.3), gridspec_kw={"wspace": 0.5})
    spec = [(axs[0], "fn", "held_out_donor", "rank_agreement", "False-negative rate", "Rank agreement", "Held-out-donor"),
            (axs[1], "fn", "held_out_donor", "containment", "False-negative rate", "Containment", "Held-out-donor"),
            (axs[2], "donor", "random_row", "rank_agreement", "Donor-effect SD", "Rank agreement", "Random-row"),
            (axs[3], "donor", "held_out_donor", "rank_agreement", "Donor-effect SD", "Rank agreement", "Held-out-donor")]
    for i, (ax, sw, pr, col, xl, yl, ttl) in enumerate(spec):
        for mech in ["linear", "interaction", "ptm_dominant"]:
            g = mt[(mt.mechanism == mech) & (mt.sweep == sw) & (mt.protocol == pr)].sort_values("level")
            ax.plot(g.level, g[col], color=pf.MECH_C[mech], marker=pf.MECH_M[mech], ls=pf.MECH_LS[mech], ms=3.5, lw=1.1)
        pf.refline(ax, y=0.95 if col == "containment" else 0)
        ax.set_xlabel(xl), ax.set_ylabel(yl)
        if sw == "donor":
            ax.set_xticks([0.5, 1.0, 1.5, 2.0, 3.0])
        pf.panel(ax, "ABCD"[i], ttl, x=-0.32)
    mh = [Line2D([], [], color=pf.MECH_C[m], marker=pf.MECH_M[m], ls=pf.MECH_LS[m], ms=3.5, lw=1.1, label=pf.MECH_L[m] + " mechanism") for m in pf.MECH_C]
    fig.legend(handles=mh, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=3)
    made["fig6"] = pf.save(fig, OUT / "Fig6_outcome_mechanisms")


# ------------------------------------------------------------------ Fig 7 number of donors (raw sweep runs + primary grid at 40 donors)
def fig7():
    frames = []
    for n in [20, 80, 160]:
        d = pd.read_csv(SENS / f"nsize_{n}/scores.csv")
        frames.append(per_rep(d).assign(n=n))
    frames.append(P[P.scenario.isin(["clean_linear", "donor_only"])].assign(n=40))
    d = pd.concat(frames)
    ns = [20, 40, 80, 160]
    fig, axs = plt.subplots(1, 3, figsize=(pf.DOUBLE, 2.35), gridspec_kw={"wspace": 0.5})
    for ax, col, lab, ref, ttl in [(axs[0], "cont", "Containment of reference", 0.95, "Containment"), (axs[1], "bias", "Bias of Brier estimate", 0, "Bias"),
                                   (axs[2], "rho", "Rank agreement", None, "Rank agreement")]:
        for sc, ls, mfc in [("clean_linear", "-", None), ("donor_only", "--", "white")]:
            for p in PR:
                mu, ci = zip(*[mci(d[(d.n == n) & (d.scenario == sc) & (d.protocol == p)][col]) for n in ns])
                ax.errorbar(ns, mu, yerr=ci, color=pf.SPLIT_C[p], marker=pf.SPLIT_M[p], mfc=mfc or pf.SPLIT_C[p], ls=ls, ms=3.5, lw=1, elinewidth=0.7)
        if ref is not None:
            pf.refline(ax, y=ref)
        ax.set_xscale("log"), ax.set_xticks(ns, [str(n) for n in ns]), ax.minorticks_off()
        ax.set_xlabel("Development donors"), ax.set_ylabel(lab)
        pf.panel(ax, "ABC"["cont bias rho".split().index(col)], ttl, x=-0.3)
    hh = SH + [Line2D([], [], color=pf.INK2, ls="-", marker="o", ms=3.5, label="clean"), Line2D([], [], color=pf.INK2, ls="--", marker="o", mfc="white", ms=3.5, label="strong donor effect (SD 2.0)")]
    fig.legend(handles=hh, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=4)
    made["fig7"] = pf.save(fig, OUT / "Fig7_donor_number")


# ------------------------------------------------------------------ Fig S1 covariate-corruption sweep; Fig S2 per-model primary grid
def figS1():
    sweep_fig(pd.read_csv(SENS / "annotation_sweep/scores.csv"), "Covariate-corruption level", "FigS2_covariate_corruption_sweep")


def figS2():
    sm = raw.groupby(["scenario", "protocol", "model"]).agg(bias=("bias", "mean"), sd=("bias", "std"), n=("bias", "size"), cont=("covers_reference", "mean")).reset_index()
    sm["ci"] = 1.96 * sm.sd / np.sqrt(sm.n)
    fig, axs = plt.subplots(2, 4, figsize=(pf.DOUBLE, 4.6), sharey="row", gridspec_kw={"hspace": 0.85, "wspace": 0.12})
    x = np.arange(len(SC))
    for j, mo in enumerate(MODELS):
        for p in PR:
            g = sm[(sm.model == mo) & (sm.protocol == p)].set_index("scenario").reindex(SC)
            axs[0, j].errorbar(x + OFF[p], g.bias, yerr=g.ci, fmt=pf.SPLIT_M[p], color=pf.SPLIT_C[p], ms=3.2, lw=0.9)
            axs[1, j].plot(x + OFF[p], g.cont, pf.SPLIT_M[p], color=pf.SPLIT_C[p], ms=3.2)
        pf.refline(axs[0, j], y=0), pf.refline(axs[1, j], y=0.95)
        for r in (0, 1):
            axs[r, j].set_xticks(x, [SCL[s].replace("\n", " ") for s in SC], rotation=55, ha="right", fontsize=5.5)
        axs[0, j].set_title(pf.MODEL_L[mo].replace("Confidence-weighted", "Conf.-weighted"), fontsize=6.5, loc="center")
    axs[0, 0].set_ylabel("Bias of Brier estimate"), axs[1, 0].set_ylabel("Containment")
    pf.panel(axs[0, 0], "A", None, x=-0.45), pf.panel(axs[1, 0], "B", None, x=-0.45)
    fig.legend(handles=SH, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2)
    made["figS2"] = pf.save(fig, OUT / "FigS1_primary_grid_by_model")


for f in [fig1, fig2, fig3, fig4, fig5, fig6, fig7, figS1, figS2]:
    f()
for k, v in made.items():
    print(k, v.relative_to(R))
