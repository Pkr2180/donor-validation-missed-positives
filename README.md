# Donor-grouped validation removes optimism but not the ranking errors caused by missed positives: a known-truth simulation of model evaluation in omics

Code, configurations, run manifests, result tables and figures for the manuscript *Donor-grouped validation removes optimism but not the ranking errors caused by missed positives: a known-truth simulation of model evaluation in omics* (submitted to BMC Bioinformatics).

## Repository layout
- `research_common/`, `topic*/`: analysis package (installable with `pip install -e .`)
- `modal_jobs/`: cloud job scripts (Modal), aggregation scripts, the manuscript builder and the independent audit
- `results/figures/`: every manuscript figure (PNG and SVG), numbered as in the manuscript and Additional file 1
- `results/tables/main`, `results/tables/supplementary`: every manuscript table as CSV, numbered as in the manuscript
- `results/tables/source_data/`: the result tables the manuscript was built from
- `results/provenance/`: run manifests (package versions, SHA-256 hashes of inputs, code and outputs) and the audit report (PASS 136 | FAIL 0)

## Data
No experimental data are used. All data are simulated and are regenerated exactly from the configuration files and seeds (primary grid: `topic2_uncertainty/configs/research.json`, seed 29; sweeps: `modal_jobs/modal_sensitivity.py` and `modal_jobs/modal_extra.py`, seeds 29 and 31).

## Reproducing the analysis
Python 3.12 was used (see the run manifests for exact package versions).

```
pip install -e .
uncertainty-audit --help                     # primary simulation grid (configuration: topic2_uncertainty/configs/research.json)
python modal_jobs/build_topic2_v4.py         # rebuild the manuscript tables/figures from result files
python modal_jobs/audit_topic2.py            # independent audit (needs the raw per-replicate outputs)
```

## Figure mapping (manuscript number: original result file)
- Fig1: figures\fig1_bias_by_scenario.png
- Fig2: figures\fig2_containment.png
- Fig3: figures\fig3_meta_evaluation.png
- Fig4: sensitivity\figures\figS8_skill_vs_brier_bias.png
- Fig5: sensitivity\figures\figS4_false_negative_sweep.png
- Fig6: sensitivity\figures\figS5_donor_effect_sweep.png
- Fig7: sensitivity\figures\figS12_mechanism_fn_sweep.png
- Fig8: sensitivity\figures\figS7_donor_count_containment.png
- FigS1: sensitivity\figures\figS6_annotation_error_sweep.png
- FigS2: sensitivity\figures\figS13_mechanism_donor_sweep.png

## Tables
- main: Paper2_Table_1.csv, Paper2_Table_2.csv, Paper2_Table_3.csv, Paper2_Table_4.csv, Paper2_Table_5.csv
- supplementary: Paper2_Table_S1.csv, Paper2_Table_S2.csv, Paper2_Table_S3.csv

## Citation
If you use this code, please cite the article (details to be added on publication) and this archive (see `CITATION.cff`).

## Licence
Code is released under the MIT licence (see `LICENSE`). Result tables and figures may be reused with attribution.
