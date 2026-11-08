# Donor-grouped validation removes optimism but not the ranking errors caused by missed positives: a known-truth simulation of model evaluation in omics

Code, configurations, run manifests, result tables and figures for the manuscript *Donor-grouped validation removes optimism but not the ranking errors caused by missed positives: a known-truth simulation of model evaluation in omics* (submitted to BMC Bioinformatics).

## Repository layout
- `research_common/`, `topic*/`: analysis package (installable with `pip install -e .`)
- `modal_jobs/`: cloud job scripts (Modal), aggregation scripts, the manuscript builder and the independent audit
- `results/figures/`: every manuscript figure (PNG and LZW TIFF at 600 dpi, vector PDF), numbered as in the manuscript and Additional file 1; redrawn by `modal_jobs/figures_paper*.py`
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
- Fig1: figures_pub\Fig1_study_design.png
- Fig2: figures_pub\Fig2_primary_grid.png
- Fig3: figures_pub\Fig3_skill_vs_brier.png
- Fig4: figures_pub\Fig4_false_negative_sweep.png
- Fig5: figures_pub\Fig5_donor_effect_sweep.png
- Fig6: figures_pub\Fig6_outcome_mechanisms.png
- Fig7: figures_pub\Fig7_donor_number.png
- FigS1: figures_pub\FigS1_primary_grid_by_model.png
- FigS2: figures_pub\FigS2_covariate_corruption_sweep.png

## Tables
- main: Paper2_Table_1.csv, Paper2_Table_2.csv, Paper2_Table_3.csv, Paper2_Table_4.csv
- supplementary: Paper2_Table_S1.csv, Paper2_Table_S2.csv, Paper2_Table_S3.csv

## Citation
If you use this code, please cite the article (details to be added on publication) and this archive (see `CITATION.cff`).

## Licence
Code is released under the MIT licence (see `LICENSE`). Result tables and figures may be reused with attribution.
