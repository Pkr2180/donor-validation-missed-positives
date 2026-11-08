# Topic 2 results (Modal run r1 + sensitivity analyses)

Final manuscript: BMC_Bioinformatics_Topic2_v4.docx (older drafts in archive_drafts/).
- Primary grid (200 replicates x 6 scenarios): summary.csv, scores.csv, meta_evaluation.csv, paired_comparisons.csv, tables/, figures/
- Sensitivity: sensitivity/tables, sensitivity/figures (S4-S8), sensitivity/runs/<sweep>/
- Large raw files (unit scores, split assignments) remain in outputs/modal_t2 and outputs/sens_t2.
- Regenerate: python modal_jobs/make_report.py; python modal_jobs/make_sensitivity_report.py 2; python modal_jobs/build_topic2_v4.py (run analyze_extra.py t2 first) (run from modal_jobs/)
- Extra analyses (Modal account yadalam1980): mechanism sweeps and pooled 200-replicate linear sweeps; see modal_jobs/analyze_extra.py.
