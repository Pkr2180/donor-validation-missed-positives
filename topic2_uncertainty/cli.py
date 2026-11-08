import argparse
from pathlib import Path

from research_common.io import managed_run, read_config, read_table, write_json, write_table
from research_common.plots import audit_scatter

from .benchmark import run_benchmark, validate_config
from .pseudobulk import aggregate_counts
from .target_ranking import rank_targets


def main():
    parser = argparse.ArgumentParser(
        description="Simulation-based meta-evaluation, with optional empirical covariates"
    )
    subs = parser.add_subparsers(dest="command", required=True)
    for command in ["demo", "plasmode"]:
        p = subs.add_parser(command)
        p.add_argument("--config", default=str(Path(__file__).parent / "configs" / "demo.json"))
        p.add_argument("--output", required=True)
        if command == "plasmode":
            p.add_argument("--covariates", required=True)
            p.add_argument("--provenance", required=True)
    p = subs.add_parser("rank-targets")
    p.add_argument("--nodes", required=True)
    p.add_argument("--edges", required=True)
    p.add_argument("--draws", type=int, default=1000)
    p.add_argument("--seed", type=int, default=19)
    p.add_argument("--top-k", type=int)
    p.add_argument("--output", required=True)
    p = subs.add_parser("pseudobulk")
    p.add_argument("--counts", required=True)
    p.add_argument("--annotations", required=True)
    p.add_argument("--minimum-confidence", type=float, default=0.0)
    p.add_argument("--output", required=True)
    args = parser.parse_args()
    if args.command == "rank-targets":
        with managed_run(
            args.output, vars(args), [args.nodes, args.edges], "EXPLORATORY_TARGET_RANKING"
        ) as out:
            result = rank_targets(
                read_table(args.nodes), read_table(args.edges), args.draws, args.seed, args.top_k
            )
            write_table(out / "target_ranking.csv", result)
        return
    if args.command == "pseudobulk":
        with managed_run(
            args.output, vars(args), [args.counts, args.annotations], "PSEUDOBULK_PREPROCESSING"
        ) as out:
            counts, qc = aggregate_counts(
                read_table(args.counts), read_table(args.annotations), args.minimum_confidence
            )
            write_table(out / "pseudobulk_counts.csv", counts)
            write_table(out / "annotation_qc.csv", qc)
        return
    config = read_config(args.config)
    validate_config(config)
    inputs = [args.config]
    source = None
    provenance = {"data_origin": "FULLY SYNTHETIC DEMONSTRATION"}
    if args.command == "plasmode":
        from research_common.provenance import validate_provenance

        provenance = read_config(args.provenance)
        validate_provenance(provenance, workflow="plasmode")
        source = read_table(args.covariates)
        inputs.extend([args.covariates, args.provenance])
    status = "SYNTHETIC_DEMONSTRATION" if source is None else "PLASMODE_IMPOSED_TRUTH"
    with managed_run(args.output, config, inputs, status, provenance) as out:
        results = run_benchmark(config, source, progress=lambda text: print(text, flush=True))
        for name, table in results.items():
            write_table(out / (name + ".csv"), table)
        audit_scatter(
            results["scores"],
            "estimated_brier",
            "reference_expected_brier",
            out / "simulation_audit.png",
            status + "\nEstimated versus imposed reference loss",
        )
        write_json(
            out / "interpretation.json",
            {
                "status": status,
                "estimand": "Equal-simulated-donor expected Brier loss on an independent finite simulated reference population",
                "interval": "95% percentile bootstrap, conditional on fitted models; random-row train/test share donors",
                "containment": "Observed-label interval containment of reference risk; not unconditional nominal coverage",
                "weighted_baseline": "Uses assumed annotation reliability and label-dependent evidence weights; not a calibrated PU estimator",
                "plasmode_limit": "Internal donor holdout and bootstrap use source donors. Reference draws can reuse source donors: conditional simulation, not independent biological validation.",
                "scope": "Imposed signals are not experimentally validated periodontal targets or PTM interactions.",
            },
        )
    print(f"Completed {status}: {out.resolve()}")


if __name__ == "__main__":
    main()
