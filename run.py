#!/usr/bin/env python
import argparse
import warnings
from src.brainmaps.config import Config
from src.brainmaps import pipeline


def main():
    # Suppress RuntimeWarnings from neuromaps spin-test matrix multiplications
    warnings.filterwarnings(
        "ignore",
        message=".*encountered in matmul",
        module="neuromaps.nulls.spins",
        category=RuntimeWarning,
    )
    
    ap = argparse.ArgumentParser()
    ap.add_argument("targets", nargs="*", default=["all"])
    ap.add_argument("--config", default="configs/config.json")
    args = ap.parse_args()
    cfg = Config.from_json(args.config)

    steps = {
        "env":   lambda: pipeline.step_env(cfg),
        "transforms": lambda: pipeline.step_transforms(cfg),
        "viz":   lambda: pipeline.step_viz(cfg, pipeline.step_transforms(cfg)),
        "stats": lambda: pipeline.step_stats(cfg, pipeline.step_transforms(cfg)),
        "fdr":   lambda: pipeline.step_fdr(cfg, pipeline.step_stats(cfg, pipeline.step_transforms(cfg))),
        "results": lambda: pipeline.step_result(cfg),
        "all":   lambda: pipeline.run_all(cfg),
        "clean": lambda: print("Add rm -r out/ here if desired"),
    }
    for t in args.targets: steps[t]()

if __name__ == "__main__":
    main()
