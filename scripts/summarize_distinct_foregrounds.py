#!/usr/bin/env python3
"""Collapse identical unrooted foreground assignments, selecting best fit per gene.
These assignments are a sensitivity analysis, not random permutations.
"""
import argparse, csv, re, json
from pathlib import Path
from collections import defaultdict
import numpy as np


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=Path, required=True)
    a = p.parse_args()
    rows = list(
        csv.DictReader((a.run / "alternative_foregrounds.tsv").open(), delimiter="\t")
    )
    groups = defaultdict(list)
    ref = next((a.run / "Background_A").iterdir())
    for i in range(16):
        tree = (ref / f"perm{i}/tree.nwk").read_text()
        signature = re.sub(r":[0-9.eE+-]+", "", tree).replace(" ", "").strip()
        groups[signature].append(i)
    chosen = []
    summ = []
    for no, (signature, configs) in enumerate(groups.items(), 1):
        subset = defaultdict(list)
        for r in rows:
            if int(r["Configuration"]) in configs:
                subset[r["Subgenome"], r["Gene_set"]].append(r)
        for (sg, gene), rs in subset.items():
            # Equivalent definitions may have different numerical local optima.
            best = max(rs, key=lambda r: float(r["LogL"]))
            chosen.append(
                dict(
                    Subgenome=sg,
                    Gene_set=gene,
                    Distinct_configuration=no,
                    Legacy_configurations=",".join(map(str, configs)),
                    Selected_legacy_configuration=best["Configuration"],
                    Omega_foreground=float(best["Omega_foreground"]),
                    Omega_background=float(best["Omega_background"]),
                    Omega_ratio=float(best["Omega_ratio"]),
                    LogL=float(best["LogL"]),
                )
            )
        for sg in ["A", "B"]:
            rs = [
                r["Omega_ratio"]
                for r in chosen
                if r["Distinct_configuration"] == no and r["Subgenome"] == sg
            ]
            summ.append(
                dict(
                    Subgenome=sg,
                    Distinct_configuration=no,
                    Legacy_configurations=",".join(map(str, configs)),
                    n=len(rs),
                    median_ratio=float(np.median(rs)),
                    mean_ratio=float(np.mean(rs)),
                    max_ratio=max(rs),
                )
            )
    for name, rs in [
        ("distinct_foreground_fits.tsv", chosen),
        ("distinct_foreground_summary.tsv", summ),
    ]:
        with (a.run / name).open("w") as f:
            w = csv.DictWriter(f, fieldnames=list(rs[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rs)
    (a.run / "foreground_equivalence.json").write_text(
        json.dumps(
            {
                "distinct_assignments": len(groups),
                "groups": list(groups.values()),
                "aggregation": "maximum likelihood fit within each gene and equivalent assignment; median and mean both retained",
                "interpretation": "predefined unrooted foreground-label sensitivity; not a permutation test",
            },
            indent=2,
        )
    )
    print("Distinct assignments:", len(groups))


if __name__ == "__main__":
    main()
