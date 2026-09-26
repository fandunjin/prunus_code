#!/usr/bin/env python3
"""Expansion scan for the six-species inflorescence-family run.

For every subfamily group, compare per-track copy density in the racemose trio
(PZI/PSP/PPA; A and B tracks separately and combined) with the single-track
outgroups (PCAM, PPE, GTR), normalised by each species' gene total.

- enrichment_A/B: observed copies / expected copies if the trio track matched
  the outgroup gene density (1.0 = ancestral complement fully retained).
- enrichment_AB: observed / expected under the same outgroup density; a fully
  retained shared WGD gives ~2.0.  Values > 2.2 with p < 0.05 flag expansion
  beyond WGD retention.
"""
from __future__ import annotations

import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path

from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "analysis" / "inflorescence_family_pipeline_crossspecies6_20260912"
TRIO = ["pzi", "psp", "ppa"]
OUTS = ["pcam", "ppe", "gtr"]


def read_tsv(path: Path):
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> int:
    out = DEFAULT_OUT
    evidence = read_tsv(out / "candidates_crossspecies.tsv")
    bg = {(r["species"], r["track"]): int(r["n"]) for r in read_tsv(out / "stats" / "genome_background_genelevel.tsv")}

    genes_trio_a = sum(bg[("pzi", "A")] + bg[("psp", "A")] + bg[("ppa", "A")] for _ in [0])
    genes_trio_b = bg[("pzi", "B")] + bg[("psp", "B")] + bg[("ppa", "B")]
    genes_out = bg[("pcam", "S")] + bg[("ppe", "S")] + bg[("gtr", "S")]

    for tier in ("placed", "strict"):
        counts: dict[tuple[str, str], dict[str, int]] = defaultdict(lambda: defaultdict(int))
        for r in evidence:
            if not r["subfamily_group"]:
                continue
            if tier == "strict" and r["evidence_tier"] != "tree_supported":
                continue
            counts[(r["family"], r["subfamily_group"])][f"{r['species']}_{r['track']}"] += 1

        rows = []
        for (family, group), c in sorted(counts.items()):
            a_list = [c.get(f"{s}_A", 0) for s in TRIO]
            b_list = [c.get(f"{s}_B", 0) for s in TRIO]
            o_list = [c.get(f"{s}_S", 0) for s in OUTS]
            if sum(a_list) + sum(b_list) == 0:
                continue
            out_n = sum(o_list)
            p_out = out_n / genes_out if genes_out else 0.0
            k_a, k_b = sum(a_list), sum(b_list)
            enrich_a = (k_a / genes_trio_a) / p_out if p_out else float("nan")
            enrich_b = (k_b / genes_trio_b) / p_out if p_out else float("nan")
            enrich_ab = ((k_a + k_b) / (genes_trio_a + genes_trio_b)) / p_out if p_out else float("nan")
            p_a = binomtest(k_a, genes_trio_a, min(p_out, 1.0)).pvalue if p_out else float("nan")
            p_b = binomtest(k_b, genes_trio_b, min(p_out, 1.0)).pvalue if p_out else float("nan")
            p_ab = binomtest(k_a + k_b, genes_trio_a + genes_trio_b, min(p_out, 1.0)).pvalue if p_out else float("nan")
            flags = []
            if enrich_a > 1.25 and p_a < 0.05:
                flags.append("A_expanded")
            if enrich_b > 1.25 and p_b < 0.05:
                flags.append("B_expanded")
            if enrich_ab > 1.15 and p_ab < 0.05 and enrich_a > 1.1 and enrich_b > 1.1:
                flags.append("AB_expanded")
            if enrich_a < 0.75 and p_a < 0.05:
                flags.append("A_reduced")
            if enrich_b < 0.75 and p_b < 0.05:
                flags.append("B_reduced")
            if 0.8 <= enrich_a <= 1.25 and 0.8 <= enrich_b <= 1.25:
                flags.append("both_tracks_retained")
            rows.append({
                "family": family, "subfamily_group": group,
                "trio_A": k_a, "trio_B": k_b, "outgroup_S": out_n,
                "per_species": ";".join(f"{s}A={c.get(s + '_A', 0)}" for s in TRIO) + ";" +
                               ";".join(f"{s}B={c.get(s + '_B', 0)}" for s in TRIO) + ";" +
                               ";".join(f"{s}={c.get(s + '_S', 0)}" for s in OUTS),
                "enrichment_A": f"{enrich_a:.3f}", "p_A": f"{p_a:.4g}",
                "enrichment_B": f"{enrich_b:.3f}", "p_B": f"{p_b:.4g}",
                "enrichment_AB": f"{enrich_ab:.3f}", "p_AB": f"{p_ab:.4g}",
                "flags": ";".join(flags),
            })
        out_path = out / "stats" / f"expansion_scan_{tier}.tsv"
        if rows:
            with out_path.open("w", encoding="utf-8", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
                w.writeheader()
                w.writerows(rows)
        print(f"=== {tier} tier ({len(rows)} groups) ===")
        for r in sorted(rows, key=lambda x: -float(x["enrichment_AB"])):
            print(f"{r['subfamily_group']:14s} A={r['enrichment_A']:>6} (p={r['p_A']:<9}) "
                  f"B={r['enrichment_B']:>6} (p={r['p_B']:<9}) AB={r['enrichment_AB']:>6} "
                  f"(p={r['p_AB']:<9}) {r['flags']}")
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
