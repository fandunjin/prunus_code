#!/usr/bin/env python3
"""Family-stratified A/B tests for confirmed inflorescence regulators.

The independent unit is a JCVI anchor-connected component.  The script also
compares an expanded confirmed matrix with the original baseline matrix using
the same exact tests, so added families are not mistaken for genome-wide
fractionation evidence.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import math
from collections import defaultdict
from pathlib import Path

from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MATRIX = ROOT / "analysis" / "inflorescence_family_pipeline" / "copy_fate" / "confirmed_ab_copy_fate.tsv"


def exact_sign_flip(differences: list[int], alternative: str = "two-sided") -> float:
    if not differences:
        return float("nan")
    observed = sum(differences)
    # Dynamic convolution of the signed-difference distribution.  This is
    # mathematically identical to enumerating all 2^n sign assignments but
    # remains fast for the expanded 20-30 component matrices.
    distribution: dict[int, int] = {0: 1}
    for difference in differences:
        next_distribution: dict[int, int] = defaultdict(int)
        for subtotal, count in distribution.items():
            next_distribution[subtotal + difference] += count
            next_distribution[subtotal - difference] += count
        distribution = dict(next_distribution)
    if alternative == "greater":
        extreme = sum(count for subtotal, count in distribution.items() if subtotal >= observed)
    elif alternative == "less":
        extreme = sum(count for subtotal, count in distribution.items() if subtotal <= observed)
    else:
        extreme = sum(count for subtotal, count in distribution.items() if abs(subtotal) >= abs(observed))
    return extreme / (2 ** len(differences))


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fmt(value: float | int | str) -> str:
    if isinstance(value, float):
        return "" if math.isnan(value) else f"{value:.8g}"
    return str(value)


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row})
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows({key: row.get(key, "") for key in fields} for row in rows)


def family_rows(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        families = [x for x in row.get("tree_family", "").split(";") if x]
        for family in families:
            grouped[family].append(row)

    output: list[dict[str, object]] = []
    for family, records in sorted(grouped.items()):
        differences = [int(r["A_copy_count"]) - int(r["B_copy_count"]) for r in records]
        nonzero = [d for d in differences if d]
        a_greater = sum(d > 0 for d in nonzero)
        b_greater = sum(d < 0 for d in nonzero)
        a_present = [int(r["A1_present"]) or int(r["A2_present"]) for r in records]
        b_present = [int(r["B1_present"]) or int(r["B2_present"]) for r in records]
        a_only = sum(a and not b for a, b in zip(a_present, b_present))
        b_only = sum(b and not a for a, b in zip(a_present, b_present))
        discordant = a_only + b_only
        output.append({
            "family": family,
            "n_components": len(records),
            "n_non_tied_components": len(nonzero),
            "A_total": sum(int(r["A_copy_count"]) for r in records),
            "B_total": sum(int(r["B_copy_count"]) for r in records),
            "A_B_ratio": fmt(sum(int(r["A_copy_count"]) for r in records) / sum(int(r["B_copy_count"]) for r in records)) if sum(int(r["B_copy_count"]) for r in records) else "",
            "mean_A_minus_B": fmt(sum(differences) / len(records)) if records else "",
            "A_greater": a_greater,
            "B_greater": b_greater,
            "p_sign_flip_two_sided": fmt(exact_sign_flip(differences)),
            "p_sign_flip_A_greater": fmt(exact_sign_flip(differences, "greater")),
            "p_sign_test_two_sided": fmt(binomtest(a_greater, len(nonzero), 0.5).pvalue) if nonzero else "",
            "p_sign_test_A_greater": fmt(binomtest(a_greater, len(nonzero), 0.5, alternative="greater").pvalue) if nonzero else "",
            "A_only_components": a_only,
            "B_only_components": b_only,
            "both_present_components": sum(a and b for a, b in zip(a_present, b_present)),
            "p_mcnemar_two_sided": fmt(binomtest(a_only, discordant, 0.5).pvalue) if discordant else "",
            "p_mcnemar_A_greater": fmt(binomtest(a_only, discordant, 0.5, alternative="greater").pvalue) if discordant else "",
            "interpretation": "component-level exact tests; missing/unconfirmed phases are unresolved",
        })
    return output


def bh(rows: list[dict[str, object]], column: str, output_column: str) -> None:
    ordered = sorted(
        ((i, float(row[column])) for i, row in enumerate(rows) if row.get(column) not in (None, "")),
        key=lambda x: x[1],
    )
    previous = 1.0
    m = len(ordered)
    for rank, (index, raw) in reversed(list(enumerate(ordered, start=1))):
        q = min(previous, raw * m / rank)
        previous = q
        rows[index][output_column] = fmt(q)


def aggregate(rows: list[dict[str, str]], label: str) -> dict[str, object]:
    differences = [int(r["A_copy_count"]) - int(r["B_copy_count"]) for r in rows]
    nonzero = [d for d in differences if d]
    a = sum(int(r["A_copy_count"]) for r in rows)
    b = sum(int(r["B_copy_count"]) for r in rows)
    a_present = [int(r["A1_present"]) or int(r["A2_present"]) for r in rows]
    b_present = [int(r["B1_present"]) or int(r["B2_present"]) for r in rows]
    a_only = sum(x and not y for x, y in zip(a_present, b_present))
    b_only = sum(y and not x for x, y in zip(a_present, b_present))
    discordant = a_only + b_only
    return {
        "matrix": label,
        "n_components": len(rows),
        "A_total": a,
        "B_total": b,
        "A_B_ratio": fmt(a / b) if b else "",
        "p_sign_flip_two_sided": fmt(exact_sign_flip(differences)),
        "p_sign_flip_A_greater": fmt(exact_sign_flip(differences, "greater")),
        "p_mcnemar_two_sided": fmt(binomtest(a_only, discordant).pvalue) if discordant else "",
        "p_mcnemar_A_greater": fmt(binomtest(a_only, discordant, alternative="greater").pvalue) if discordant else "",
        "A_only_components": a_only,
        "B_only_components": b_only,
        "confirmed_records": a + b,
        "added_records_vs_baseline": "",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    expanded = read_rows(args.matrix.resolve())
    baseline = read_rows(args.baseline.resolve()) if args.baseline.resolve().exists() else []
    stats = family_rows(expanded)
    bh(stats, "p_sign_flip_two_sided", "q_sign_flip_bh")
    bh(stats, "p_mcnemar_two_sided", "q_mcnemar_bh")
    write_tsv(args.outdir.resolve() / "confirmed_family_ab_bias_statistics.tsv", stats)

    aggregates = [aggregate(baseline, "baseline_original")]
    aggregates.append(aggregate(expanded, "expanded_confirmed"))
    if baseline:
        aggregates[1]["added_records_vs_baseline"] = (
            int(aggregates[1]["confirmed_records"]) - int(aggregates[0]["confirmed_records"])
        )
    write_tsv(args.outdir.resolve() / "confirmed_family_ab_bias_comparison.tsv", aggregates)

    print(args.outdir.resolve() / "confirmed_family_ab_bias_statistics.tsv")
    print(args.outdir.resolve() / "confirmed_family_ab_bias_comparison.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
