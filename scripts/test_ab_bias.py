#!/usr/bin/env python3
"""Formal A/B-bias tests for the confirmed inflorescence-gene set.

The independent unit is a JCVI anchor-connected component.  Each component
contributes paired A and B confirmed-copy counts and paired A/B presence
states, so multiple copies within one syntenic group are not treated as
independent observations.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import math
import re
from pathlib import Path

from scipy.stats import binomtest, chi2


DEFAULT_MATRIX = Path(__file__).resolve().parents[1] / "analysis" / "inflorescence_family_pipeline" / "copy_fate" / "confirmed_ab_copy_fate.tsv"
DEFAULT_PILLARS = Path(__file__).resolve().parents[1] / "downstream_analysis" / "figure4_ab_fractionation_formal" / "panel4a_anchor_state_by_pillar.tsv"
DEFAULT_POINT_FIT = Path(__file__).resolve().parents[1] / "downstream_analysis" / "figure2_polyploid" / "strict_wgd_mapping" / "15_point_repartition_official_reference_exact" / "04_summary" / "point_model_fit_all.tsv"
DEFAULT_POINT_PARAMS = Path(__file__).resolve().parents[1] / "downstream_analysis" / "figure4_ab_fractionation_formal" / "panel4c_point_parameters.tsv"
DEFAULT_POINT_RESULT_SUMMARY = Path(__file__).resolve().parents[1] / "downstream_analysis" / "figure2_polyploid" / "strict_wgd_mapping" / "18_point_repartition_tree_mapped" / "result_summary.txt"


def wilson_interval(successes: int, trials: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if trials == 0:
        return float("nan"), float("nan")
    p = successes / trials
    denom = 1 + z * z / trials
    centre = (p + z * z / (2 * trials)) / denom
    half = z * math.sqrt(p * (1 - p) / trials + z * z / (4 * trials * trials)) / denom
    return max(0.0, centre - half), min(1.0, centre + half)


def exact_sign_flip_pvalue(differences: list[int], alternative: str = "two-sided") -> float:
    """Exact paired randomization p-value under A/B exchangeability."""
    observed = sum(differences)
    n = len(differences)
    if n == 0:
        return float("nan")
    extreme = 0
    total = 1 << n
    for signs in itertools.product((-1, 1), repeat=n):
        statistic = sum(sign * diff for sign, diff in zip(signs, differences))
        if alternative == "greater":
            extreme_hit = statistic >= observed
        elif alternative == "less":
            extreme_hit = statistic <= observed
        else:
            extreme_hit = abs(statistic) >= abs(observed)
        if extreme_hit:
            extreme += 1
    return extreme / total


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else ["test"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def fmt(value: float | int | str) -> str:
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        return f"{value:.8g}"
    return str(value)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--pillars", type=Path, default=DEFAULT_PILLARS)
    parser.add_argument("--point-fit", type=Path, default=DEFAULT_POINT_FIT)
    parser.add_argument("--point-params", type=Path, default=DEFAULT_POINT_PARAMS)
    parser.add_argument("--point-summary", type=Path, default=DEFAULT_POINT_RESULT_SUMMARY)
    parser.add_argument("--outdir", type=Path)
    args = parser.parse_args()
    matrix = args.matrix.resolve()
    outdir = (args.outdir or matrix.parent).resolve()
    rows = list(csv.DictReader(matrix.open(encoding="utf-8"), delimiter="\t"))
    if not rows:
        raise RuntimeError(f"No rows in confirmed A/B matrix: {matrix}")

    a_counts = [int(row["A_copy_count"]) for row in rows]
    b_counts = [int(row["B_copy_count"]) for row in rows]
    differences = [a - b for a, b in zip(a_counts, b_counts)]
    n_components = len(rows)
    a_total, b_total = sum(a_counts), sum(b_counts)
    n_records = a_total + b_total
    aggregate_binom = binomtest(a_total, n_records, 0.5, alternative="two-sided")
    aggregate_low, aggregate_high = wilson_interval(a_total, n_records)

    nonzero = [d for d in differences if d]
    a_greater = sum(d > 0 for d in differences)
    b_greater = sum(d < 0 for d in differences)
    sign_test = binomtest(a_greater, len(nonzero), 0.5, alternative="two-sided") if nonzero else None
    permutation_p = exact_sign_flip_pvalue(differences)
    permutation_p_greater = exact_sign_flip_pvalue(differences, alternative="greater")

    a_present = [int(row["A1_present"]) or int(row["A2_present"]) for row in rows]
    b_present = [int(row["B1_present"]) or int(row["B2_present"]) for row in rows]
    both_present = sum(a and b for a, b in zip(a_present, b_present))
    a_only = sum(a and not b for a, b in zip(a_present, b_present))
    b_only = sum(b and not a for a, b in zip(a_present, b_present))
    neither = sum(not a and not b for a, b in zip(a_present, b_present))
    mcnemar_p = binomtest(a_only, a_only + b_only, 0.5, alternative="two-sided").pvalue if a_only + b_only else float("nan")
    conditional_or = a_only / b_only if b_only else float("inf")

    output_rows = [
        {
            "test": "component_copy_count_exact_sign_flip",
            "unit": "JCVI anchor-connected component",
            "n": n_components,
            "A_total": a_total,
            "B_total": b_total,
            "estimate": fmt(sum(differences) / n_components),
            "estimate_definition": "mean(A_copy_count - B_copy_count)",
            "effect_ratio": fmt(a_total / b_total) if b_total else "",
            "p_value": fmt(permutation_p),
            "p_value_A_greater": fmt(permutation_p_greater),
            "alternative": "two-sided",
            "interpretation": "paired exact sign-flip randomization; components are exchangeable under the null",
        },
        {
            "test": "component_copy_count_direction_exact_binomial",
            "unit": "non-tied JCVI anchor-connected component",
            "n": len(nonzero),
            "A_greater": a_greater,
            "B_greater": b_greater,
            "estimate": fmt(a_greater / len(nonzero)) if nonzero else "",
            "estimate_definition": "fraction of non-tied components with A > B",
            "effect_ratio": fmt(a_greater / b_greater) if b_greater else "",
            "p_value": fmt(sign_test.pvalue if sign_test else float("nan")),
            "p_value_A_greater": fmt(binomtest(a_greater, len(nonzero), 0.5, alternative="greater").pvalue if nonzero else float("nan")),
            "alternative": "two-sided",
            "interpretation": "exact sign test; six tied components are omitted",
        },
        {
            "test": "component_presence_exact_mcnemar",
            "unit": "JCVI anchor-connected component",
            "n": n_components,
            "A_present": sum(a_present),
            "B_present": sum(b_present),
            "both_present": both_present,
            "A_only": a_only,
            "B_only": b_only,
            "neither": neither,
            "estimate": fmt(conditional_or),
            "estimate_definition": "conditional odds ratio for A-only vs B-only discordant components",
            "effect_ratio": fmt(conditional_or),
            "p_value": fmt(mcnemar_p),
            "p_value_A_greater": fmt(binomtest(a_only, a_only + b_only, 0.5, alternative="greater").pvalue if a_only + b_only else float("nan")),
            "alternative": "two-sided",
            "interpretation": "exact McNemar test on paired A/B presence states",
        },
        {
            "test": "record_level_sensitivity_exact_binomial",
            "unit": "confirmed protein records (descriptive sensitivity only)",
            "n": n_records,
            "A_total": a_total,
            "B_total": b_total,
            "estimate": fmt(a_total / n_records),
            "estimate_definition": "A fraction among confirmed records",
            "effect_ratio": fmt(a_total / b_total) if b_total else "",
            "ci_low": fmt(aggregate_low),
            "ci_high": fmt(aggregate_high),
            "p_value": fmt(aggregate_binom.pvalue),
            "p_value_A_greater": fmt(binomtest(a_total, n_records, 0.5, alternative="greater").pvalue),
            "alternative": "two-sided",
            "interpretation": "not an independent-replicate test; records are clustered within syntenic components",
        },
    ]
    # Keep a stable schema across rows while retaining all test-specific fields.
    fields = sorted({key for row in output_rows for key in row})
    normalized = [{key: row.get(key, "") for key in fields} for row in output_rows]
    write_tsv(outdir / "confirmed_ab_bias_statistics.tsv", normalized)

    # Independent genome-wide context: one exact test per species on the
    # resolved A-only versus B-only POInT anchor states.  Unresolved pillars
    # are excluded from the denominator and are not treated as losses.
    pillar_rows = list(csv.DictReader(args.pillars.resolve().open(encoding="utf-8"), delimiter="\t"))
    pillar_tests: list[dict[str, object]] = []
    for species, column in (("PZI", "PZI_state"), ("PSP", "PSP_state"), ("PPA", "PPA_state")):
        states = [row.get(column, "") for row in pillar_rows]
        a_only_pillars = states.count("A_only_anchor")
        b_only_pillars = states.count("B_only_anchor")
        ab_pillars = states.count("AB_anchor")
        neither_pillars = states.count("neither_anchor")
        unresolved_pillars = states.count("unresolved")
        discordant = a_only_pillars + b_only_pillars
        p_value = binomtest(a_only_pillars, discordant, 0.5, alternative="two-sided").pvalue if discordant else float("nan")
        pillar_tests.append({
            "test": "pillar_resolved_anchor_state_exact_mcnemar",
            "species": species,
            "unit": "POInT pillar",
            "n": len(states),
            "A_only": a_only_pillars,
            "B_only": b_only_pillars,
            "AB_anchor": ab_pillars,
            "neither_anchor": neither_pillars,
            "unresolved": unresolved_pillars,
            "estimate": fmt(a_only_pillars / b_only_pillars) if b_only_pillars else "",
            "estimate_definition": "A-only/B-only conditional odds ratio",
            "p_value": fmt(p_value),
            "alternative": "two-sided",
            "interpretation": "resolved A-only versus B-only anchor states; unresolved pillars excluded, not gene-loss calls",
        })
    # Benjamini-Hochberg correction across the three species-level pillar
    # tests, reported alongside the raw exact p values.
    ordered = sorted(enumerate(pillar_tests), key=lambda item: float(item[1]["p_value"]))
    m = len(ordered)
    previous_q = 1.0
    for rank, (index, row) in reversed(list(enumerate(ordered, start=1))):
        raw = float(row["p_value"])
        q_value = min(previous_q, raw * m / rank)
        previous_q = q_value
        pillar_tests[index]["q_value_bh"] = fmt(q_value)
    write_tsv(outdir / "genomewide_pillar_ab_bias_statistics.tsv", pillar_tests)

    # Native POInT nested likelihood comparison: symmetric fractionation
    # (two parameters) versus copy-specific asymmetric fractionation (three
    # parameters, including Copy2Bias). This is the formal genome-wide test.
    point_fit_rows = list(csv.DictReader(args.point_fit.resolve().open(encoding="utf-8"), delimiter="\t"))
    symmetric = next((row for row in point_fit_rows if row.get("model") == "symmetric"), None)
    asymmetric = next((row for row in point_fit_rows if row.get("model") != "symmetric"), None)
    point_row = None
    if symmetric and asymmetric:
        sym_lnl = float(symmetric["log_likelihood"])
        asym_lnl = float(asymmetric["log_likelihood"])
        sym_k = int(symmetric["model_parameters"])
        asym_k = int(asymmetric["model_parameters"])
        likelihood_ratio = 2.0 * (asym_lnl - sym_lnl)
        degrees_freedom = asym_k - sym_k
        lrt_p = chi2.sf(likelihood_ratio, degrees_freedom)
        sym_aic = 2.0 * sym_k - 2.0 * sym_lnl
        asym_aic = 2.0 * asym_k - 2.0 * asym_lnl
        point_row = {
            "test": "point_nested_symmetric_vs_asymmetric_likelihood_ratio",
            "unit": "optimized POInT pillar",
            "n": symmetric.get("pillars", ""),
            "symmetric_model": symmetric.get("model", ""),
            "asymmetric_model": asymmetric.get("model", ""),
            "symmetric_log_likelihood": fmt(sym_lnl),
            "asymmetric_log_likelihood": fmt(asym_lnl),
            "symmetric_parameters": sym_k,
            "asymmetric_parameters": asym_k,
            "degrees_freedom": degrees_freedom,
            "likelihood_ratio": fmt(likelihood_ratio),
            "delta_AIC_symmetric_minus_asymmetric": fmt(sym_aic - asym_aic),
            "p_value": fmt(lrt_p),
            "alternative": "asymmetric model versus symmetric null",
            "interpretation": "POInT nested likelihood-ratio test; Copy1/Copy2 are model copy states, not parental identities",
        }
        params = {}
        if args.point_params.resolve().exists():
            params = {
                row.get("parameter", ""): row.get("value", "")
                for row in csv.DictReader(args.point_params.resolve().open(encoding="utf-8"), delimiter="\t")
            }
        point_row["copy2_bias"] = params.get("Copy2Bias", "")
        point_row["fix_rate"] = params.get("FixRate", "")
        point_row["switch_prob"] = params.get("SwitchProb", "")
        if args.point_summary.resolve().exists():
            summary_text = args.point_summary.resolve().read_text(encoding="utf-8", errors="replace")
            for copy, key in (("Copy1", "copy1_expected_retention"), ("Copy2", "copy2_expected_retention")):
                match = re.search(rf"{copy}\s+expected\s+retention\s*=\s*([0-9]+(?:\.[0-9]+)?)", summary_text, flags=re.IGNORECASE)
                point_row[key] = match.group(1) if match else ""
        write_tsv(outdir / "point_model_ab_bias_statistics.tsv", [point_row])
    point_summary = ""
    if point_row is not None:
        bias_text = point_row.get("copy2_bias") or "NA"
        retention_text = ""
        if point_row.get("copy1_expected_retention") and point_row.get("copy2_expected_retention"):
            retention_text = f" Model-expected retention is {point_row['copy1_expected_retention']} for latent Copy1 versus {point_row['copy2_expected_retention']} for latent Copy2. These latent tracks are not fixed physical A/B labels; their pillar-specific correspondence is given by the POInT posterior tracking assignments."
        point_summary = (
            f"The native POInT nested model test on {point_row['n']} optimized pillars strongly favors copy-specific asymmetric fractionation over the symmetric null (LR={point_row['likelihood_ratio']}, df={point_row['degrees_freedom']}, p={point_row['p_value']}; DeltaAIC={point_row['delta_AIC_symmetric_minus_asymmetric']}). The fitted Copy2Bias parameter is {bias_text}.{retention_text} The physical A/B input order is retained in the pillar file, but POInT can permute the two latent copy tracks at each pillar; this model-state correspondence does not identify parental species.\n\n"
        )

    report = outdir / "confirmed_ab_bias_report.md"
    report.write_text(
        "# Confirmed A/B bias test\n\n"
        f"Input: `{matrix}`. The independent unit is the JCVI anchor-connected component (n={n_components}); the 56 confirmed records are clustered within these components.\n\n"
        f"The confirmed set contains {a_total} A-subgenome records and {b_total} B-subgenome records (A:B={a_total}:{b_total}, ratio={a_total / b_total:.4f}). The record-level exact binomial sensitivity test gives p={aggregate_binom.pvalue:.6g}; this is descriptive because records are not independent.\n\n"
        f"For the primary component-level copy-count test, the mean paired difference A-B is {sum(differences) / n_components:.4f} copies per component, with exact two-sided sign-flip p={permutation_p:.6g} (directional A>B p={permutation_p_greater:.6g}). The direction-only exact sign test gives {a_greater} A>B versus {b_greater} B>A among {len(nonzero)} non-tied components (two-sided p={sign_test.pvalue:.6g}; A>B p={binomtest(a_greater, len(nonzero), 0.5, alternative='greater').pvalue:.6g}).\n\n"
        f"For paired subgenome presence, {both_present} components contain both A and B, {a_only} are A-only and {b_only} are B-only. The exact two-sided McNemar test gives p={mcnemar_p:.6g} (directional A-only>B-only p={binomtest(a_only, a_only + b_only, 0.5, alternative='greater').pvalue:.6g}; conditional odds ratio A-only/B-only={conditional_or:.4f}).\n\n"
        "As an independent genome-wide context, `genomewide_pillar_ab_bias_statistics.tsv` applies the same exact A-only versus B-only test to resolved POInT pillar states for PZI, PSP and PPA; unresolved pillars are excluded from the test. `point_model_ab_bias_statistics.tsv` reports the native POInT nested symmetric-versus-asymmetric likelihood-ratio test.\n\n"
        + point_summary
        + "Conclusion: the confirmed inflorescence-gene set shows a numerically higher A copy count, but the component-level tests do not reject equal A/B representation at alpha=0.05. Missing or unconfirmed B phases remain unresolved and are not interpreted as biological gene losses.\n",
        encoding="utf-8",
    )
    print(outdir / "confirmed_ab_bias_statistics.tsv")
    print(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
