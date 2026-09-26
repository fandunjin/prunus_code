#!/usr/bin/env python3
"""Block-level Ks summaries, clustered bootstrap, and mixture stability.

The input is the all-pairs YN result produced by ks_02.  Each JCVI block is
reduced to one median Ks before inference so large blocks do not dominate the
comparison.  Gaussian mixtures are fitted on log(Ks) and are treated as a
descriptive stability analysis, not as automatic WGD calls.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


FOCAL = [
    "pzi_hapA1_vs_pzi_hapB1",
    "psp_hapA1_vs_psp_hapB1",
    "ppa_hapA_vs_ppa_hapB",
]

LABELS = {
    "pzi_hapA1_vs_pzi_hapB1": "P. zippeliana A1-B1",
    "psp_hapA1_vs_psp_hapB1": "P. spinulosa A1-B1",
    "ppa_hapA_vs_ppa_hapB": "P. padus A-B",
}

COLORS = {
    "pzi_hapA1_vs_pzi_hapB1": "#B9473F",
    "psp_hapA1_vs_psp_hapB1": "#185A9D",
    "ppa_hapA_vs_ppa_hapB": "#3D929B",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--comparisons", nargs="+", default=FOCAL)
    parser.add_argument("--min-ks", type=float, default=0.0)
    parser.add_argument("--max-ks", type=float, default=3.0)
    parser.add_argument("--min-pairs-per-block", type=int, default=5)
    parser.add_argument("--min-pairs-sensitivity", nargs="+", type=int, default=[3, 5, 10])
    parser.add_argument("--bootstrap-replicates", type=int, default=1000)
    parser.add_argument("--mixture-bootstrap-replicates", type=int, default=1000)
    parser.add_argument("--max-components", type=int, default=3)
    parser.add_argument("--mixture-initializations", type=int, default=30)
    parser.add_argument("--bootstrap-initializations", type=int, default=8)
    parser.add_argument("--min-component-blocks", type=int, default=5)
    parser.add_argument("--min-component-fraction", type=float, default=0.10)
    parser.add_argument("--kde-bandwidth", type=float, default=0.02)
    parser.add_argument("--seed", type=int, default=20260815)
    return parser.parse_args()


def quantile(values: Sequence[float], probability: float) -> float:
    if not values:
        return float("nan")
    return float(np.quantile(np.asarray(values, dtype=float), probability))


def fmt(value: float) -> str:
    return "NA" if not math.isfinite(value) else "{:.8f}".format(value)


def write_tsv(path: Path, rows: Iterable[dict], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_blocks(
    path: Path,
    comparisons: Sequence[str],
    min_ks: float,
    max_ks: float,
) -> Dict[str, Dict[str, dict]]:
    selected = set(comparisons)
    blocks: Dict[str, Dict[str, dict]] = defaultdict(dict)
    required = {"comparison", "block_id", "Ks", "query_seqid", "ref_seqid"}
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise SystemExit("Missing required input columns: {}".format(", ".join(missing)))
        for row in reader:
            comparison = row["comparison"]
            if comparison not in selected:
                continue
            block_id = row.get("block_id", "").strip()
            if not block_id:
                continue
            record = blocks[comparison].setdefault(
                block_id,
                {
                    "values": [],
                    "n_rows": 0,
                    "query_seqids": set(),
                    "ref_seqids": set(),
                    "comparison_type": row.get("comparison_type", ""),
                },
            )
            record["n_rows"] += 1
            if row.get("query_seqid"):
                record["query_seqids"].add(row["query_seqid"])
            if row.get("ref_seqid"):
                record["ref_seqids"].add(row["ref_seqid"])
            try:
                value = float(row.get("Ks", ""))
            except ValueError:
                continue
            if math.isfinite(value) and min_ks < value < max_ks:
                record["values"].append(value)
    return blocks


def summarize_blocks(blocks: Dict[str, Dict[str, dict]]) -> List[dict]:
    rows: List[dict] = []
    for comparison in sorted(blocks):
        for block_id, record in sorted(blocks[comparison].items()):
            values = sorted(record["values"])
            if not values:
                continue
            rows.append(
                {
                    "comparison": comparison,
                    "display_label": LABELS.get(comparison, comparison),
                    "comparison_type": record["comparison_type"],
                    "block_id": block_id,
                    "n_pairs_total": record["n_rows"],
                    "n_valid_ks": len(values),
                    "ks_median": float(np.median(values)),
                    "ks_q25": quantile(values, 0.25),
                    "ks_q75": quantile(values, 0.75),
                    "ks_min": values[0],
                    "ks_max": values[-1],
                    "query_seqids": ",".join(sorted(record["query_seqids"])),
                    "ref_seqids": ",".join(sorted(record["ref_seqids"])),
                }
            )
    return rows


def retained_values(rows: Sequence[dict], comparisons: Sequence[str], min_pairs: int) -> Dict[str, List[float]]:
    grouped: Dict[str, List[float]] = {comparison: [] for comparison in comparisons}
    for row in rows:
        if row["comparison"] in grouped and int(row["n_valid_ks"]) >= min_pairs:
            grouped[row["comparison"]].append(float(row["ks_median"]))
    return grouped


def bootstrap_medians(values: Sequence[float], replicates: int, rng: np.random.Generator) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.size == 0:
        return np.asarray([], dtype=float)
    indexes = rng.integers(0, array.size, size=(replicates, array.size))
    return np.median(array[indexes], axis=1)


def normal_logpdf(x: np.ndarray, mean: float, sd: float) -> np.ndarray:
    return -0.5 * math.log(2.0 * math.pi) - math.log(sd) - 0.5 * ((x - mean) / sd) ** 2


def fit_gmm_logks(
    values: Sequence[float],
    components: int,
    initializations: int,
    rng: np.random.Generator,
    min_component_blocks: int,
    min_component_fraction: float,
) -> dict:
    x = np.log(np.asarray(values, dtype=float))
    n = x.size
    if n < max(8, components * 5):
        return {"valid": False, "reason": "too_few_blocks"}
    global_sd = max(float(np.std(x, ddof=1)), 0.05)
    best = None
    for init in range(initializations):
        if init == 0:
            probabilities = np.linspace(0.15, 0.85, components)
            means = np.quantile(x, probabilities)
        else:
            means = np.sort(rng.choice(x, size=components, replace=False))
            means = means + rng.normal(0.0, global_sd * 0.05, size=components)
        sds = np.repeat(global_sd, components)
        weights = np.repeat(1.0 / components, components)
        previous = -float("inf")
        valid = True
        for _ in range(600):
            log_terms = np.column_stack(
                [np.log(max(weights[k], 1e-12)) + normal_logpdf(x, means[k], sds[k]) for k in range(components)]
            )
            row_max = np.max(log_terms, axis=1)
            log_sum = row_max + np.log(np.sum(np.exp(log_terms - row_max[:, None]), axis=1))
            loglik = float(np.sum(log_sum))
            responsibilities = np.exp(log_terms - log_sum[:, None])
            nk = np.sum(responsibilities, axis=0)
            if np.any(nk < 2.0):
                valid = False
                break
            weights = nk / n
            means = np.sum(responsibilities * x[:, None], axis=0) / nk
            variances = np.sum(responsibilities * (x[:, None] - means) ** 2, axis=0) / nk
            sds = np.sqrt(np.maximum(variances, 0.0025))
            order = np.argsort(means)
            weights, means, sds = weights[order], means[order], sds[order]
            if abs(loglik - previous) < 1e-8 * (1.0 + abs(loglik)):
                break
            previous = loglik
        if not valid:
            continue
        parameter_count = 3 * components - 1
        candidate = {
            "valid": True,
            "log_likelihood": loglik,
            "aic": 2.0 * parameter_count - 2.0 * loglik,
            "bic": math.log(n) * parameter_count - 2.0 * loglik,
            "weights": weights,
            "log_means": means,
            "log_sds": sds,
            "effective_counts": weights * n,
        }
        if best is None or candidate["log_likelihood"] > best["log_likelihood"]:
            best = candidate
    if best is None:
        return {"valid": False, "eligible": False, "reason": "all_initializations_failed"}
    minimum_effective = max(float(min_component_blocks), min_component_fraction * n)
    best["eligible"] = bool(np.all(best["effective_counts"] >= minimum_effective))
    best["reason"] = "" if best["eligible"] else "component_below_minimum_effective_blocks"
    return best


def mixture_row(comparison: str, components: int, model: dict, n_blocks: int) -> dict:
    row = {
        "comparison": comparison,
        "display_label": LABELS.get(comparison, comparison),
        "n_blocks": n_blocks,
        "components": components,
        "valid": "yes" if model.get("valid") else "no",
        "eligible": "yes" if model.get("eligible") else "no",
        "log_likelihood": "NA",
        "aic": "NA",
        "bic": "NA",
        "component_centers_ks": "",
        "component_weights": "",
        "component_log_sd": "",
        "component_effective_blocks": "",
        "reason": model.get("reason", ""),
    }
    if model.get("valid"):
        row.update(
            {
                "log_likelihood": fmt(float(model["log_likelihood"])),
                "aic": fmt(float(model["aic"])),
                "bic": fmt(float(model["bic"])),
                "component_centers_ks": ",".join("{:.8f}".format(math.exp(x)) for x in model["log_means"]),
                "component_weights": ",".join("{:.8f}".format(x) for x in model["weights"]),
                "component_log_sd": ",".join("{:.8f}".format(x) for x in model["log_sds"]),
                "component_effective_blocks": ",".join("{:.4f}".format(x) for x in model["effective_counts"]),
            }
        )
    return row


def reflected_kde(values: Sequence[float], grid: np.ndarray, bandwidth: float) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    if x.size == 0:
        return np.zeros_like(grid)
    z1 = (grid[:, None] - x[None, :]) / bandwidth
    z2 = (grid[:, None] + x[None, :]) / bandwidth
    return np.sum(np.exp(-0.5 * z1 * z1) + np.exp(-0.5 * z2 * z2), axis=1) / (
        x.size * bandwidth * math.sqrt(2.0 * math.pi)
    )


def save_figure(
    outdir: Path,
    grouped: Dict[str, List[float]],
    bootstrap_summary: Dict[str, Tuple[float, float, float]],
    mixture_models: Dict[str, Dict[int, dict]],
    comparisons: Sequence[str],
    bandwidth: float,
) -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
            "font.size": 7,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.7))
    ax_density, ax_ecdf, ax_ci, ax_bic = axes.ravel()
    display_max = 0.5
    grid = np.linspace(0.0, display_max, 800)
    for comparison in comparisons:
        values = grouped.get(comparison, [])
        if not values:
            continue
        color = COLORS.get(comparison, "#555555")
        label = "{}  n={}".format(LABELS.get(comparison, comparison), len(values))
        density = reflected_kde(values, grid, bandwidth)
        ax_density.plot(grid, density, color=color, lw=1.5, label=label)
        ax_density.fill_between(grid, density, color=color, alpha=0.08)
        ordered = np.sort(np.asarray(values))
        ecdf = np.arange(1, ordered.size + 1) / ordered.size
        ax_ecdf.step(ordered, ecdf, where="post", color=color, lw=1.4, label=label)
    ax_density.set(xlim=(0, display_max), xlabel=r"Block median $K_s$", ylabel="Probability density")
    ax_density.legend(loc="upper right", fontsize=6.2)
    ax_ecdf.set(xlim=(0, display_max), ylim=(0, 1.01), xlabel=r"Block median $K_s$", ylabel="Cumulative probability")

    positions = np.arange(len(comparisons))[::-1]
    for position, comparison in zip(positions, comparisons):
        estimate, low, high = bootstrap_summary[comparison]
        color = COLORS.get(comparison, "#555555")
        ax_ci.errorbar(
            estimate,
            position,
            xerr=[[estimate - low], [high - estimate]],
            fmt="o",
            ms=4,
            color=color,
            ecolor=color,
            elinewidth=1.2,
            capsize=2,
        )
    ax_ci.set_yticks(positions)
    ax_ci.set_yticklabels([LABELS.get(x, x) for x in comparisons])
    ax_ci.set_xlabel(r"Median of block medians $K_s$ (95% CI)")
    ax_ci.grid(axis="x", color="#E5E5E5", lw=0.6)

    for comparison in comparisons:
        models = mixture_models.get(comparison, {})
        fitted = [(k, model) for k, model in sorted(models.items()) if model.get("valid")]
        if not fitted:
            continue
        best_bic = min(float(model["bic"]) for _, model in fitted)
        color = COLORS.get(comparison, "#555555")
        x_values = [k for k, _ in fitted]
        y_values = [float(model["bic"]) - best_bic for _, model in fitted]
        ax_bic.plot(x_values, y_values, lw=1.0, color=color, alpha=0.75)
        for k, model in fitted:
            ax_bic.plot(
                k,
                float(model["bic"]) - best_bic,
                marker="o",
                ms=4,
                markerfacecolor=color if model.get("eligible") else "white",
                markeredgecolor=color,
                markeredgewidth=1.0,
            )
    ax_bic.set_xticks(range(1, max(max(mixture_models[x]) for x in mixture_models) + 1))
    ax_bic.set(
        xlabel="Gaussian mixture components",
        ylabel=r"Raw $\Delta$BIC (open: occupancy filter failed)",
    )
    ax_bic.axhline(0, color="#777777", lw=0.7, ls=":")

    for label, axis in zip("abcd", axes.ravel()):
        axis.text(-0.13, 1.05, label, transform=axis.transAxes, fontsize=9, fontweight="bold", va="top")
    fig.tight_layout(w_pad=1.4, h_pad=1.5)
    stem = outdir / "blockwise_median_ks"
    fig.savefig(str(stem) + ".svg", bbox_inches="tight")
    fig.savefig(str(stem) + ".pdf", bbox_inches="tight")
    fig.savefig(str(stem) + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(str(stem) + ".tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    if args.min_pairs_per_block < 1:
        raise SystemExit("--min-pairs-per-block must be positive")
    if args.max_components < 1 or args.max_components > 3:
        raise SystemExit("--max-components must be between 1 and 3")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    blocks = read_blocks(Path(args.input), args.comparisons, args.min_ks, args.max_ks)
    rows = summarize_blocks(blocks)
    if not rows:
        raise SystemExit("No valid block-level Ks values were found")
    block_fields = [
        "comparison", "display_label", "comparison_type", "block_id", "n_pairs_total", "n_valid_ks",
        "ks_median", "ks_q25", "ks_q75", "ks_min", "ks_max", "query_seqids", "ref_seqids",
    ]
    write_tsv(outdir / "block_medians.all.tsv", rows, block_fields)

    grouped = retained_values(rows, args.comparisons, args.min_pairs_per_block)
    retained_rows = [row for row in rows if int(row["n_valid_ks"]) >= args.min_pairs_per_block]
    write_tsv(outdir / "block_medians.retained.tsv", retained_rows, block_fields)

    summary_rows = []
    bootstrap_by_comparison: Dict[str, np.ndarray] = {}
    bootstrap_summary: Dict[str, Tuple[float, float, float]] = {}
    for comparison in args.comparisons:
        values = grouped.get(comparison, [])
        if len(values) < 10:
            raise SystemExit("Fewer than 10 retained blocks for {}".format(comparison))
        boot = bootstrap_medians(values, args.bootstrap_replicates, rng)
        bootstrap_by_comparison[comparison] = boot
        estimate = float(np.median(values))
        low, high = float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))
        bootstrap_summary[comparison] = (estimate, low, high)
        summary_rows.append(
            {
                "comparison": comparison,
                "display_label": LABELS.get(comparison, comparison),
                "n_blocks_retained": len(values),
                "n_valid_gene_pairs": sum(int(row["n_valid_ks"]) for row in retained_rows if row["comparison"] == comparison),
                "median_of_block_medians": fmt(estimate),
                "ci025": fmt(low),
                "ci975": fmt(high),
                "mean_of_block_medians": fmt(float(np.mean(values))),
                "bootstrap_replicates": args.bootstrap_replicates,
                "bootstrap_unit": "JCVI_block_id",
                "min_pairs_per_block": args.min_pairs_per_block,
                "ks_filter": "{}<Ks<{}".format(args.min_ks, args.max_ks),
            }
        )
    write_tsv(
        outdir / "blockwise_bootstrap_summary.tsv",
        summary_rows,
        [
            "comparison", "display_label", "n_blocks_retained", "n_valid_gene_pairs",
            "median_of_block_medians", "ci025", "ci975", "mean_of_block_medians",
            "bootstrap_replicates", "bootstrap_unit", "min_pairs_per_block", "ks_filter",
        ],
    )

    contrast_rows = []
    for left_index, left in enumerate(args.comparisons):
        for right in args.comparisons[left_index + 1 :]:
            differences = bootstrap_by_comparison[left] - bootstrap_by_comparison[right]
            observed = bootstrap_summary[left][0] - bootstrap_summary[right][0]
            contrast_rows.append(
                {
                    "comparison_left": left,
                    "comparison_right": right,
                    "display_contrast": "{} minus {}".format(LABELS.get(left, left), LABELS.get(right, right)),
                    "observed_difference": fmt(observed),
                    "bootstrap_median_difference": fmt(float(np.median(differences))),
                    "ci025": fmt(float(np.quantile(differences, 0.025))),
                    "ci975": fmt(float(np.quantile(differences, 0.975))),
                    "two_sided_p": fmt(
                        min(
                            1.0,
                            2.0
                            * (
                                min(int(np.sum(differences <= 0)), int(np.sum(differences >= 0))) + 1
                            )
                            / (len(differences) + 1),
                        )
                    ),
                    "bootstrap_replicates": args.bootstrap_replicates,
                }
            )
    write_tsv(
        outdir / "blockwise_bootstrap_contrasts.tsv",
        contrast_rows,
        [
            "comparison_left", "comparison_right", "display_contrast", "observed_difference",
            "bootstrap_median_difference", "ci025", "ci975", "two_sided_p", "bootstrap_replicates",
        ],
    )

    sensitivity_rows = []
    for minimum in sorted(set(args.min_pairs_sensitivity)):
        sensitivity = retained_values(rows, args.comparisons, minimum)
        for comparison in args.comparisons:
            values = sensitivity.get(comparison, [])
            sensitivity_rows.append(
                {
                    "comparison": comparison,
                    "display_label": LABELS.get(comparison, comparison),
                    "min_pairs_per_block": minimum,
                    "n_blocks": len(values),
                    "median_of_block_medians": fmt(float(np.median(values))) if values else "NA",
                    "q25": fmt(quantile(values, 0.25)),
                    "q75": fmt(quantile(values, 0.75)),
                }
            )
    write_tsv(
        outdir / "minimum_block_size_sensitivity.tsv",
        sensitivity_rows,
        ["comparison", "display_label", "min_pairs_per_block", "n_blocks", "median_of_block_medians", "q25", "q75"],
    )

    mixture_models: Dict[str, Dict[int, dict]] = defaultdict(dict)
    mixture_rows = []
    for comparison in args.comparisons:
        values = grouped[comparison]
        for components in range(1, args.max_components + 1):
            model = fit_gmm_logks(
                values,
                components,
                args.mixture_initializations,
                rng,
                args.min_component_blocks,
                args.min_component_fraction,
            )
            mixture_models[comparison][components] = model
            mixture_rows.append(mixture_row(comparison, components, model, len(values)))
    write_tsv(
        outdir / "mixture_model_selection.tsv",
        mixture_rows,
        [
            "comparison", "display_label", "n_blocks", "components", "valid", "eligible", "log_likelihood",
            "aic", "bic", "component_centers_ks", "component_weights", "component_log_sd",
            "component_effective_blocks", "reason",
        ],
    )

    mixture_boot_rows = []
    selection_counts: Dict[Tuple[str, int], int] = defaultdict(int)
    for comparison in args.comparisons:
        values = np.asarray(grouped[comparison], dtype=float)
        for replicate in range(1, args.mixture_bootstrap_replicates + 1):
            sample = values[rng.integers(0, len(values), size=len(values))]
            fitted = {}
            for components in range(1, args.max_components + 1):
                fitted[components] = fit_gmm_logks(
                    sample,
                    components,
                    args.bootstrap_initializations,
                    rng,
                    args.min_component_blocks,
                    args.min_component_fraction,
                )
            valid = [(k, model) for k, model in fitted.items() if model.get("valid") and model.get("eligible")]
            if not valid:
                continue
            selected_k, selected = min(valid, key=lambda item: float(item[1]["bic"]))
            selection_counts[(comparison, selected_k)] += 1
            mixture_boot_rows.append(
                {
                    "comparison": comparison,
                    "replicate": replicate,
                    "selected_components_bic": selected_k,
                    "selected_bic": fmt(float(selected["bic"])),
                    "component_centers_ks": ",".join("{:.8f}".format(math.exp(x)) for x in selected["log_means"]),
                    "component_weights": ",".join("{:.8f}".format(x) for x in selected["weights"]),
                }
            )
    write_tsv(
        outdir / "mixture_bootstrap_replicates.tsv",
        mixture_boot_rows,
        ["comparison", "replicate", "selected_components_bic", "selected_bic", "component_centers_ks", "component_weights"],
    )
    stability_rows = []
    for comparison in args.comparisons:
        total = sum(selection_counts[(comparison, k)] for k in range(1, args.max_components + 1))
        for components in range(1, args.max_components + 1):
            count = selection_counts[(comparison, components)]
            stability_rows.append(
                {
                    "comparison": comparison,
                    "display_label": LABELS.get(comparison, comparison),
                    "components": components,
                    "selected_count": count,
                    "valid_bootstrap_replicates": total,
                    "selection_fraction_bic": fmt(count / total) if total else "NA",
                }
            )
    write_tsv(
        outdir / "mixture_bootstrap_stability.tsv",
        stability_rows,
        ["comparison", "display_label", "components", "selected_count", "valid_bootstrap_replicates", "selection_fraction_bic"],
    )

    save_figure(outdir, grouped, bootstrap_summary, mixture_models, args.comparisons, args.kde_bandwidth)
    (outdir / "analysis_parameters.txt").write_text(
        "\n".join(
            [
                "input={}".format(Path(args.input).resolve()),
                "ks_filter={}<Ks<{}".format(args.min_ks, args.max_ks),
                "min_pairs_per_block={}".format(args.min_pairs_per_block),
                "min_pairs_sensitivity={}".format(",".join(map(str, sorted(set(args.min_pairs_sensitivity))))),
                "bootstrap_replicates={}".format(args.bootstrap_replicates),
                "mixture_bootstrap_replicates={}".format(args.mixture_bootstrap_replicates),
                "mixture_scale=natural_log_Ks",
                "mixture_selection=BIC_primary_AIC_reported",
                "max_components={}".format(args.max_components),
                "min_component_blocks={}".format(args.min_component_blocks),
                "min_component_fraction={}".format(args.min_component_fraction),
                "seed={}".format(args.seed),
                "interpretation=descriptive_components_not_automatic_WGD_calls",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (outdir / "COMPLETE").write_text("blockwise_median_ks_complete\n", encoding="utf-8")
    print("[OK] Block-wise median Ks analysis completed: {}".format(outdir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
