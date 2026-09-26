#!/usr/bin/env python3
"""Test whether strict POInT anchor states form contiguous blocks.

The unit of analysis is an ordered GTR pillar.  A-only and B-only states are
treated as anchor-state labels, not validated biological gene-loss calls.
Permutations are performed independently within each GTR chromosome while
preserving the observed counts of A-only, B-only and all other states.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


SPECIES = ("PZI", "PSP", "PPA")
STATE_A = "A_only_anchor"
STATE_B = "B_only_anchor"
STATE_COL = {"PZI": "PZI_state", "PSP": "PSP_state", "PPA": "PPA_state"}
STATE_LABEL = {STATE_A: "A-only", STATE_B: "B-only"}
COLORS = {STATE_A: "#2b6cb0", STATE_B: "#c05621"}


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def chrom_key(seqid: str):
    digits = "".join(ch for ch in seqid if ch.isdigit())
    return (int(digits) if digits else 10**9, seqid)


def split_chromosomes(rows, state_col):
    groups = defaultdict(list)
    for row in rows:
        groups[row["gtr_seqid"]].append((int(row["gtr_global_index"]), row[state_col]))
    result = {}
    for chrom, values in groups.items():
        values.sort(key=lambda x: x[0])
        result[chrom] = {"index": np.array([x[0] for x in values], dtype=np.int64),
                         "state": np.array([x[1] for x in values], dtype=object)}
    return dict(sorted(result.items(), key=lambda item: chrom_key(item[0])))


def runs_for_mask(mask, indices):
    """Return run lengths, requiring adjacent GTR global indices."""
    if len(mask) == 0:
        return []
    starts = np.flatnonzero(mask & np.r_[True, ~mask[:-1] | (np.diff(indices) != 1)])
    if len(starts) == 0:
        return []
    ends = []
    for start in starts:
        end = int(start)
        while end + 1 < len(mask) and mask[end + 1] and indices[end + 1] == indices[end] + 1:
            end += 1
        ends.append(end)
    return [end - int(start) + 1 for start, end in zip(starts, ends)]


def metric_for_label(groups, label):
    run_lengths = []
    per_chrom = {}
    for chrom, data in groups.items():
        lengths = runs_for_mask(data["state"] == label, data["index"])
        run_lengths.extend(lengths)
        per_chrom[chrom] = lengths
    if not run_lengths:
        return {"n_state": 0, "n_runs": 0, "mean_run": math.nan, "median_run": math.nan,
                "max_run": 0.0, "run_lengths": [], "per_chrom": per_chrom}
    return {"n_state": int(sum(run_lengths)), "n_runs": len(run_lengths),
            "mean_run": float(np.mean(run_lengths)), "median_run": float(np.median(run_lengths)),
            "max_run": float(np.max(run_lengths)), "run_lengths": run_lengths,
            "per_chrom": per_chrom}


def random_category_metrics(groups, rng, n_a, n_b):
    """Shuffle three categories within chromosomes, keeping A/B counts."""
    a_lengths = []
    b_lengths = []
    for chrom, data in groups.items():
        n = len(data["state"])
        n_a_chrom = int(np.sum(data["state"] == STATE_A))
        n_b_chrom = int(np.sum(data["state"] == STATE_B))
        if n_a_chrom + n_b_chrom > n:
            raise ValueError("A/B counts exceed chromosome length")
        codes = np.zeros(n, dtype=np.int8)
        if n_a_chrom:
            codes[:n_a_chrom] = 1
        if n_b_chrom:
            codes[n_a_chrom:n_a_chrom + n_b_chrom] = 2
        rng.shuffle(codes)
        a_lengths.extend(runs_for_mask(codes == 1, data["index"]))
        b_lengths.extend(runs_for_mask(codes == 2, data["index"]))
    def summarize(values):
        return float(np.mean(values)) if values else math.nan, float(np.max(values)) if values else 0.0
    a_mean, a_max = summarize(a_lengths)
    b_mean, b_max = summarize(b_lengths)
    return a_mean, a_max, b_mean, b_max


def safe_z(value, null):
    sd = float(np.std(null, ddof=1)) if len(null) > 1 else 0.0
    return (value - float(np.mean(null))) / sd if sd > 0 else 0.0


def empirical_upper(value, null):
    return (1.0 + float(np.sum(null >= value))) / (len(null) + 1.0)


def empirical_two_sided(value, null):
    center = float(np.mean(null))
    return (1.0 + float(np.sum(np.abs(null - center) >= abs(value - center)))) / (len(null) + 1.0)


def permutation_ci(values):
    return float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))


def format_p(value):
    return "p<0.001" if value < 0.001 else f"p={value:.3f}"


def make_figure(summary_rows, chrom_rows, effect_rows, outdir: Path):
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8, "svg.fonttype": "none", "pdf.fonttype": 42,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.8, "legend.frameon": False,
    })
    fig = plt.figure(figsize=(10.8, 6.4), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05], width_ratios=[1.1, 1.0])
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, :])

    # Panel a: observed mean run length against chromosome-preserving null.
    xbase = np.arange(len(SPECIES))
    offsets = {STATE_A: -0.16, STATE_B: 0.16}
    for state in (STATE_A, STATE_B):
        subset = [r for r in summary_rows if r["state"] == state]
        observed = np.array([float(r["observed_mean_run"]) for r in subset])
        lo = np.array([float(r["null_mean_run_ci_low"]) for r in subset])
        hi = np.array([float(r["null_mean_run_ci_high"]) for r in subset])
        x = xbase + offsets[state]
        for xx, yy, ll, hh in zip(x, observed, lo, hi):
            ax_a.vlines(xx, ll, hh, color=COLORS[state], lw=2.0, alpha=0.85)
            ax_a.plot(xx, yy, "o", ms=5, color=COLORS[state], zorder=3)
        for xx, yy, row in zip(x, observed, subset):
            ax_a.text(xx, yy + 0.022, format_p(float(row["mean_run_empirical_p_upper"])),
                      ha="center", va="bottom", fontsize=6, color=COLORS[state])
        ax_a.plot([], [], "o", ms=5, color=COLORS[state], label=STATE_LABEL[state])
    ax_a.set_xticks(xbase, SPECIES)
    ax_a.set_ylabel("Mean contiguous run length\n(GTR pillars)")
    ax_a.set_title("Observed blockiness versus chromosome-preserving null", loc="left", fontsize=9, weight="bold")
    ax_a.axhline(1, color="#777777", lw=0.7, ls=":")
    ax_a.legend(ncol=2, loc="upper right")

    # Panel b: direct B-minus-A effect with a joint null interval.
    effects = [next(r for r in effect_rows if r["species"] == s) for s in SPECIES]
    obs = np.array([float(r["observed_B_minus_A_mean_run"]) for r in effects])
    lo = np.array([float(r["null_B_minus_A_mean_ci_low"]) for r in effects])
    hi = np.array([float(r["null_B_minus_A_mean_ci_high"]) for r in effects])
    y = np.arange(len(SPECIES))
    for yy, xx, ll, hh in zip(y, obs, lo, hi):
        ax_b.hlines(yy, ll, hh, color="#4a5568", lw=2.0, alpha=0.85)
        ax_b.plot(xx, yy, "o", ms=5, color="#4a5568", zorder=3)
    for yy, row in zip(y, effects):
        p = float(row["effect_empirical_p_two_sided"])
        ax_b.text(float(row["observed_B_minus_A_mean_run"]) - 0.015, yy - 0.12,
                  format_p(p), ha="right", va="bottom", fontsize=6, color="#4a5568")
    ax_b.axvline(0, color="#777777", lw=0.8, ls="--")
    ax_b.set_yticks(y, SPECIES)
    ax_b.set_xlabel("B-only minus A-only mean run length")
    ax_b.set_title("Direct A/B clustering contrast", loc="left", fontsize=9, weight="bold")
    ax_b.invert_yaxis()

    # Panel c: chromosome-level z scores, with observed state labels only.
    chroms = sorted({r["gtr_seqid"] for r in chrom_rows}, key=chrom_key)
    labels = []
    values = []
    for species in SPECIES:
        for state in (STATE_A, STATE_B):
            labels.append(f"{species} {STATE_LABEL[state]}")
            row = [r for r in chrom_rows if r["species"] == species and r["state"] == state]
            by_chrom = {r["gtr_seqid"]: float(r["z_mean_run"]) for r in row}
            values.append([by_chrom.get(c, 0.0) for c in chroms])
    im = ax_c.imshow(np.asarray(values), aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3)
    ax_c.set_xticks(np.arange(len(chroms)), [c.replace("GTR__", "") for c in chroms])
    ax_c.set_yticks(np.arange(len(labels)), labels)
    ax_c.set_xlabel("Ancestral GTR chromosome")
    ax_c.set_title("Chromosome-level clustering z scores", loc="left", fontsize=9, weight="bold")
    cbar = fig.colorbar(im, ax=ax_c, fraction=0.018, pad=0.012)
    cbar.set_label("z score versus within-chromosome permutation", rotation=90)

    for ax, letter in ((ax_a, "a"), (ax_b, "b"), (ax_c, "c")):
        ax.text(-0.13, 1.08, letter, transform=ax.transAxes, fontsize=12, weight="bold", va="top")
    for suffix, kwargs in (("png", {"dpi": 600}), ("tiff", {"dpi": 600}), ("pdf", {}), ("svg", {})):
        fig.savefig(outdir / f"figure5_anchor_state_clustering.{suffix}", bbox_inches="tight", **kwargs)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pillars", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--permutations", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260828)
    args = parser.parse_args()
    if args.permutations < 100:
        raise SystemExit("ERROR: use at least 100 permutations")
    rows = read_tsv(args.pillars)
    if not rows:
        raise SystemExit("ERROR: empty pillar table")
    required = {"gtr_seqid", "gtr_global_index", *STATE_COL.values()}
    missing = required - set(rows[0])
    if missing:
        raise SystemExit("ERROR: missing columns: " + ",".join(sorted(missing)))
    outdir = args.outdir
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit(f"ERROR: output directory is not empty: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    summary_rows = []
    chrom_rows = []
    effect_rows = []
    observed_runs_rows = []
    for species in SPECIES:
        groups = split_chromosomes(rows, STATE_COL[species])
        observed = {state: metric_for_label(groups, state) for state in (STATE_A, STATE_B)}
        null = {state: {"mean": [], "max": []} for state in (STATE_A, STATE_B)}
        effect_mean = []
        effect_max = []
        per_chrom_null = {state: defaultdict(lambda: {"mean": [], "max": []}) for state in (STATE_A, STATE_B)}
        for chrom, data in groups.items():
            for state in (STATE_A, STATE_B):
                obs_lengths = observed[state]["per_chrom"].get(chrom, [])
                for i, length in enumerate(obs_lengths, 1):
                    observed_runs_rows.append({"species": species, "state": STATE_LABEL[state],
                                               "gtr_seqid": chrom, "run_index": i, "run_length": length})
            for _ in range(args.permutations):
                # Joint shuffling preserves A/B counts and their mutual exclusivity.
                a_mean, a_max, b_mean, b_max = random_category_metrics({chrom: data}, rng, 0, 0)
                null[STATE_A]["mean"].append(a_mean)
                null[STATE_A]["max"].append(a_max)
                null[STATE_B]["mean"].append(b_mean)
                null[STATE_B]["max"].append(b_max)
                # Per-chromosome values are recomputed below using a compact one-chrom call.
                a_lengths = runs_for_mask(data["state"] == STATE_A, data["index"])
                b_lengths = runs_for_mask(data["state"] == STATE_B, data["index"])
                per_chrom_null[STATE_A][chrom]["mean"].append(float(np.mean(a_lengths)) if a_lengths else math.nan)
                per_chrom_null[STATE_A][chrom]["max"].append(float(np.max(a_lengths)) if a_lengths else 0.0)
                per_chrom_null[STATE_B][chrom]["mean"].append(float(np.mean(b_lengths)) if b_lengths else math.nan)
                per_chrom_null[STATE_B][chrom]["max"].append(float(np.max(b_lengths)) if b_lengths else 0.0)
        # The loop above stores one-chromosome null values; the global joint contrast
        # is generated independently so category counts are preserved per chromosome.
        null_global = {"A_mean": [], "B_mean": [], "A_max": [], "B_max": []}
        for _ in range(args.permutations):
            a_mean, a_max, b_mean, b_max = random_category_metrics(groups, rng, 0, 0)
            null_global["A_mean"].append(a_mean); null_global["A_max"].append(a_max)
            null_global["B_mean"].append(b_mean); null_global["B_max"].append(b_max)
        a_arr = np.asarray(null_global["A_mean"], dtype=float)
        b_arr = np.asarray(null_global["B_mean"], dtype=float)
        a_max_arr = np.asarray(null_global["A_max"], dtype=float)
        b_max_arr = np.asarray(null_global["B_max"], dtype=float)
        for state, key in ((STATE_A, "A"), (STATE_B, "B")):
            arr = a_arr if state == STATE_A else b_arr
            max_arr = a_max_arr if state == STATE_A else b_max_arr
            obs = observed[state]
            ci_lo, ci_hi = permutation_ci(arr)
            summary_rows.append({
                "species": species, "state": state, "state_label": STATE_LABEL[state],
                "n_state": obs["n_state"], "observed_n_runs": obs["n_runs"],
                "observed_mean_run": f"{obs['mean_run']:.6f}", "observed_median_run": f"{obs['median_run']:.6f}",
                "observed_max_run": f"{obs['max_run']:.6f}", "null_mean_run": f"{np.mean(arr):.6f}",
                "null_mean_run_ci_low": f"{ci_lo:.6f}", "null_mean_run_ci_high": f"{ci_hi:.6f}",
                "mean_run_z": f"{safe_z(obs['mean_run'], arr):.6f}",
                "mean_run_empirical_p_upper": f"{empirical_upper(obs['mean_run'], arr):.6f}",
                "max_run_null_ci_low": f"{permutation_ci(max_arr)[0]:.6f}",
                "max_run_null_ci_high": f"{permutation_ci(max_arr)[1]:.6f}",
                "max_run_empirical_p_upper": f"{empirical_upper(obs['max_run'], max_arr):.6f}",
                "permutations": args.permutations,
            })
        # Direct effect row and its joint null distribution.
        observed_diff = observed[STATE_B]["mean_run"] - observed[STATE_A]["mean_run"]
        null_diff = b_arr - a_arr
        effect_rows.append({
            "species": species,
            "observed_B_minus_A_mean_run": f"{observed_diff:.6f}",
            "null_B_minus_A_mean_run": f"{np.mean(null_diff):.6f}",
            "null_B_minus_A_mean_ci_low": f"{permutation_ci(null_diff)[0]:.6f}",
            "null_B_minus_A_mean_ci_high": f"{permutation_ci(null_diff)[1]:.6f}",
            "effect_z": f"{safe_z(observed_diff, null_diff):.6f}",
            "effect_empirical_p_two_sided": f"{empirical_two_sided(observed_diff, null_diff):.6f}",
            "observed_B_minus_A_max_run": f"{observed[STATE_B]['max_run'] - observed[STATE_A]['max_run']:.6f}",
            "permutations": args.permutations,
        })
        # Chromosome-level z scores use independent within-chromosome shuffles.
        for state in (STATE_A, STATE_B):
            for chrom, data in groups.items():
                values = []
                for _ in range(args.permutations):
                    n_a = int(np.sum(data["state"] == STATE_A)); n_b = int(np.sum(data["state"] == STATE_B))
                    codes = np.zeros(len(data["state"]), dtype=np.int8)
                    codes[:n_a] = 1; codes[n_a:n_a + n_b] = 2; rng.shuffle(codes)
                    lengths = runs_for_mask(codes == (1 if state == STATE_A else 2), data["index"])
                    values.append(float(np.mean(lengths)) if lengths else 0.0)
                obs_lengths = observed[state]["per_chrom"].get(chrom, [])
                obs_value = float(np.mean(obs_lengths)) if obs_lengths else 0.0
                chrom_rows.append({"species": species, "state": state, "state_label": STATE_LABEL[state],
                                   "gtr_seqid": chrom, "n_state": int(np.sum(data["state"] == state)),
                                   "observed_mean_run": f"{obs_value:.6f}",
                                   "null_mean_run": f"{np.mean(values):.6f}",
                                   "z_mean_run": f"{safe_z(obs_value, np.asarray(values)):.6f}",
                                   "empirical_p_upper": f"{empirical_upper(obs_value, np.asarray(values)):.6f}",
                                   "permutations": args.permutations})

    write_tsv(outdir / "anchor_state_clustering_summary.tsv", list(summary_rows[0]), summary_rows)
    write_tsv(outdir / "anchor_state_clustering_effect.tsv", list(effect_rows[0]), effect_rows)
    write_tsv(outdir / "anchor_state_clustering_by_chromosome.tsv", list(chrom_rows[0]), chrom_rows)
    write_tsv(outdir / "anchor_state_observed_runs.tsv", list(observed_runs_rows[0]), observed_runs_rows)
    make_figure(summary_rows, chrom_rows, effect_rows, outdir)
    (outdir / "method_contract.txt").write_text(
        "Question: are A-only and B-only strict POInT anchor states spatially clustered along GTR gene order?\n"
        "Unit: ordered GTR pillars; runs require adjacent gtr_global_index values.\n"
        "Null: independent within-chromosome shuffling preserving A-only, B-only and other-state counts.\n"
        "Interpretation: anchor-state blockiness, not validated biological gene loss.\n"
        f"Permutations: {args.permutations}; seed: {args.seed}.\n",
        encoding="ascii",
    )
    (outdir / "COMPLETE").write_text("Anchor-state clustering analysis completed.\n", encoding="ascii")
    print(f"[OK] Wrote anchor-state clustering analysis to {outdir}")


if __name__ == "__main__":
    main()
