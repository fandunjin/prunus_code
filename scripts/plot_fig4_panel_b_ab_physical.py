#!/usr/bin/env python3
"""Plot Fig. 4b from the resolved physical A/B subgenomes.

The input is restricted to the 0903 GTR-centred pillar table.  No POInT
Copy1/Copy2 state is read or relabelled.  The right-hand panel reports
parsimony-compatible presence/absence patterns under the topology
PPA,(PSP,PZI); patterns requiring more than one loss are kept as unresolved
multi-loss/older-absence classes.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7.0,
        "svg.fonttype": "none",
        "axes.linewidth": 0.7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
    }
)


ROOT = Path(__file__).resolve().parents[1]
PILLARS = (
    ROOT
    / "0903"
    / "downstream_analysis"
    / "figure2_polyploid"
    / "strict_wgd_mapping"
    / "14_point_repartition_reference_exact"
    / "01_pillars"
    / "gtr_centered_pillars.tsv"
)
FIGURE_DIR = ROOT / "0903" / "figures"
OUT_BASE = FIGURE_DIR / "figure4_panel_b_ab_physical"

SPECIES = ["PZI", "PSP", "PPA"]
TRACKS = ["PZI-A", "PZI-B", "PSP-A", "PSP-B", "PPA-A", "PPA-B"]
TRACK_COLORS = {"A": "#C96561", "B": "#4E78A7"}

PATTERN_LABELS = {
    "111": "All three retained",
    "011": "PPA-branch attrition",
    "101": "PSP-branch attrition",
    "110": "PZI-branch attrition",
    "100": "PSP+PZI shared-branch attrition",
    "other": "Multiple-loss or older absence",
}
PATTERN_COLORS = {
    "111": "#365C73",
    "011": "#D8913A",
    "101": "#7B6BA8",
    "110": "#4F8A6E",
    "100": "#B56B8C",
    "other": "#B9BDC2",
}


def wilson(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    den = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / den
    half = z * np.sqrt((p * (1.0 - p) + z * z / (4.0 * n)) / n) / den
    return max(0.0, center - half), min(1.0, center + half)


def normalize_states(data: pd.DataFrame) -> pd.DataFrame:
    required = [
        "pillar_id",
        "gtr_global_index",
        "gtr_seqid",
        "gtr_gene",
        *[f"{sp}_{track}_state" for sp in SPECIES for track in ("A", "B")],
    ]
    missing = [col for col in required if col not in data.columns]
    if missing:
        raise ValueError(f"Missing required pillar columns: {missing}")
    data = data.sort_values("gtr_global_index", kind="stable").reset_index(drop=True)
    if len(data) != 13204:
        raise ValueError(f"Expected 13,204 GTR-centred pillars, found {len(data):,}")
    for sp in SPECIES:
        for track in ("A", "B"):
            col = f"{sp}_{track}"
            data[col] = data[f"{col}_state"].fillna(0).astype(int).clip(0, 1)
    return data


def pattern_class(row: pd.Series, track: str) -> str:
    # The tuple order follows the dated species topology PPA,(PSP,PZI).
    pattern = "".join(str(int(row[f"{sp}_{track}"])) for sp in ("PPA", "PSP", "PZI"))
    return pattern if pattern in PATTERN_LABELS else "other"


def write_source_tables(data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    source_columns = [
        "pillar_id",
        "gtr_global_index",
        "gtr_seqid",
        "gtr_gene",
        *[f"{sp}_{track}" for sp in SPECIES for track in ("A", "B")],
    ]
    data[source_columns].to_csv(
        FIGURE_DIR / "figure4_panel_b_ab_physical_source.tsv",
        sep="\t",
        index=False,
    )

    retention_rows: list[dict[str, object]] = []
    n = len(data)
    for sp in SPECIES:
        for track in ("A", "B"):
            k = int(data[f"{sp}_{track}"].sum())
            low, high = wilson(k, n)
            retention_rows.append(
                {
                    "species": sp,
                    "subgenome": track,
                    "n_pillars": n,
                    "retained_count": k,
                    "retained_fraction": k / n,
                    "retained_percent": 100.0 * k / n,
                    "wilson95_low": low,
                    "wilson95_high": high,
                    "wilson95_low_percent": 100.0 * low,
                    "wilson95_high_percent": 100.0 * high,
                    "input_semantics": "physical A/B pillar presence; not POInT Copy1/Copy2",
                }
            )
    retention = pd.DataFrame(retention_rows)
    retention.to_csv(
        FIGURE_DIR / "figure4_panel_b_ab_retention_summary.tsv",
        sep="\t",
        index=False,
    )

    pattern_rows: list[dict[str, object]] = []
    for track in ("A", "B"):
        labels = data.apply(lambda row: pattern_class(row, track), axis=1)
        for pattern in PATTERN_LABELS:
            k = int((labels == pattern).sum())
            low, high = wilson(k, n)
            pattern_rows.append(
                {
                    "subgenome": track,
                    "presence_pattern": pattern,
                    "pattern_label": PATTERN_LABELS[pattern],
                    "n_pillars": n,
                    "count": k,
                    "fraction": k / n,
                    "percent": 100.0 * k / n,
                    "wilson95_low_percent": 100.0 * low,
                    "wilson95_high_percent": 100.0 * high,
                    "pattern_order": "PPA,PSP,PZI",
                    "input_semantics": "physical A/B pillar presence; parsimony-compatible annotation",
                }
            )
    patterns = pd.DataFrame(pattern_rows)
    patterns.to_csv(
        FIGURE_DIR / "figure4_panel_b_ab_branch_attrition.tsv",
        sep="\t",
        index=False,
    )
    return retention, patterns


def draw_figure(retention: pd.DataFrame, patterns: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(8.25, 4.85), facecolor="white")
    gs = fig.add_gridspec(
        1,
        2,
        width_ratios=(0.95, 1.55),
        left=0.085,
        right=0.985,
        bottom=0.22,
        top=0.78,
        wspace=0.38,
    )
    ax_ret = fig.add_subplot(gs[0, 0])
    ax_pat = fig.add_subplot(gs[0, 1])

    # Left: exact retention proportions with Wilson intervals.
    x = np.arange(len(SPECIES))
    width = 0.28
    for offset, track in [(-width / 2.0, "A"), (width / 2.0, "B")]:
        subset = retention[retention["subgenome"] == track].set_index("species").loc[SPECIES]
        y = subset["retained_percent"].to_numpy()
        low = y - subset["wilson95_low_percent"].to_numpy()
        high = subset["wilson95_high_percent"].to_numpy() - y
        ax_ret.bar(
            x + offset,
            y,
            width=width,
            color=TRACK_COLORS[track],
            edgecolor="white",
            linewidth=0.45,
            label=f"{track} subgenome",
            zorder=3,
        )
        ax_ret.errorbar(
            x + offset,
            y,
            yerr=np.vstack([low, high]),
            fmt="none",
            ecolor="#252525",
            elinewidth=0.65,
            capsize=2.0,
            capthick=0.65,
            zorder=4,
        )
        for xi, yi in zip(x + offset, y):
            ax_ret.text(xi, yi + 0.55, f"{yi:.1f}", ha="center", va="bottom", fontsize=5.8)
    ax_ret.set_ylim(80, 100.8)
    ax_ret.set_xticks(x, SPECIES)
    ax_ret.set_ylabel("GTR-anchored pillars retained (%)")
    ax_ret.set_title("Resolved A/B retention", loc="left", fontsize=8.2, pad=8)
    ax_ret.grid(axis="y", color="#E1E4E6", linewidth=0.5, zorder=0)
    ax_ret.set_axisbelow(True)
    ax_ret.text(0.70, 1.055, "A", transform=ax_ret.transAxes, color=TRACK_COLORS["A"], fontsize=6.5, fontweight="bold", va="bottom")
    ax_ret.text(0.78, 1.055, "B", transform=ax_ret.transAxes, color=TRACK_COLORS["B"], fontsize=6.5, fontweight="bold", va="bottom")
    ax_ret.text(0.0, -0.19, "Bars: 13,204 GTR-centred pillars; error bars: Wilson 95% CI", transform=ax_ret.transAxes, fontsize=5.7, color="#58616A")

    # Right: physical A/B presence patterns, with no Copy nomenclature.  The
    # patterns are defined across the three-species topology and therefore
    # appear once per resolved subgenome, rather than being repeated by species.
    order = ["111", "011", "101", "110", "100", "other"]
    track_labels = ["A", "B"]
    track_keys = ["A", "B"]
    y_positions = np.arange(len(track_labels))
    left = np.zeros(len(track_labels))
    for pattern in order:
        vals = []
        for track in track_keys:
            vals.append(float(patterns.loc[patterns["subgenome"] == track].set_index("presence_pattern").loc[pattern, "fraction"]))
        vals = np.asarray(vals)
        ax_pat.barh(y_positions, vals, left=left, height=0.60, color=PATTERN_COLORS[pattern], edgecolor="white", linewidth=0.35, label=pattern)
        for yi, x0, value in zip(y_positions, left, vals):
            if value >= 0.065:
                ax_pat.text(x0 + value / 2.0, yi, f"{100 * value:.1f}", ha="center", va="center", fontsize=5.2, color="white" if pattern in {"111", "101", "110", "100"} else "#303030", fontweight="bold")
        left += vals
    ax_pat.set_xlim(0, 1)
    ax_pat.set_xticks([0, 0.25, 0.5, 0.75, 1.0], ["0", "25", "50", "75", "100"])
    ax_pat.set_xlabel("Pillar presence pattern (%)")
    ax_pat.set_yticks(y_positions, track_labels)
    ax_pat.invert_yaxis()
    ax_pat.tick_params(axis="y", length=0, pad=4)
    ax_pat.grid(axis="x", color="#E1E4E6", linewidth=0.5, zorder=0)
    ax_pat.set_axisbelow(True)
    ax_pat.set_title("Branch-associated attrition under PPA,(PSP,PZI)", loc="left", fontsize=8.2, pad=8)
    ax_pat.text(0.0, -0.19, "Pattern order: PPA, PSP, PZI; 'other' requires multiple losses or predates the sampled split", transform=ax_pat.transAxes, fontsize=5.7, color="#58616A")

    handles = [Patch(facecolor=PATTERN_COLORS[p], label=PATTERN_LABELS[p]) for p in order]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.57, 0.965), ncol=3, fontsize=5.8, handlelength=1.0, columnspacing=1.0)
    fig.text(0.010, 0.985, "b", ha="left", va="top", fontsize=10, fontweight="bold")
    fig.text(0.085, 0.875, "Resolved physical A/B chromosome assignments; POInT Copy1/Copy2 states are not used.", fontsize=6.1, color="#58616A")

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{OUT_BASE}.svg", bbox_inches="tight")
    fig.savefig(f"{OUT_BASE}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def write_boundary() -> None:
    candidate_like = []
    for path in (ROOT / "0903" / "karyotype").rglob("*"):
        if path.is_file() and any(token in path.name.lower() for token in ("breakpoint", "event", "rearrangement")):
            candidate_like.append(str(path.relative_to(ROOT / "0903")))
    boundary = (
        "Fig. 4b uses only the 0903 GTR-centred pillar table and its explicit physical A/B state columns.\n"
        "POInT Copy1/Copy2 outputs are not read, renamed or interpreted as A/B.\n"
        "Fig. 4b retention is observed physical A/B pillar presence; branch classes are parsimony-compatible presence/absence patterns under PPA,(PSP,PZI), not independently validated gene-loss calls.\n"
        "Fig. 4c and Fig. 4d were not generated because no explicit RCT/EEJ/inversion breakpoint-candidate table was found under 0903/karyotype.\n"
        "WGDI collinearity/blockinfo files are evidence inputs, but were not auto-labelled as breakpoint mechanisms.\n"
        f"Candidate-like filenames found under 0903/karyotype: {len(candidate_like)}\n"
    )
    (FIGURE_DIR / "figure4_panel_c_d_status.txt").write_text(boundary, encoding="ascii")


def main() -> None:
    if not PILLARS.is_file() or PILLARS.stat().st_size == 0:
        raise SystemExit(f"Missing or empty input: {PILLARS}")
    data = normalize_states(pd.read_csv(PILLARS, sep="\t"))
    retention, patterns = write_source_tables(data)
    draw_figure(retention, patterns)
    write_boundary()
    print(f"[OK] Fig. 4b written to {OUT_BASE}.svg and {OUT_BASE}.png")
    print("[OK] Fig. 4c/d withheld: no explicit breakpoint candidate table under 0903/karyotype")


if __name__ == "__main__":
    main()
