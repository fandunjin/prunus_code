#!/usr/bin/env python3
"""Visualize the current PSP-to-repartitioned-PZI chromosome assignment.

Inputs are restricted to the current PSP assignment table and the PAF
generated against the current C1A/C1B repartitioned PZI reference. The script
does not read any legacy haploid_rename or analysis/figure2 outputs.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


TARGET_ORDER = [
    *(f"C{i}A_A1" for i in range(1, 9)),
    *(f"C{i}A_A2" for i in range(1, 9)),
    *(f"C{i}B_B1" for i in range(1, 9)),
    *(f"C{i}B_B2" for i in range(1, 9)),
]
PHASE_COLORS = {
    "A1": "#0072B2",
    "A2": "#56B4E9",
    "B1": "#D55E00",
    "B2": "#E69F00",
}
TARGET_RE = re.compile(r"^C[1-8][AB]_([AB][12])$")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise SystemExit(f"ERROR: TSV has no header: {path}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    name = None
    length = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if raw.startswith(">"):
                if name is not None:
                    lengths[name] = length
                name = raw[1:].split(None, 1)[0]
                length = 0
            else:
                length += len(raw.strip())
    if name is not None:
        lengths[name] = length
    return lengths


def parse_paf(path: Path, min_mapq: int, min_aln_len: int):
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 12:
                continue
            qname, qlen, qstart, qend = fields[0], int(fields[1]), int(fields[2]), int(fields[3])
            strand = fields[4]
            tname, tlen, tstart, tend = fields[5], int(fields[6]), int(fields[7]), int(fields[8])
            mapq = int(fields[11])
            aln_len = qend - qstart
            if aln_len < min_aln_len or mapq < min_mapq:
                continue
            tp = "?"
            for item in fields[12:]:
                if item.startswith("tp:A:"):
                    tp = item[-1]
                    break
            yield {
                "qname": qname,
                "qlen": qlen,
                "qstart": qstart,
                "qend": qend,
                "strand": strand,
                "tname": tname,
                "tlen": tlen,
                "tstart": tstart,
                "tend": tend,
                "mapq": mapq,
                "aln_len": aln_len,
                "tp": tp,
            }


def union_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    total = 0
    start, end = sorted(intervals)[0]
    for left, right in sorted(intervals)[1:]:
        if left <= end:
            end = max(end, right)
        else:
            total += end - start
            start, end = left, right
    return total + end - start


def phase_of(target: str) -> str:
    match = TARGET_RE.fullmatch(target)
    if not match:
        raise ValueError(f"Unexpected PZI target name: {target}")
    return match.group(1)


def cumulative_offsets(names: list[str], lengths: dict[str, int], gap: int = 2_000_000):
    offsets = {}
    boundaries = []
    cursor = 0
    for name in names:
        length = lengths.get(name, 0)
        if length <= 0:
            raise SystemExit(f"Missing length for sequence {name}")
        offsets[name] = cursor
        cursor += length
        boundaries.append((cursor, name))
        cursor += gap
    return offsets, boundaries, cursor - gap


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 7,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "axes.linewidth": 0.7,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def save_figure(fig, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{stem}.svg", bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paf", type=Path, required=True)
    parser.add_argument("--assignment", type=Path, required=True)
    parser.add_argument("--pzi-reference", type=Path, required=True)
    parser.add_argument("--psp-chr32", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument("--min-aln-len", type=int, default=10000)
    parser.add_argument("--max-points", type=int, default=250000)
    args = parser.parse_args()

    assignment = read_tsv(args.assignment)
    required = {"old_seqid", "new_seqid", "hap_set", "pzi_target"}
    if len(assignment) != 32 or not required.issubset(assignment[0]):
        raise SystemExit(f"ERROR: expected 32-row current PSP assignment table: {args.assignment}")
    assignment_by_query = {row["old_seqid"]: row for row in assignment}
    if len(assignment_by_query) != 32:
        raise SystemExit("ERROR: PSP assignment contains duplicate old_seqid values")
    if set(row["new_seqid"] for row in assignment) != set(TARGET_ORDER):
        raise SystemExit("ERROR: PSP assignment does not contain exactly the current C1A/C1B phase targets")

    pzi_lengths = fasta_lengths(args.pzi_reference)
    psp_lengths = fasta_lengths(args.psp_chr32)
    if set(pzi_lengths) != set(TARGET_ORDER):
        raise SystemExit("ERROR: PZI reference does not contain the current 32 phase-qualified targets")
    if set(psp_lengths) != set(assignment_by_query):
        missing = sorted(set(assignment_by_query) - set(psp_lengths))
        raise SystemExit(f"ERROR: PSP chr32 FASTA and assignment disagree; missing={missing}")

    hits = [hit for hit in parse_paf(args.paf, args.min_mapq, args.min_aln_len)
            if hit["tname"] in pzi_lengths and hit["qname"] in psp_lengths]
    if not hits:
        raise SystemExit("ERROR: no PAF alignments remain after filtering")

    intervals: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    for hit in hits:
        intervals[(hit["qname"], hit["tname"])].append((hit["qstart"], hit["qend"]))

    query_order = sorted(
        assignment_by_query,
        key=lambda name: TARGET_ORDER.index(assignment_by_query[name]["new_seqid"]),
    )
    matrix = []
    matrix_rows = []
    for query in query_order:
        row = []
        for target in TARGET_ORDER:
            covered = union_length(intervals.get((query, target), []))
            row.append(covered)
            if covered:
                matrix_rows.append(
                    {
                        "psp_seqid": query,
                        "assigned_target": assignment_by_query[query]["new_seqid"],
                        "pzi_target": target,
                        "aligned_union_bp": covered,
                        "psp_length_bp": psp_lengths[query],
                        "pzi_length_bp": pzi_lengths[target],
                        "psp_coverage_fraction": f"{covered / psp_lengths[query]:.6f}",
                        "pzi_coverage_fraction": f"{covered / pzi_lengths[target]:.6f}",
                        "phase": phase_of(target),
                    }
                )
        matrix.append(row)

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    fields = [
        "psp_seqid", "assigned_target", "pzi_target", "aligned_union_bp",
        "psp_length_bp", "pzi_length_bp", "psp_coverage_fraction",
        "pzi_coverage_fraction", "phase",
    ]
    with (out_dir / "psp_pzi_alignment_matrix.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(matrix_rows)

    audit_fields = [
        "psp_seqid", "assigned_target", "psp_length_bp", "assigned_aligned_union_bp",
        "assigned_coverage_fraction", "best_target", "best_aligned_union_bp",
        "best_coverage_fraction", "assignment_is_best",
    ]
    with (out_dir / "psp_pzi_alignment_audit.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=audit_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for query in query_order:
            values = {target: union_length(intervals.get((query, target), [])) for target in TARGET_ORDER}
            best = max(values, key=values.get)
            assigned = assignment_by_query[query]["new_seqid"]
            assigned_bp = values.get(assigned, 0)
            writer.writerow(
                {
                    "psp_seqid": query,
                    "assigned_target": assigned,
                    "psp_length_bp": psp_lengths[query],
                    "assigned_aligned_union_bp": assigned_bp,
                    "assigned_coverage_fraction": f"{assigned_bp / psp_lengths[query]:.6f}",
                    "best_target": best,
                    "best_aligned_union_bp": values[best],
                    "best_coverage_fraction": f"{values[best] / psp_lengths[query]:.6f}",
                    "assignment_is_best": "yes" if best == assigned else "no",
                }
            )

    configure_matplotlib()
    matrix_arr = [[math.log10(1 + value) for value in row] for row in matrix]
    vmax = max(1, max(max(row) for row in matrix_arr))
    fig, ax = plt.subplots(figsize=(13.5, 9.5))
    image = ax.imshow(
        matrix_arr, aspect="auto", interpolation="nearest", cmap="viridis",
        norm=Normalize(vmin=0, vmax=vmax),
    )
    ax.set_xticks(range(32), TARGET_ORDER, rotation=90, fontsize=5.5)
    ax.set_yticks(range(32), [assignment_by_query[q]["new_seqid"] for q in query_order], fontsize=5.5)
    ax.set_xlabel("PZI target chromosome")
    ax.set_ylabel("PSP chr32 chromosome (assigned label)")
    ax.set_title("PSP-to-PZI chromosome assignment from the current asm10 PAF")
    for boundary in (7.5, 15.5, 23.5):
        ax.axvline(boundary, color="white", lw=0.7)
        ax.axhline(boundary, color="white", lw=0.7)
    cbar = fig.colorbar(image, ax=ax, pad=0.02, fraction=0.035)
    cbar.set_label("log10(1 + union-aligned bp)")
    save_figure(fig, out_dir / "psp_vs_pzi_assignment_heatmap")

    t_offsets, t_bounds, t_total = cumulative_offsets(TARGET_ORDER, pzi_lengths)
    q_offsets, q_bounds, q_total = cumulative_offsets(query_order, psp_lengths)
    step = max(1, math.ceil(len(hits) / args.max_points))
    sampled = hits[::step]
    fig, ax = plt.subplots(figsize=(13.5, 8.0))
    for hit in sampled:
        x = t_offsets[hit["tname"]] + (hit["tstart"] + hit["tend"]) / 2
        y = q_offsets[hit["qname"]] + (hit["qstart"] + hit["qend"]) / 2
        ax.scatter(x, y, s=1.2, alpha=0.22, color=PHASE_COLORS.get(phase_of(hit["tname"]), "#777777"), linewidths=0)
    for position, _ in t_bounds[:-1]:
        ax.axvline(position, color="#bdbdbd", lw=0.35)
    for position, _ in q_bounds[:-1]:
        ax.axhline(position, color="#bdbdbd", lw=0.35)
    ax.set_xlim(0, t_total)
    ax.set_ylim(0, q_total)
    ax.invert_yaxis()
    ax.set_xlabel("PZI reference coordinate (phase-qualified chromosomes)")
    ax.set_ylabel("PSP chr32 coordinate (assignment order)")
    ax.set_title("Coordinate-level PSP-to-PZI alignments")
    ax.set_xticks([pos - pzi_lengths[name] / 2 for pos, name in t_bounds], [name for _, name in t_bounds], rotation=90, fontsize=5)
    ax.set_yticks([pos - psp_lengths[name] / 2 for pos, name in q_bounds], [assignment_by_query[name]["new_seqid"] for _, name in q_bounds], fontsize=5)
    save_figure(fig, out_dir / "psp_vs_pzi_coordinate_dotplot")

    fig, axes = plt.subplots(2, 1, figsize=(13.5, 6.5), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]})
    assigned_cov = []
    best_cov = []
    assignment_best = []
    for query in query_order:
        values = {target: union_length(intervals.get((query, target), [])) for target in TARGET_ORDER}
        assigned = assignment_by_query[query]["new_seqid"]
        assigned_cov.append(values[assigned] / psp_lengths[query])
        best_cov.append(max(values.values()) / psp_lengths[query])
        assignment_best.append(1 if assigned == max(values, key=values.get) else 0)
    x = range(32)
    axes[0].bar(x, best_cov, color="#777777", width=0.78, label="best target")
    axes[0].bar(x, assigned_cov, color="#0072B2", width=0.48, label="assigned target")
    axes[0].set_ylabel("PSP covered fraction")
    axes[0].set_ylim(0, 1.05)
    axes[0].legend(frameon=False, ncol=2, loc="upper right")
    axes[1].bar(x, assignment_best, color="#009E73", width=0.78)
    axes[1].set_yticks([0, 1], ["no", "yes"])
    axes[1].set_ylabel("Assigned = best")
    axes[1].set_xticks(x, [assignment_by_query[q]["new_seqid"] for q in query_order], rotation=90, fontsize=5.5)
    axes[1].set_xlabel("Current PSP assignment label")
    axes[0].set_title("Assignment coverage and best-target consistency")
    save_figure(fig, out_dir / "psp_vs_pzi_assignment_consistency")

    (out_dir / "README.txt").write_text(
        "These figures use only the current PSP asm10 PAF, current PSP assignment table, "
        "current repartitioned PZI reference, and current PSP chr32 FASTA. The heatmap "
        "uses query-interval unions per PSP/PZI pair to avoid double-counting secondary "
        "alignments. The dotplot retains all filtered PAF records for structural context.\n",
        encoding="utf-8",
    )
    print(f"[OK] Wrote PSP-PZI visualizations: {out_dir}")


if __name__ == "__main__":
    main()
