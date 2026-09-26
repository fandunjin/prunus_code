#!/usr/bin/env python3
"""Compute non-redundant TE and TRF density in aligned 1-Mb windows.

TE coordinates must come from the same frozen pan-TE RepeatMasker run used
for the composition table.  RepeatMasker hits are merged before coverage is
assigned to windows, preventing nested hits from inflating the signal.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


TARGETS = [
    "pziA1", "pziA2", "pziB1", "pziB2",
    "pspA1", "pspA2", "pspB1", "pspB2",
    "ppaA", "ppaB", "pcam", "ppe",
]
GROUP = {t: ("PZI" if t.startswith("pzi") else "PSP" if t.startswith("psp")
             else "PPA" if t.startswith("ppa") else t.upper()) for t in TARGETS}


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    fields = fields or (list(rows[0]) if rows else [])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def target_from_path(path):
    text = str(path).lower()
    for target in TARGETS:
        if target.lower() in text:
            return target
    aliases = {"pzi_hapa1": "pziA1", "pzi_hapa2": "pziA2", "pzi_hapb1": "pziB1", "pzi_hapb2": "pziB2",
               "psp_hapa1": "pspA1", "psp_hapa2": "pspA2", "psp_hapb1": "pspB1", "psp_hapb2": "pspB2"}
    for alias, target in aliases.items():
        if alias in text:
            return target
    return ""


def parse_out(path):
    rows = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 11 or not fields[0].isdigit():
                continue
            try:
                seqid = fields[4]
                start = int(fields[5]) - 1
                end = int(fields[6])
            except (ValueError, IndexError):
                continue
            if end > start:
                rows.append((seqid, start, end, fields[9], fields[10]))
    return rows


def merge_intervals(intervals):
    intervals = sorted(intervals)
    merged = []
    for start, end in intervals:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        elif end > merged[-1][1]:
            merged[-1][1] = end
    return merged


def coverage_by_window(intervals, window_bp):
    merged = merge_intervals(intervals)
    coverage = defaultdict(int)
    for start, end in merged:
        while start < end:
            window = (start // window_bp) * window_bp
            stop = min(end, window + window_bp)
            coverage[window] += stop - start
            start = stop
    return coverage, merged


def rank(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        midpoint = (i + j - 1) / 2.0 + 1.0
        for k in order[i:j]:
            ranks[k] = midpoint
        i = j
    return ranks


def spearman(a, b):
    if len(a) < 3:
        return ""
    ra, rb = rank(a), rank(b)
    ma, mb = sum(ra) / len(ra), sum(rb) / len(rb)
    num = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    den_a = sum((x - ma) ** 2 for x in ra) ** 0.5
    den_b = sum((y - mb) ** 2 for y in rb) ** 0.5
    return f"{num / (den_a * den_b):.8f}" if den_a and den_b else ""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeatmasker-root", type=Path, required=True)
    parser.add_argument("--trf-density", type=Path, required=True)
    parser.add_argument("--te-composition", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--window-bp", type=int, default=1_000_000)
    args = parser.parse_args()

    trf_rows = read_tsv(args.trf_density)
    genome_rows = {row["target"]: int(float(row["genome_bp"])) for row in read_tsv(args.te_composition)}
    trf_by_target = defaultdict(dict)
    for row in trf_rows:
        target = row["target"]
        key = (row["seqid"], int(row["start"]), int(row["end"]))
        trf_by_target[target][key] = int(float(row["tandem_union_bp"]))

    files_by_target = defaultdict(list)
    for path in sorted(args.repeatmasker_root.rglob("*.out")):
        target = target_from_path(path)
        if target:
            files_by_target[target].append(path)
    missing = [target for target in TARGETS if not files_by_target[target]]
    if missing:
        raise SystemExit("ERROR: missing target-resolvable RepeatMasker .out files: " + ", ".join(missing))

    output = []
    corr = []
    for target in TARGETS:
        intervals_by_seq = defaultdict(set)
        for path in files_by_target[target]:
            for seqid, start, end, name, cls in parse_out(path):
                intervals_by_seq[seqid].add((start, end))
        te_window_by_seq = {}
        for seqid, intervals in intervals_by_seq.items():
            te_window_by_seq[seqid], _ = coverage_by_window(intervals, args.window_bp)

        target_trf = trf_by_target[target]
        te_values, trf_values = [], []
        for (seqid, start, end), trf_bp in sorted(target_trf.items()):
            window_bp = end - start
            te_bp = te_window_by_seq.get(seqid, {}).get(start, 0)
            output.append({"target": target, "group": GROUP[target], "seqid": seqid,
                           "start": start, "end": end, "window_bp": window_bp,
                           "te_union_bp": te_bp, "te_fraction": f"{te_bp / window_bp:.8f}",
                           "trf_union_bp": trf_bp, "trf_fraction": f"{trf_bp / window_bp:.8f}",
                           "te_trf_overlap_not_estimated": "1"})
            te_values.append(te_bp / window_bp)
            trf_values.append(trf_bp / window_bp)
        corr.append({"target": target, "group": GROUP[target], "window_count": len(te_values),
                     "spearman_te_trf_density": spearman(te_values, trf_values),
                     "mean_te_fraction": f"{sum(te_values) / len(te_values):.8f}" if te_values else "",
                     "mean_trf_fraction": f"{sum(trf_values) / len(trf_values):.8f}" if trf_values else ""})

    out = args.output_dir
    write_tsv(out / "te_trf_1mb_density.tsv", output)
    write_tsv(out / "te_trf_1mb_correlation.tsv", corr)
    (out / "te_trf_density_manifest.json").write_text(
        '{\n  "window_bp": %d,\n  "te_coverage": "merged RepeatMasker union per window",\n'
        '  "trf_coverage": "existing non-redundant TRF union per window",\n'
        '  "overlap": "not computed in this table; do not add TE and TRF fractions"\n}\n' % args.window_bp,
        encoding="utf-8",
    )
    print(f"Wrote {out / 'te_trf_1mb_density.tsv'}")


if __name__ == "__main__":
    main()
