#!/usr/bin/env python3
"""Quantify non-redundant overlap between TRF arrays and RepeatMasker hits.

This is a classification-diagnostic analysis.  A high overlap says that a
TRF-defined tandem array is covered by the current pan-TE RepeatMasker library;
it does not independently establish that the array is a transposable element.
"""

import argparse
import collections
import csv
from pathlib import Path


TARGETS = ["pziA1", "pziA2", "pziB1", "pziB2", "pspA1", "pspA2", "pspB1", "pspB2", "ppaA", "ppaB", "pcam", "ppe"]
GROUP = {**{x: "PZI" for x in TARGETS[:4]}, **{x: "PSP" for x in TARGETS[4:8]}, "ppaA": "PPA", "ppaB": "PPA", "pcam": "PCAM", "ppe": "PPE"}


def merge(values):
    if not values:
        return []
    result = []
    for start, end in sorted(values):
        if not result or start > result[-1][1]:
            result.append([start, end])
        else:
            result[-1][1] = max(result[-1][1], end)
    return [(a, b) for a, b in result]


def overlap(a, b):
    i = j = total = 0
    while i < len(a) and j < len(b):
        total += max(0, min(a[i][1], b[j][1]) - max(a[i][0], b[j][0]))
        if a[i][1] <= b[j][1]:
            i += 1
        else:
            j += 1
    return total


def trf_intervals(path):
    result = collections.defaultdict(list)
    with open(path, encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            try:
                start, end = int(row["start"]), int(row["end"])
            except (ValueError, KeyError):
                continue
            if end > start:
                result[row["seqid"]].append((start, end))
    return {seq: merge(values) for seq, values in result.items()}


def repeatmasker_intervals(path):
    result = collections.defaultdict(list)
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 11 or not fields[0].isdigit():
                continue
            try:
                start, end = int(fields[5]) - 1, int(fields[6])
            except ValueError:
                continue
            if end > start:
                result[fields[4]].append((start, end))
    return {seq: merge(values) for seq, values in result.items()}


def write(path, rows):
    fields = ["target", "group", "tandem_union_bp", "repeatmasker_overlap_bp", "repeatmasker_overlap_fraction", "tandem_only_bp", "tandem_only_fraction"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); w.writeheader(); w.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tandem-root", type=Path, required=True)
    p.add_argument("--mask-root", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    rows = []
    for target in TARGETS:
        trf = args.tandem_root / target / f"{target}.trf.tsv"
        rmout = args.mask_root / target / f"{target}.fa.out"
        if not trf.is_file() or not rmout.is_file():
            raise SystemExit(f"ERROR: missing TRF or RepeatMasker output for {target}")
        tandem, repeats = trf_intervals(trf), repeatmasker_intervals(rmout)
        tandem_bp = sum(sum(b - a for a, b in values) for values in tandem.values())
        shared_bp = sum(overlap(values, repeats.get(seq, [])) for seq, values in tandem.items())
        rows.append({"target": target, "group": GROUP[target], "tandem_union_bp": tandem_bp, "repeatmasker_overlap_bp": shared_bp, "repeatmasker_overlap_fraction": f"{shared_bp/tandem_bp:.8f}" if tandem_bp else "0", "tandem_only_bp": tandem_bp - shared_bp, "tandem_only_fraction": f"{(tandem_bp-shared_bp)/tandem_bp:.8f}" if tandem_bp else "0"})
    write(args.output, rows)
    print(f"Wrote {len(rows)} target rows to {args.output}")


if __name__ == "__main__":
    main()
