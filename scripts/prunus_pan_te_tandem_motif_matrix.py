#!/usr/bin/env python3
"""Build a complete cross-genome matrix for short tandem-repeat motifs.

The earlier visualisation treated a motif that was outside a target's local
top-N list as zero.  This program resolves every TRF-overlap base once, then
constructs one shared motif catalogue before calculating species means.
Only units <= ``--max-period-bp`` are considered sequence motifs here.  Longer
TRF periods are retained in the period-class and array analyses because their
raw consensus strings are not comparable satellite-family identifiers.
"""

import argparse
import collections
import csv
import json
import statistics
from pathlib import Path


TARGETS = ["pziA1", "pziA2", "pziB1", "pziB2", "pspA1", "pspA2", "pspB1", "pspB2", "ppaA", "ppaB", "pcam", "ppe"]
GROUPS = {**{x: "PZI" for x in TARGETS[:4]}, **{x: "PSP" for x in TARGETS[4:8]}, "ppaA": "PPA", "ppaB": "PPA", "pcam": "PCAM", "ppe": "PPE"}
ORDER = ["PZI", "PSP", "PPA", "PCAM", "PPE"]


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def shortest_unit(sequence):
    for width in range(1, len(sequence) + 1):
        if len(sequence) % width == 0 and sequence == sequence[:width] * (len(sequence) // width):
            return sequence[:width]
    return sequence


def canonical(sequence):
    sequence = "".join(base for base in sequence.upper() if base in "ACGT")
    if not sequence:
        return "unknown"
    sequence = shortest_unit(sequence)
    complement = str.maketrans("ACGT", "TGCA")
    reverse = sequence.translate(complement)[::-1]
    candidates = [sequence[i:] + sequence[:i] for i in range(len(sequence))]
    candidates += [reverse[i:] + reverse[:i] for i in range(len(reverse))]
    return min(candidates)


def read_calls(path, min_score, max_period):
    output = collections.defaultdict(list)
    with open(path, encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            try:
                start, end = int(row["start"]), int(row["end"])
                score, period = int(row["score"]), int(row["period_bp"])
            except (KeyError, ValueError):
                continue
            if end <= start or score < min_score or period > max_period:
                continue
            output[row["seqid"]].append((start, end, score, period, canonical(row.get("consensus_motif", ""))))
    return output


def resolve(calls):
    events = collections.defaultdict(list)
    for index, row in enumerate(calls):
        events[row[0]].append((1, index)); events[row[1]].append((-1, index))
    active, total, last = set(), collections.Counter(), None
    for position in sorted(events):
        if last is not None and position > last and active:
            best = min(active, key=lambda i: (-calls[i][2], -(calls[i][1] - calls[i][0]), calls[i][4]))
            total[(calls[best][3], calls[best][4])] += position - last
        for sign, index in events[position]:
            if sign == 1: active.add(index)
            else: active.discard(index)
        last = position
    return total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tandem-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-score", type=int, default=50)
    parser.add_argument("--max-period-bp", type=int, default=100)
    parser.add_argument("--min-total-union-bp", type=int, default=100000)
    args = parser.parse_args()
    values = {}
    for target in TARGETS:
        path = args.tandem_root / target / f"{target}.trf.tsv"
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"ERROR: missing {path}")
        counts = collections.Counter()
        for calls in read_calls(path, args.min_score, args.max_period).values():
            counts.update(resolve(calls))
        values[target] = counts
    catalogue = sorted(set().union(*(set(v) for v in values.values())))
    catalogue = [key for key in catalogue if sum(values[target].get(key, 0) for target in TARGETS) >= args.min_total_union_bp]
    target_rows = []
    for period, motif in catalogue:
        for target in TARGETS:
            target_rows.append({"target": target, "group": GROUPS[target], "period_bp": period, "canonical_motif": motif, "exclusive_union_bp": values[target].get((period, motif), 0)})
    group_rows = []
    for period, motif in catalogue:
        out = {"period_bp": period, "canonical_motif": motif}
        for group in ORDER:
            targets = [target for target in TARGETS if GROUPS[target] == group]
            out[f"{group.lower()}_mean_bp"] = round(statistics.mean(values[target].get((period, motif), 0) for target in targets))
        out["pzi_minus_psp_mean_bp"] = out["pzi_mean_bp"] - out["psp_mean_bp"]
        out["total_union_bp"] = sum(values[target].get((period, motif), 0) for target in TARGETS)
        group_rows.append(out)
    group_rows.sort(key=lambda row: row["pzi_minus_psp_mean_bp"], reverse=True)
    write_tsv(args.output_dir / "tandem_short_motif_target_matrix.tsv", ["target", "group", "period_bp", "canonical_motif", "exclusive_union_bp"], target_rows)
    write_tsv(args.output_dir / "tandem_short_motif_group_matrix.tsv", ["period_bp", "canonical_motif", "pzi_mean_bp", "psp_mean_bp", "ppa_mean_bp", "pcam_mean_bp", "ppe_mean_bp", "pzi_minus_psp_mean_bp", "total_union_bp"], group_rows)
    with open(args.output_dir / "tandem_short_motif_manifest.json", "w", encoding="utf-8") as handle:
        json.dump({"trf_min_score": args.min_score, "max_period_bp": args.max_period_bp, "min_total_union_bp": args.min_total_union_bp, "definition": "complete cross-genome motif catalogue; exact canonical rotational/reverse-complement motif matching for short TRF units only"}, handle, indent=2)
    print(json.dumps({"motifs": len(group_rows), "output_dir": str(args.output_dir)}))


if __name__ == "__main__":
    main()
