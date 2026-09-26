#!/usr/bin/env python3
"""Map POInT_ances_order rows back to their original GTR-centered loci."""

from __future__ import print_function

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--initial-pillars", required=True)
    parser.add_argument("--initial-map", required=True)
    parser.add_argument("--optimized-pillars", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    initial_lines = Path(args.initial_pillars).read_text(encoding="ascii").splitlines()
    optimized_lines = Path(args.optimized_pillars).read_text(encoding="ascii").splitlines()
    if not initial_lines or not optimized_lines:
        raise SystemExit("ERROR: empty POInT pillar file")
    initial_map = read_tsv(args.initial_map)
    if len(initial_map) != len(initial_lines) - 1:
        raise SystemExit("ERROR: initial pillar map does not match initial pillar count")
    lookup = {}
    for index, line in enumerate(initial_lines[1:]):
        key = tuple(line.split("\t"))
        if key in lookup:
            raise SystemExit("ERROR: duplicate initial POInT pillar key at {}".format(index))
        lookup[key] = initial_map[index]
    rows = []
    for index, line in enumerate(optimized_lines[1:]):
        key = tuple(line.split("\t"))
        if key not in lookup:
            raise SystemExit("ERROR: optimized pillar {} cannot be mapped to initial GTR locus".format(index))
        row = dict(lookup[key])
        row["optimized_index"] = str(index)
        rows.append(row)
    if len(rows) != len(initial_map):
        raise SystemExit("ERROR: optimized pillar count differs from initial pillar count")
    write_tsv(args.out, ["optimized_index", "pillar_id", "gtr_global_index", "gtr_seqid", "gtr_gene"], rows)
    print("[OK] Mapped {} optimized POInT pillars to GTR loci".format(len(rows)))


if __name__ == "__main__":
    main()
