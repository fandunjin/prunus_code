#!/usr/bin/env python3
"""Extract focal A/B anchor pairs from i-ADHoRe output."""

from __future__ import print_function

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--iadhore-output", required=True, type=Path)
    parser.add_argument("--stable-map", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    multi = args.iadhore_output / "multiplicons.txt"
    anchors = args.iadhore_output / "anchorpoints.txt"
    if not multi.is_file() or not anchors.is_file():
        raise SystemExit("ERROR: i-ADHoRe output lacks multiplicons.txt or anchorpoints.txt")
    stable_to_original = {}
    for row in read_tsv(args.stable_map):
        stable_to_original[row["stable_id"]] = row["original_id"]

    cross = {}
    for row in read_tsv(multi):
        lx, ly = row["list_x"], row["list_y"]
        if row["genome_x"] != row["genome_y"]:
            continue
        if lx.startswith("A_") and ly.startswith("B_"):
            cross[row["id"]] = (row["genome_x"], False)
        elif lx.startswith("B_") and ly.startswith("A_"):
            cross[row["id"]] = (row["genome_x"], True)
    if not cross:
        raise SystemExit("ERROR: no within-species A/B multiplicons found")

    pairs = defaultdict(set)
    audit = []
    with anchors.open("r", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"multiplicon", "gene_x", "gene_y", "is_real_anchorpoint"}
        if not required.issubset(reader.fieldnames or []):
            raise SystemExit("ERROR: unexpected anchorpoints columns: {}".format(reader.fieldnames))
        for row in reader:
            info = cross.get(row["multiplicon"])
            if info is None or row["is_real_anchorpoint"] in ("0", "false", "False"):
                continue
            species, reverse = info
            left, right = row["gene_x"], row["gene_y"]
            if reverse:
                left, right = right, left
            if left not in stable_to_original or right not in stable_to_original:
                raise SystemExit("ERROR: anchor gene absent from stable map: {} {}".format(left, right))
            pairs[species].add((left, right))
            audit.append({"species": species, "multiplicon": row["multiplicon"],
                          "gene_a_stable": left, "gene_b_stable": right,
                          "gene_a_original": stable_to_original[left],
                          "gene_b_original": stable_to_original[right]})

    args.outdir.mkdir(parents=True, exist_ok=True)
    name_map = {"PZI": "pzi_hapA1_vs_pzi_hapB1", "PSP": "psp_hapA1_vs_psp_hapB1",
                "PPA": "ppa_hapA_vs_ppa_hapB"}
    for species, rows in pairs.items():
        if species not in name_map or not rows:
            continue
        folder = args.outdir / name_map[species]
        folder.mkdir(parents=True, exist_ok=True)
        stem = name_map[species].replace("_vs_", ".")
        with (folder / (stem + ".anchors.new")).open("w", encoding="utf-8", newline="\n") as handle:
            for left, right in sorted(rows):
                handle.write("{}\t{}\n".format(stable_to_original[left], stable_to_original[right]))
    with (args.outdir / "iadhore_ab_anchor_audit.tsv").open("w", encoding="utf-8", newline="") as handle:
        fields = ["species", "multiplicon", "gene_a_stable", "gene_b_stable", "gene_a_original", "gene_b_original"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(audit)
    with (args.outdir / "anchor_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        handle.write("species\tunique_ab_pairs\n")
        for species in ("PZI", "PSP", "PPA"):
            handle.write("{}\t{}\n".format(species, len(pairs.get(species, set()))))
    (args.outdir / "COMPLETE").write_text(
        "A/B pairs extracted from i-ADHoRe anchorpoints; these are not JCVI anchors.\n", encoding="ascii")
    print("[OK] Extracted i-ADHoRe A/B pairs: {}".format(", ".join(
        "{}={}".format(k, len(v)) for k, v in sorted(pairs.items()))))


if __name__ == "__main__":
    main()
