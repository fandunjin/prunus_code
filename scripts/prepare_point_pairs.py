#!/usr/bin/env python3
"""Prepare two-taxon POInT inputs from the formal six-track pillars."""

from __future__ import print_function

import argparse
import csv
from pathlib import Path


SPECIES = ("PZI", "PSP", "PPA")
PAIRS = (("PZI", "PSP"), ("PZI", "PPA"), ("PSP", "PPA"))


def read_rows(path):
    """Read the native POInT optimized file (header plus six copy columns)."""
    with Path(path).open("r", encoding="ascii") as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]
    if not lines or lines[0].split("\t") != list(SPECIES):
        raise ValueError("unexpected optimized POInT header")
    rows = []
    for line_number, line in enumerate(lines[1:], 2):
        fields = line.split("\t")
        if len(fields) != 2 * len(SPECIES):
            raise ValueError(
                "optimized pillar line {} has {} fields".format(
                    line_number, len(fields)
                )
            )
        rows.append(
            {
                species + "_A_gene": fields[index]
                for index, species in enumerate(SPECIES)
            }
        )
        rows[-1].update(
            {
                species + "_B_gene": fields[len(SPECIES) + index]
                for index, species in enumerate(SPECIES)
            }
        )
    return rows


def copy_geneorder(source, destination, name):
    lines = Path(source).read_text(encoding="ascii").splitlines()
    if not lines:
        raise ValueError("empty gene-order file: {}".format(source))
    lines[0] = name
    Path(destination).write_text("\n".join(lines) + "\n", encoding="ascii")


def write_tree(path, names):
    Path(path).write_text(
        "#NEXUS\n\nBegin trees;\nTranslate\n"
        "  1 {},\n  2 {}\n;\n"
        "tree fixed = [&R] (1:0.1,2:0.1):0.1;\nEnd;\n".format(*names),
        encoding="ascii",
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", required=True, help="shared POInT 01_inputs directory")
    ap.add_argument(
        "--optimized", required=True, help="six-track optimized pillar file"
    )
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit("ERROR: output directory is not empty: {}".format(outdir))
    outdir.mkdir(parents=True, exist_ok=True)
    rows = read_rows(args.optimized)
    if not rows:
        raise SystemExit("ERROR: optimized pillar file is empty")

    manifest = []
    for left, right in PAIRS:
        pair_name = left + "_" + right
        pdir = outdir / pair_name
        pdir.mkdir()
        names = (left, right)
        for species in names:
            copy_geneorder(
                Path(args.inputs) / (species + "_POInT_geneorder.txt"),
                pdir / (species + "_POInT_geneorder.txt"),
                species,
            )
        # POInT consumes copy 0 for all taxa followed by copy 1 for all taxa.
        with (pdir / "pillars.txt").open("w", encoding="ascii") as handle:
            handle.write("\t".join(names) + "\n")
            for row in rows:
                tokens = [
                    row[left + "_A_gene"] or "NONE",
                    row[right + "_A_gene"] or "NONE",
                    row[left + "_B_gene"] or "NONE",
                    row[right + "_B_gene"] or "NONE",
                ]
                handle.write("\t".join(tokens) + "\n")
        write_tree(pdir / "species_tree.nex", names)
        (pdir / "README.txt").write_text(
            "Two-taxon shared-WGD fit for {} and {}.\n"
            "The root Dupl state is shared by the pair; terminal branches are fitted separately.\n"
            "The resulting likelihood is combined with the independent singleton likelihood\n"
            "for the third species in the model comparison.\n".format(left, right),
            encoding="ascii",
        )
        manifest.append((pair_name, len(rows), left, right))
    with (outdir / "input_manifest.tsv").open(
        "w", encoding="ascii", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["pair", "pillars", "left", "right"])
        writer.writerows(manifest)
    (outdir / "PREPARED").write_text(
        "Prepared pairwise shared POInT inputs.\n", encoding="ascii"
    )
    print(
        "[OK] Prepared {} pairwise POInT inputs: {} pillars".format(
            len(PAIRS), len(rows)
        )
    )


if __name__ == "__main__":
    main()
