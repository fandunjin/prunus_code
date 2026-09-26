#!/usr/bin/env python3
"""Prepare replicated-terminal POInT inputs for independent WGD fits.

POInT 1.63 expects a multi-terminal rooted tree and is unstable for a single
terminal. For each species we therefore make three identical terminal copies.
The internal and root branches are fixed to zero, so the three terminal
likelihoods are independent conditional on the fixed Dupl root and the total
lnL is exactly three times the one-species lnL (up to floating-point error).
Only one terminal branch and one set of model parameters are counted for the
effective one-species fit.
"""

from __future__ import print_function

import argparse
import csv
import re
import shutil
from pathlib import Path

SPECIES = ("PZI", "PSP", "PPA")


def read_pillars(path):
    with Path(path).open("r", encoding="ascii") as handle:
        rows = [line.rstrip("\n").split("\t") for line in handle if line.strip()]
    if not rows or rows[0] != list(SPECIES):
        raise ValueError("unexpected pillar header: {}".format(rows[0] if rows else []))
    expected = 2 * len(SPECIES)
    if any(len(row) != expected for row in rows[1:]):
        raise ValueError("pillar rows must contain {} fields".format(expected))
    return rows[1:]


def copy_geneorder(source, destination, name):
    lines = Path(source).read_text(encoding="ascii").splitlines()
    if not lines:
        raise ValueError("empty gene order: {}".format(source))
    lines[0] = name
    Path(destination).write_text("\n".join(lines) + "\n", encoding="ascii")


def write_tree(path, names):
    Path(path).write_text(
        "#NEXUS\n\nBegin trees;\nTranslate\n"
        "  1 {},\n  2 {},\n  3 {}\n;\n"
        "tree fixed = [&R] (1:0.1,(2:0.1,3:0.1):0.0):0.0;\nEnd;\n".format(*names),
        encoding="ascii",
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", required=True, help="shared POInT 01_inputs directory")
    ap.add_argument("--optimized", required=True, help="six-track optimized pillar file")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit("ERROR: output directory is not empty: {}".format(outdir))
    outdir.mkdir(parents=True, exist_ok=True)
    rows = read_pillars(args.optimized)
    manifest = []
    for species_index, species in enumerate(SPECIES):
        sdir = outdir / species
        sdir.mkdir()
        names = [species + "_R{}".format(i) for i in (1, 2, 3)]
        source_geneorder = Path(args.inputs) / (species + "_POInT_geneorder.txt")
        for name in names:
            copy_geneorder(source_geneorder, sdir / (name + "_geneorder.txt"), name)
        # POInT consumes all first-copy genes for all taxa, then all second-copy
        # genes. Repeat this species' A/B pair across the three terminals.
        a_col = species_index
        b_col = len(SPECIES) + species_index
        with (sdir / "pillars.txt").open("w", encoding="ascii") as handle:
            handle.write("\t".join(names) + "\n")
            for row in rows:
                a, b = row[a_col], row[b_col]
                handle.write("\t".join([a, a, a, b, b, b]) + "\n")
        write_tree(sdir / "species_tree.nex", names)
        (sdir / "README.txt").write_text(
            "Three identical terminals represent one independent {} WGD fit.\n"
            "Root and internal branches are fixed at zero; divide total lnL by 3.\n"
            "Effective parameter count is model parameters plus one unique terminal branch.\n".format(species),
            encoding="ascii",
        )
        manifest.append((species, len(rows), ",".join(names)))
    with (outdir / "input_manifest.tsv").open("w", encoding="ascii", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["species", "pillars", "replicated_terminals"])
        writer.writerows(manifest)
    (outdir / "PREPARED").write_text("Prepared replicated-terminal independent POInT inputs.\n", encoding="ascii")
    print("[OK] Prepared replicated independent inputs: {} pillars x 3 species".format(len(rows)))


if __name__ == "__main__":
    main()
