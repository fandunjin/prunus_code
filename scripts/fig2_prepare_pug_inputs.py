#!/usr/bin/env python3
"""Prepare traceable PUG inputs from the strict Fig. 2b family set.

PUG expects exact gene identifiers in its paralog file and species tokens that
are substrings of the gene identifiers.  This script copies only the formal
family trees into a clean directory, checks their taxon content, and writes
pooled and taxon-specific A-B paralog tables.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
from pathlib import Path


DEFAULT_EXPECTED_LABELS = {"pzi_A", "pzi_B", "psp_A", "psp_B", "ppa_A", "ppa_B", "pcam", "gtr", "rch", "ppe"}
PAIR_COLUMNS = {
    "pzi_ab": "pzi_A__pzi_B",
    "psp_ab": "psp_A__psp_B",
    "ppa_ab": "ppa_A__ppa_B",
}


def read_rows(path: Path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def split_pair(value: str, family: str, column: str):
    fields = value.split("|")
    if len(fields) != 2 or not all(fields):
        raise ValueError(f"{family}: malformed pair in {column}: {value}")
    return tuple(f"{part.split('__', 1)[0]}|{part}" for part in fields)


def leaf_tokens(newick: str):
    text = re.sub(r"\[.*?\]", "", newick)
    tokens = re.split(r"[(),;:]", text)
    leaves = []
    for token in tokens:
        token = token.strip()
        if not token or re.fullmatch(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", token):
            continue
        leaves.append(token)
    return leaves


def write_pairs(path: Path, rows, groups):
    with path.open("w", encoding="utf-8", newline="") as handle:
        for row in rows:
            for group in groups:
                left, right = split_pair(row[PAIR_COLUMNS[group]], row["orthogroup"], PAIR_COLUMNS[group])
                handle.write(f"{left}\t{right}\tshared_candidate\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--families", required=True, type=Path)
    parser.add_argument("--trees", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--species-tree", default="(rch,(gtr,((ppa,(psp,pzi)),(pcam,ppe))));")
    parser.add_argument(
        "--expected-label", action="append", default=None,
        help="Expected family-tree taxon label; repeat for the complete taxon set.",
    )
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    expected_labels = set(args.expected_label or DEFAULT_EXPECTED_LABELS)

    rows = read_rows(args.families)
    if not rows:
        raise SystemExit("No family rows found")
    required_columns = {"orthogroup", *PAIR_COLUMNS.values()}
    missing = required_columns - set(rows[0])
    if missing:
        raise SystemExit(f"Missing columns in family table: {','.join(sorted(missing))}")

    if args.outdir.exists() and any(args.outdir.iterdir()):
        if not args.force:
            raise SystemExit(f"Output directory is not empty: {args.outdir}; use --force to rebuild")
        shutil.rmtree(args.outdir)
    tree_out = args.outdir / "trees"
    tree_out.mkdir(parents=True, exist_ok=True)

    manifest = []
    for row in rows:
        og = row["orthogroup"]
        source = args.trees / f"{og}.tre"
        if not source.is_file() or source.stat().st_size == 0:
            raise SystemExit(f"Missing formal tree for {og}: {source}")
        text = source.read_text(encoding="utf-8")
        leaves = leaf_tokens(text)
        labels = {leaf.split("|", 1)[0] for leaf in leaves}
        if labels != expected_labels:
            raise SystemExit(
                f"{og}: expected labels {sorted(expected_labels)}, observed {sorted(labels)}"
            )
        if len(leaves) != len(set(leaves)):
            raise SystemExit(f"{og}: duplicate leaf identifier in {source}")
        target = tree_out / source.name
        shutil.copy2(source, target)
        manifest.append({"orthogroup": og, "source_tree": str(source), "prepared_tree": str(target), "tip_count": len(leaves)})

    with (args.outdir / "tree_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest[0], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest)

    (args.outdir / "species_tree_for_pug.tre").write_text(args.species_tree + "\n", encoding="ascii")
    for group, column in PAIR_COLUMNS.items():
        path = args.outdir / f"paralogs_{group}.tsv"
        with path.open("w", encoding="utf-8", newline="") as handle:
            for row in rows:
                left, right = split_pair(row[column], row["orthogroup"], column)
                handle.write(f"{left}\t{right}\t{group}\n")
    write_pairs(args.outdir / "paralogs_shared_candidate.tsv", rows, list(PAIR_COLUMNS))

    (args.outdir / "README.txt").write_text(
        "PUG input prepared from the strict Fig. 2b family selection.\n"
        f"families={len(rows)}\n"
        "tree_inference=FastTree_default; support values are local support values in 0-1.\n"
        "PUG tree_type=Bayesian is used only to multiply 0-1 support values by 100; values are not Bayesian posterior probabilities.\n"
        "PUG outgroup=rch; species tree labels are lowercase substrings of gene identifiers.\n"
        f"expected_family_tree_labels={','.join(sorted(expected_labels))}\n"
        "The pooled table tests a shared candidate event; taxon-specific tables are diagnostic and are not independent biological replicates.\n",
        encoding="ascii",
    )
    print(f"[OK] Prepared {len(rows)} PUG families in {args.outdir}")


if __name__ == "__main__":
    main()
