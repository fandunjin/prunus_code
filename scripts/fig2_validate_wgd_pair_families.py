#!/usr/bin/env python3
"""Validate pair-family CDS with the same Biopython rules used by WGD v2."""

from __future__ import print_function

import argparse
import csv
from pathlib import Path

from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError


def read_families(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if "family_id" not in fields:
            raise ValueError("Missing family_id column in {}".format(path))
        gene_fields = [field for field in fields if field != "family_id"]
        if len(gene_fields) != 2:
            raise ValueError("Expected exactly two gene columns in {}; found {}".format(path, gene_fields))
        rows = list(reader)
    if not rows:
        raise ValueError("No pair families in {}".format(path))
    return fields, gene_fields, rows


def write_tsv(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--families", type=Path, required=True)
    parser.add_argument("--cds", type=Path, action="append", required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    fields, gene_fields, families = read_families(args.families)
    required = set()
    for row in families:
        required.update(row[field] for field in gene_fields)

    gene_status = {}
    for cds_path in args.cds:
        if not cds_path.is_file():
            raise SystemExit("[ERROR] Missing CDS FASTA: {}".format(cds_path))
        for record in SeqIO.parse(str(cds_path), "fasta"):
            gene = record.id
            if gene not in required:
                continue
            if gene in gene_status:
                raise ValueError("Duplicate required CDS identifier across inputs: {}".format(gene))
            try:
                record.translate(to_stop=False, cds=True, stop_symbol="")
                gene_status[gene] = ("valid", "PASS")
            except TranslationError as error:
                gene_status[gene] = ("invalid_cds", str(error).replace("\t", " "))

    for gene in required - set(gene_status):
        gene_status[gene] = ("missing_from_cds", "Identifier was not found in supplied CDS FASTA")

    retained = []
    pair_audit = []
    for row in families:
        gene_1 = row[gene_fields[0]]
        gene_2 = row[gene_fields[1]]
        status_1, detail_1 = gene_status[gene_1]
        status_2, detail_2 = gene_status[gene_2]
        pair_status = "retained" if status_1 == status_2 == "valid" else "excluded_invalid_or_missing_cds"
        if pair_status == "retained":
            retained.append(row)
        pair_audit.append({
            "dataset": args.dataset,
            "family_id": row["family_id"],
            "gene_1": gene_1,
            "gene_1_status": status_1,
            "gene_1_detail": detail_1,
            "gene_2": gene_2,
            "gene_2_status": status_2,
            "gene_2_detail": detail_2,
            "pair_status": pair_status,
        })

    if not retained:
        raise SystemExit("[ERROR] No CDS-valid pair families retained for {}".format(args.dataset))

    gene_audit = []
    for gene in sorted(required):
        status, detail = gene_status[gene]
        gene_audit.append({"dataset": args.dataset, "gene": gene, "status": status, "detail": detail})

    write_tsv(args.outdir / (args.dataset + ".families.cds_valid.tsv"), fields, retained)
    write_tsv(
        args.outdir / (args.dataset + ".pair_cds_audit.tsv"),
        [
            "dataset", "family_id", "gene_1", "gene_1_status", "gene_1_detail",
            "gene_2", "gene_2_status", "gene_2_detail", "pair_status",
        ],
        pair_audit,
    )
    write_tsv(args.outdir / (args.dataset + ".gene_cds_audit.tsv"), ["dataset", "gene", "status", "detail"], gene_audit)

    valid_genes = sum(1 for status, _ in gene_status.values() if status == "valid")
    invalid_genes = sum(1 for status, _ in gene_status.values() if status == "invalid_cds")
    missing_genes = sum(1 for status, _ in gene_status.values() if status == "missing_from_cds")
    summary = [{
        "dataset": args.dataset,
        "families_input": len(families),
        "families_retained": len(retained),
        "families_excluded": len(families) - len(retained),
        "required_genes": len(required),
        "valid_genes": valid_genes,
        "invalid_cds_genes": invalid_genes,
        "missing_genes": missing_genes,
    }]
    write_tsv(
        args.outdir / (args.dataset + ".cds_validation_summary.tsv"),
        [
            "dataset", "families_input", "families_retained", "families_excluded",
            "required_genes", "valid_genes", "invalid_cds_genes", "missing_genes",
        ],
        summary,
    )
    (args.outdir / "COMPLETE").write_text("WGD CDS pair-family validation completed.\n", encoding="ascii")
    print("[OK] {}: retained {}/{} CDS-valid pair families; invalid genes={}, missing genes={}".format(
        args.dataset, len(retained), len(families), invalid_genes, missing_genes))


if __name__ == "__main__":
    main()
