#!/usr/bin/env python3
"""Build a one-genome i-ADHoRe 3.0 self-synteny input from MCL families."""

from __future__ import annotations

import argparse
import csv
from collections import OrderedDict
from pathlib import Path


def read_bed(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise ValueError("Malformed BED at {}:{}".format(path, line_number))
            rows.append({"seqid": fields[0], "gene": fields[3], "strand": fields[5]})
    if not rows:
        raise ValueError("No BED rows in {}".format(path))
    return rows


def read_mcl(path: Path) -> tuple[dict[str, str], dict[str, int]]:
    gene_to_family: dict[str, str] = {}
    family_sizes: dict[str, int] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for index, line in enumerate(handle, 1):
            genes = [gene for gene in line.rstrip("\n").split("\t") if gene]
            if not genes:
                continue
            family = "MCL_{:08d}".format(index)
            family_sizes[family] = len(genes)
            for gene in genes:
                if gene in gene_to_family:
                    raise ValueError("Gene occurs in more than one MCL family: {}".format(gene))
                gene_to_family[gene] = family
    if not family_sizes:
        raise ValueError("No MCL families in {}".format(path))
    return gene_to_family, family_sizes


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label", required=True)
    parser.add_argument("--bed", type=Path, required=True)
    parser.add_argument("--mcl", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=12)
    args = parser.parse_args()
    if args.outdir.exists() and any(args.outdir.iterdir()):
        raise SystemExit("[ERROR] Output directory is not empty: {}".format(args.outdir))
    if args.threads < 1:
        raise SystemExit("[ERROR] --threads must be positive")
    for path in (args.bed, args.mcl):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit("[ERROR] Missing input: {}".format(path))
    args.outdir.mkdir(parents=True, exist_ok=True)
    lists_dir = args.outdir / "lists"
    lists_dir.mkdir()

    bed_rows = read_bed(args.bed)
    gene_to_family, family_sizes = read_mcl(args.mcl)
    duplicate_families = {family for family, size in family_sizes.items() if size >= 2}
    by_seqid: OrderedDict[str, list[tuple[str, str]]] = OrderedDict()
    retained = []
    for row in bed_rows:
        family = gene_to_family.get(row["gene"], "")
        if family not in duplicate_families:
            continue
        by_seqid.setdefault(row["seqid"], []).append((row["gene"], row["strand"]))
        retained.append((row["gene"], family))
    if len(by_seqid) < 2 or len(retained) < 100:
        raise SystemExit("[ERROR] Insufficient cross-chromosome self-synteny input for {}".format(args.label))

    list_rows = []
    for index, (seqid, genes) in enumerate(by_seqid.items(), 1):
        list_name = "L{:03d}".format(index)
        path = lists_dir / (list_name + ".lst")
        with path.open("w", encoding="utf-8", newline="\n") as handle:
            for gene, strand in genes:
                handle.write("{}{}\n".format(gene, strand if strand in {"+", "-"} else "+"))
        list_rows.append({"label": args.label, "list_name": list_name, "seqid": seqid, "genes": len(genes), "path": str(path.resolve())})
    with (args.outdir / "families.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        for gene, family in sorted(retained):
            handle.write("{}\t{}\n".format(gene, family))

    output_dir = args.outdir / "iadhore_output"
    with (args.outdir / "iadhore.ini").open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("genome= {}\n".format(args.label))
        for row in list_rows:
            handle.write("{} {}\n".format(row["list_name"], row["path"]))
        handle.write("\nblast_table= {}\n".format((args.outdir / "families.tsv").resolve()))
        handle.write("table_type= family\n")
        handle.write("output_path= {}\n".format(output_dir.resolve()))
        handle.write("gap_size= 30\ncluster_gap= 35\nq_value= 0.75\nprob_cutoff= 0.01\n")
        handle.write("anchor_points= 3\nalignment_method= gg2\nlevel_2_only= false\n")
        handle.write("multiple_hypothesis_correction= FDR\nvisualizeGHM= false\nnumber_of_threads= {}\n".format(args.threads))
    write_tsv(args.outdir / "list_manifest.tsv", ["label", "list_name", "seqid", "genes", "path"], list_rows)
    (args.outdir / "method_contract.txt").write_text(
        "i-ADHoRe 3.0 self-synteny input using BLASTP/MCL family labels.\n"
        "Parameters: gap_size=30; cluster_gap=35; q_value=0.75; prob_cutoff=0.01; anchor_points=3; gg2; FDR.\n"
        "Only MCL families containing at least two genes enter i-ADHoRe.\n", encoding="ascii")
    (args.outdir / "COMPLETE").write_text("i-ADHoRe self-synteny input prepared.\n", encoding="ascii")
    print("[OK] {}: {} cross-chromosome lists; {} genes in {} duplicate MCL families".format(
        args.label, len(list_rows), len(retained), len(duplicate_families)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
