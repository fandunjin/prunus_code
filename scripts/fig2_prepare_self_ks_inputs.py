#!/usr/bin/env python3
"""Prepare one genome for conservative self-syntenic Ks analysis.

The script preserves only BED-ordered gene models with a matched, non-empty
protein and an in-frame CDS.  Strict biological CDS validation is performed
again by ``fig2_validate_wgd_pair_families.py`` after i-ADHoRe pairing.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import OrderedDict
from pathlib import Path


def read_fasta(path: Path) -> OrderedDict[str, str]:
    records: OrderedDict[str, str] = OrderedDict()
    identifier = None
    sequence: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    records[identifier] = "".join(sequence).upper()
                identifier = line[1:].split(None, 1)[0]
                if not identifier:
                    raise ValueError("Empty FASTA identifier at {}:{}".format(path, line_number))
                if identifier in records:
                    raise ValueError("Duplicate FASTA identifier {} in {}".format(identifier, path))
                sequence = []
            else:
                if identifier is None:
                    raise ValueError("FASTA sequence before header in {}".format(path))
                sequence.append(re.sub(r"\s+", "", line))
    if identifier is not None:
        records[identifier] = "".join(sequence).upper()
    if not records:
        raise ValueError("No FASTA records in {}".format(path))
    return records


def read_bed(path: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    seen: set[str] = set()
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise ValueError("Malformed BED at {}:{}".format(path, line_number))
            gene = fields[3]
            if gene in seen:
                raise ValueError("Duplicate BED gene {} in {}".format(gene, path))
            seen.add(gene)
            rows.append(fields[:6])
    if not rows:
        raise ValueError("No gene records in {}".format(path))
    return rows


def write_fasta(path: Path, records: OrderedDict[str, str]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for identifier, sequence in records.items():
            handle.write(">{}\n".format(identifier))
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label", required=True)
    parser.add_argument("--bed", type=Path, action="append", required=True,
                        help="One or more component BED files in genome order.")
    parser.add_argument("--pep", type=Path, action="append", required=True,
                        help="Protein FASTA files corresponding one-to-one with --bed.")
    parser.add_argument("--cds", type=Path, action="append", required=True,
                        help="CDS FASTA files corresponding one-to-one with --bed.")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    if args.outdir.exists() and any(args.outdir.iterdir()):
        raise SystemExit("[ERROR] Output directory is not empty: {}".format(args.outdir))
    if not (len(args.bed) == len(args.pep) == len(args.cds)):
        raise SystemExit("[ERROR] --bed, --pep and --cds must have equal numbers of components")
    for path in args.bed + args.pep + args.cds:
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit("[ERROR] Missing input: {}".format(path))
    args.outdir.mkdir(parents=True, exist_ok=True)

    bed_rows: list[list[str]] = []
    peptides: OrderedDict[str, str] = OrderedDict()
    cds: OrderedDict[str, str] = OrderedDict()
    component_rows: list[dict[str, object]] = []
    for index, (bed_path, pep_path, cds_path) in enumerate(zip(args.bed, args.pep, args.cds), 1):
        component_bed = read_bed(bed_path)
        component_pep = read_fasta(pep_path)
        component_cds = read_fasta(cds_path)
        previous_genes = {row[3] for row in bed_rows}
        duplicate_genes = previous_genes & {row[3] for row in component_bed}
        if duplicate_genes:
            raise ValueError("Duplicate gene IDs across input components: {}".format(
                ",".join(sorted(duplicate_genes)[:5])))
        duplicate_pep = set(peptides) & set(component_pep)
        duplicate_cds = set(cds) & set(component_cds)
        if duplicate_pep or duplicate_cds:
            raise ValueError("Duplicate FASTA IDs across input components")
        bed_rows.extend(component_bed)
        peptides.update(component_pep)
        cds.update(component_cds)
        component_rows.append({"component": index, "bed": str(bed_path), "pep": str(pep_path),
                               "cds": str(cds_path), "bed_models": len(component_bed),
                               "peptide_records": len(component_pep), "cds_records": len(component_cds)})

    retained_bed: list[list[str]] = []
    retained_pep: OrderedDict[str, str] = OrderedDict()
    retained_cds: OrderedDict[str, str] = OrderedDict()
    status_rows: list[dict[str, object]] = []
    counts: dict[str, int] = {}
    for fields in bed_rows:
        gene = fields[3]
        protein = peptides.get(gene, "")
        coding = cds.get(gene, "")
        if not protein:
            status = "excluded_missing_or_empty_peptide"
        elif not coding:
            status = "excluded_missing_or_empty_cds"
        elif len(coding) % 3 != 0:
            status = "excluded_nontriplet_cds"
        else:
            status = "retained_pre_iadhore"
            retained_bed.append(fields)
            retained_pep[gene] = protein
            retained_cds[gene] = coding
        counts[status] = counts.get(status, 0) + 1
        status_rows.append({"label": args.label, "gene": gene, "status": status})

    if len(retained_bed) < 100:
        raise SystemExit("[ERROR] Too few usable gene models retained for {}".format(args.label))

    with (args.outdir / "{}.bed".format(args.label)).open("w", encoding="utf-8", newline="\n") as handle:
        for fields in retained_bed:
            handle.write("\t".join(fields) + "\n")
    write_fasta(args.outdir / "{}.pep.fa".format(args.label), retained_pep)
    write_fasta(args.outdir / "{}.cds.fa".format(args.label), retained_cds)
    write_tsv(args.outdir / "gene_input_audit.tsv", ["label", "gene", "status"], status_rows)
    write_tsv(args.outdir / "component_manifest.tsv",
              ["component", "bed", "pep", "cds", "bed_models", "peptide_records", "cds_records"], component_rows)
    summary = [{"label": args.label, "bed_models_input": len(bed_rows), "peptide_records_input": len(peptides),
                "cds_records_input": len(cds), "models_retained_pre_iadhore": len(retained_bed),
                "excluded_missing_or_empty_peptide": counts.get("excluded_missing_or_empty_peptide", 0),
                "excluded_missing_or_empty_cds": counts.get("excluded_missing_or_empty_cds", 0),
                "excluded_nontriplet_cds": counts.get("excluded_nontriplet_cds", 0)}]
    write_tsv(args.outdir / "input_summary.tsv", list(summary[0]), summary)
    (args.outdir / "method_contract.txt").write_text(
        "Genome-native self-syntenic Ks prefilter. One or more explicitly listed genome components are concatenated without altering gene or chromosome identifiers.\n"
        "Only BED-ordered models with matching non-empty protein and in-frame CDS are retained.\n"
        "Strict CDS translation validation occurs after i-ADHoRe anchor extraction.\n", encoding="ascii")
    (args.outdir / "COMPLETE").write_text("Self-syntenic Ks inputs prepared.\n", encoding="ascii")
    print("[OK] {}: retained {}/{} BED-ordered genes".format(args.label, len(retained_bed), len(bed_rows)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
