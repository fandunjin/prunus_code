#!/usr/bin/env python3
"""Extract high-occupancy, single-copy orthogroups from OrthoFinder output."""

import argparse
import csv
import math
from collections import Counter
from pathlib import Path


FOCAL = set(["pzi_A", "pzi_B", "psp_A", "psp_B", "ppa_A", "ppa_B"])


def read_fasta(path):
    records = {}
    name = None
    chunks = []
    with path.open("r", encoding="ascii") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        records[name] = "".join(chunks)
    return records


def split_genes(cell):
    return [item.strip() for item in cell.split(",") if item.strip()]


def read_orthogroups(path):
    """Read OrthoFinder 3.x TSV or legacy colon-separated Orthogroups.txt."""
    first_line = path.open("r", encoding="utf-8").readline()
    if first_line.startswith("Orthogroup\t"):
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = reader.fieldnames or []
            rows = []
            for row in reader:
                rows.append((row[fields[0]], dict((field, split_genes(row.get(field, ""))) for field in fields[1:])))
            return fields[1:], rows

    # OrthoFinder 3.1 writes Orthogroups.txt as `OGxxxx: gene gene ...`.
    rows = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if ":" not in line:
                raise SystemExit("ERROR: unsupported Orthogroups format: %s" % path)
            og, genes_text = line.split(":", 1)
            grouped = {}
            for gene in genes_text.split():
                if "__" not in gene:
                    raise SystemExit("ERROR: cannot infer taxon from gene ID %s" % gene)
                species = gene.split("__", 1)[0]
                grouped.setdefault(species, []).append(gene)
            rows.append((og.strip(), grouped))
    species = sorted({species for _, grouped in rows for species in grouped})
    return species, rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--orthogroups", required=True)
    parser.add_argument("--pep-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--min-occupancy", type=float, default=1.0)
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--required", default="pzi_A,pzi_B,psp_A,psp_B,ppa_A,ppa_B,gtr")
    args = parser.parse_args()

    table = Path(args.orthogroups)
    pep_dir = Path(args.pep_dir)
    out_dir = Path(args.out_dir)
    unaligned = out_dir / "unaligned"
    unaligned.mkdir(parents=True, exist_ok=False)

    pep_paths = sorted(pep_dir.glob("*.fa"))
    species = [path.stem for path in pep_paths]
    indexes = dict((sp, read_fasta(pep_dir / (sp + ".fa"))) for sp in species)
    required = set(item for item in args.required.split(",") if item)
    missing_required = required - set(species)
    if missing_required:
        raise SystemExit("ERROR: required taxa absent from protein directory: %s" % ",".join(sorted(missing_required)))
    minimum_taxa = int(math.ceil(args.min_occupancy * len(species)))

    selected = []
    rejected = Counter()
    table_species, rows = read_orthogroups(table)
    if set(species) != set(table_species):
        raise SystemExit("ERROR: OrthoFinder taxa and frozen protein files differ")
    for og, grouped in rows:
        genes = dict((sp, grouped.get(sp, [])) for sp in species)
        if any(len(values) > 1 for values in genes.values()):
            rejected["multicopy"] += 1
            continue
        present = set(sp for sp, values in genes.items() if len(values) == 1)
        if not required.issubset(present):
            rejected["missing_required_taxon"] += 1
            continue
        if len(present) < minimum_taxa:
            rejected["low_occupancy"] += 1
            continue

        output = unaligned / (og + ".fa")
        with output.open("w", encoding="ascii", newline="\n") as fasta:
            for sp in species:
                if sp not in present:
                    continue
                gene = genes[sp][0]
                if gene not in indexes[sp]:
                    raise SystemExit("ERROR: %s from %s is absent from frozen FASTA" % (gene, sp))
                fasta.write(">%s\n%s\n" % (sp, indexes[sp][gene]))
        selected.append((og, len(present), ",".join(sorted(set(species) - present))))

    with (out_dir / "selected_orthogroups.tsv").open("w", encoding="utf-8", newline="") as handle:
        handle.write("orthogroup\ttaxa_present\tmissing_taxa\n")
        for og, count, missing in selected:
            handle.write("%s\t%d\t%s\n" % (og, count, missing))
    with (out_dir / "selection_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        handle.write("metric\tvalue\n")
        handle.write("total_taxa\t%d\n" % len(species))
        handle.write("minimum_taxa\t%d\n" % minimum_taxa)
        handle.write("selected\t%d\n" % len(selected))
        for reason in sorted(rejected):
            handle.write("rejected_%s\t%d\n" % (reason, rejected[reason]))
    if len(selected) < args.min_genes:
        raise SystemExit("ERROR: retained %d orthogroups; minimum is %d" % (len(selected), args.min_genes))
    (out_dir / "COMPLETE").write_text("Orthogroup extraction completed.\n", encoding="ascii")


if __name__ == "__main__":
    main()
