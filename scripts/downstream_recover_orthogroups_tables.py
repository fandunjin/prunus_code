#!/usr/bin/env python3
"""Recover OrthoFinder orthogroup tables from Orthogroups.txt.

This is a rescue utility for runs stopped after OrthoFinder inferred
orthogroups but before it finished the expensive MSA/gene-tree post-processing.
It uses the prepared OrthoFinder protein FASTA files and, when available,
OrthoFinder's WorkingDirectory/SequenceIDs.txt and SpeciesIDs.txt to map each
gene ID back to its species label and writes the tabular outputs required by
downstream scripts.
"""

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta_ids(path):
    ids = []
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                header = line[1:].strip()
                seq_id = header.split()[0]
                desc = header[len(seq_id):].strip()
                ids.append((seq_id, desc))
    return ids


def key_variants(value):
    variants = {value}
    variants.add(value.replace("|", "_"))
    variants.add(value.replace("|", "."))
    variants.add(re.sub(r"[^A-Za-z0-9_.:-]", "_", value))
    for suffix in (".p", ".protein", ".cds"):
        if value.endswith(suffix):
            variants.add(value[: -len(suffix)])
    if ".v" in value:
        variants.add(value.rsplit(".v", 1)[0])
    if "|" in value:
        variants.add(value.split("|", 1)[0])
        variants.add(value.split("|", 1)[1])
    return variants


def candidate_ids(seq_id, desc):
    ids = {seq_id}
    text = "%s %s" % (seq_id, desc)
    for pattern in [
        r"(?:transcript|transcript_id|protein_id|polypeptide|gene|gene_id|locus|locus_tag|Parent|ID)=([^;\s]+)",
        r"\[protein_id=([^\]]+)\]",
        r"\[gene=([^\]]+)\]",
        r"\[locus_tag=([^\]]+)\]",
    ]:
        for hit in re.findall(pattern, text):
            ids.add(hit.split()[0].strip())
    if "." in seq_id:
        ids.add(seq_id.rsplit(".", 1)[0])
    return ids


def parse_id_mapping(path):
    mapping = {}
    path = Path(path)
    if not path.is_file():
        return mapping
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if ":" in line:
                key, value = line.split(":", 1)
            else:
                parts = line.split(None, 1)
                if len(parts) != 2:
                    continue
                key, value = parts
            mapping[key.strip()] = value.strip()
    return mapping


def strip_common_suffixes(name):
    path = Path(name)
    candidates = {name, path.name, path.stem}
    current = path.name
    changed = True
    while changed:
        changed = False
        for suffix in (
            ".pep.fa", ".protein.fa", ".faa", ".fasta", ".fa", ".pep", ".protein",
            ".primary.pep.fa", ".primary.fa",
        ):
            if current.endswith(suffix):
                current = current[: -len(suffix)]
                candidates.add(current)
                changed = True
    return candidates


def load_gene_species_map(pep_manifest, working_directory=None):
    species = []
    gene_to_species = {}
    gene_alias = {}
    duplicates = set()
    label_lookup = {}
    for row in read_tsv(pep_manifest):
        label = row["label"]
        species.append(label)
        pep = Path(row["orthofinder_file"])
        for key in {label, *strip_common_suffixes(str(pep)), *strip_common_suffixes(pep.name)}:
            label_lookup[key] = label
        if not pep.is_file():
            raise SystemExit("ERROR: missing OrthoFinder protein file for %s: %s" % (label, pep))
        for seq_id, desc in read_fasta_ids(pep):
            for candidate in candidate_ids(seq_id, desc):
                for key in key_variants(candidate):
                    if key in gene_to_species and gene_to_species[key] != label:
                        duplicates.add(key)
                    gene_to_species.setdefault(key, label)
                    gene_alias.setdefault(key, seq_id)

    if working_directory:
        workdir = Path(working_directory)
        sequence_ids = parse_id_mapping(workdir / "SequenceIDs.txt")
        species_ids = parse_id_mapping(workdir / "SpeciesIDs.txt")
        internal_species = {}
        for internal_id, raw_name in species_ids.items():
            label = None
            for candidate in strip_common_suffixes(raw_name):
                if candidate in label_lookup:
                    label = label_lookup[candidate]
                    break
            if label:
                internal_species[internal_id] = label
        for internal_id, original_id in sequence_ids.items():
            sp = internal_species.get(internal_id.split("_", 1)[0])
            if not sp:
                sp = find_species_by_manifest(original_id, gene_to_species)
            if not sp:
                continue
            for key in key_variants(internal_id):
                gene_to_species.setdefault(key, sp)
                gene_alias.setdefault(key, original_id)
            for key in key_variants(original_id):
                gene_to_species.setdefault(key, sp)
                gene_alias.setdefault(key, original_id)
    for key in duplicates:
        gene_to_species.pop(key, None)
        gene_alias.pop(key, None)
    return species, gene_to_species, gene_alias


def parse_orthogroups_txt(path):
    rows = []
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line or ":" not in line:
                continue
            og, members = line.split(":", 1)
            genes = [part.strip() for part in members.replace(",", " ").split() if part.strip()]
            rows.append((og.strip(), genes))
    return rows


def find_species(gene, gene_to_species):
    for key in key_variants(gene):
        if key in gene_to_species:
            return gene_to_species[key]
    for key in key_variants(gene.split("_", 1)[-1]):
        if key in gene_to_species:
            return gene_to_species[key]
    return None


def find_species_by_manifest(gene, gene_to_species):
    return find_species(gene, gene_to_species)


def resolve_gene_id(gene, gene_alias):
    for key in key_variants(gene):
        if key in gene_alias:
            return gene_alias[key]
    for key in key_variants(gene.split("_", 1)[-1]):
        if key in gene_alias:
            return gene_alias[key]
    return gene


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orthogroups-txt", required=True)
    parser.add_argument("--pep-manifest", required=True)
    parser.add_argument("--working-directory", help="OrthoFinder WorkingDirectory containing SequenceIDs.txt and SpeciesIDs.txt.")
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    species, gene_to_species, gene_alias = load_gene_species_map(args.pep_manifest, args.working_directory)
    og_rows = parse_orthogroups_txt(args.orthogroups_txt)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    unmapped = 0
    table_rows = []
    count_rows = []
    single_copy = []
    for og, genes in og_rows:
        by_species = defaultdict(list)
        for gene in genes:
            sp = find_species(gene, gene_to_species)
            if sp is None:
                unmapped += 1
                continue
            by_species[sp].append(resolve_gene_id(gene, gene_alias))

        table = {"Orthogroup": og}
        counts = {"Orthogroup": og}
        total = 0
        is_single = True
        for sp in species:
            genes_for_sp = by_species.get(sp, [])
            table[sp] = ", ".join(genes_for_sp)
            counts[sp] = str(len(genes_for_sp))
            total += len(genes_for_sp)
            if len(genes_for_sp) != 1:
                is_single = False
        counts["Total"] = str(total)
        table_rows.append(table)
        count_rows.append(counts)
        if is_single:
            single_copy.append(og)

    with (out_dir / "Orthogroups.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["Orthogroup"] + species, lineterminator="\n")
        writer.writeheader()
        writer.writerows(table_rows)
    with (out_dir / "Orthogroups.GeneCount.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["Orthogroup"] + species + ["Total"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(count_rows)
    with (out_dir / "Orthogroups_SingleCopyOrthologues.txt").open("w", encoding="utf-8") as handle:
        for og in single_copy:
            handle.write(og + "\n")

    print("orthogroups=%d single_copy=%d unmapped_genes=%d out_dir=%s" % (len(table_rows), len(single_copy), unmapped, out_dir))


if __name__ == "__main__":
    main()
