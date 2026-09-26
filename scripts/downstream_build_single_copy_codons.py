#!/usr/bin/env python3
"""Build single-copy CDS/codon alignments from OrthoFinder outputs.

Inputs:
- OrthoFinder Orthogroups.tsv
- OrthoFinder MultipleSequenceAlignments directory
- canonical_manifest.tsv from downstream_prepare_inputs.py

For each orthogroup with exactly one gene in every requested species, the script
finds the protein MSA and back-translates it using the canonical primary CDS.
"""

import argparse
import csv
import re
from pathlib import Path


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path):
    order = []
    records = {}
    current = None
    desc = ""
    chunks = []
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if current is not None:
                    records[current] = (desc, "".join(chunks))
                header = line[1:].strip()
                current = header.split()[0]
                desc = header[len(current):].strip()
                order.append(current)
                chunks = []
            else:
                chunks.append(line.strip())
        if current is not None:
            records[current] = (desc, "".join(chunks))
    return order, records


def write_fasta(path, order, records):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as out:
        for seq_id in order:
            out.write(">%s\n" % seq_id)
            seq = records[seq_id]
            for idx in range(0, len(seq), 80):
                out.write(seq[idx:idx + 80] + "\n")


def key_variants(value):
    variants = {value}
    variants.add(value.replace("|", "_"))
    variants.add(value.replace("|", "."))
    variants.add(re.sub(r"[^A-Za-z0-9_.:-]", "_", value))
    if "|" in value:
        variants.add(value.split("|", 1)[0])
        variants.add(value.split("|", 1)[1])
    return variants


def load_cds_indexes(canonical_manifest, species):
    by_species = {}
    for row in read_tsv(canonical_manifest):
        label = row["label"]
        if label not in species:
            continue
        cds_path = Path(row["canonical_dir"]) / "primary.cds.fa"
        if not cds_path.is_file():
            raise SystemExit("ERROR: missing primary CDS for %s: %s" % (label, cds_path))
        _, records = read_fasta(cds_path)
        index = {}
        for seq_id, (_, seq) in records.items():
            for key in key_variants(seq_id):
                index.setdefault(key, (seq_id, seq.upper()))
        by_species[label] = index
    missing = [sp for sp in species if sp not in by_species]
    if missing:
        raise SystemExit("ERROR: species absent from canonical manifest: %s" % ",".join(missing))
    return by_species


def parse_gene_cell(value):
    value = value.strip()
    if not value:
        return []
    genes = []
    for part in value.replace(",", " ").split():
        part = part.strip()
        if part:
            genes.append(part)
    return genes


def find_cds(index, gene_id):
    for key in key_variants(gene_id):
        if key in index:
            return index[key]
    # OrthoFinder alignments may contain species-prefixed identifiers.
    for key in key_variants(gene_id.split("_", 1)[-1]):
        if key in index:
            return index[key]
    return None


def map_alignment_id_to_species(seq_id, species, single_copy_genes):
    for sp in species:
        expected = single_copy_genes[sp]
        if seq_id == expected or seq_id in key_variants(expected):
            return sp
        if seq_id.startswith(sp + "_") or seq_id.startswith(sp + "|"):
            return sp
    matches = [sp for sp in species if single_copy_genes[sp] in seq_id or seq_id in key_variants(single_copy_genes[sp])]
    return matches[0] if len(matches) == 1 else None


def backtranslate(protein_aln, cds_seq):
    cds_seq = cds_seq.upper().replace("-", "")
    codon_pos = 0
    out = []
    for aa in protein_aln:
        if aa == "-":
            out.append("---")
        else:
            codon = cds_seq[codon_pos:codon_pos + 3]
            if len(codon) != 3:
                raise ValueError("CDS exhausted while back-translating")
            out.append(codon)
            codon_pos += 3
    return "".join(out)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--canonical-manifest", required=True)
    parser.add_argument("--orthogroups-tsv", required=True)
    parser.add_argument("--protein-alignment-dir", required=True)
    parser.add_argument("--species", required=True, help="Comma-separated species labels in desired alignment order.")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--min-aa-aligned", type=int, default=60)
    args = parser.parse_args()

    species = [x.strip() for x in args.species.split(",") if x.strip()]
    out_dir = Path(args.out_dir)
    codon_dir = out_dir / "codon_alignments"
    protein_dir = out_dir / "protein_alignments"
    codon_dir.mkdir(parents=True, exist_ok=True)
    protein_dir.mkdir(parents=True, exist_ok=True)

    cds_indexes = load_cds_indexes(args.canonical_manifest, species)
    og_rows = read_tsv(args.orthogroups_tsv)
    aln_dir = Path(args.protein_alignment_dir)

    stats = []
    gene_tree_taxa = []
    used = 0
    skipped = 0

    for row in og_rows:
        og = row.get("Orthogroup") or row.get("HOG") or row.get("OG") or ""
        if not og:
            continue
        single = {}
        ok = True
        for sp in species:
            genes = parse_gene_cell(row.get(sp, ""))
            if len(genes) != 1:
                ok = False
                break
            single[sp] = genes[0]
        if not ok:
            skipped += 1
            stats.append({"orthogroup": og, "status": "skip_not_single_copy_all_species"})
            continue
        aln_path = aln_dir / ("%s.fa" % og)
        if not aln_path.is_file():
            aln_path = aln_dir / ("%s.fasta" % og)
        if not aln_path.is_file():
            skipped += 1
            stats.append({"orthogroup": og, "status": "skip_missing_protein_msa"})
            continue

        prot_order, prot_records = read_fasta(aln_path)
        codon_records = {}
        protein_records = {}
        try:
            for seq_id in prot_order:
                sp = map_alignment_id_to_species(seq_id, species, single)
                if not sp:
                    continue
                found = find_cds(cds_indexes[sp], single[sp])
                if not found:
                    raise ValueError("missing_cds_for_%s_%s" % (sp, single[sp]))
                _, cds_seq = found
                _, prot_aln = prot_records[seq_id]
                ungapped_aa = len(prot_aln.replace("-", ""))
                if ungapped_aa < args.min_aa_aligned:
                    raise ValueError("too_short_%s" % sp)
                codon_records[sp] = backtranslate(prot_aln.upper(), cds_seq)
                protein_records[sp] = prot_aln.upper()
        except ValueError as exc:
            skipped += 1
            stats.append({"orthogroup": og, "status": "skip_%s" % str(exc)})
            continue
        if set(codon_records) != set(species):
            skipped += 1
            stats.append({"orthogroup": og, "status": "skip_taxon_mismatch"})
            continue
        write_fasta(codon_dir / ("%s.fa" % og), species, codon_records)
        write_fasta(protein_dir / ("%s.fa" % og), species, protein_records)
        gene_tree_taxa.append(og)
        used += 1
        stats.append({"orthogroup": og, "status": "used"})

    with (out_dir / "single_copy_build_stats.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["orthogroup", "status"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(stats)
    with (out_dir / "single_copy_orthogroups.txt").open("w", encoding="utf-8") as handle:
        for og in gene_tree_taxa:
            handle.write(og + "\n")
    print("used_single_copy=%d skipped=%d" % (used, skipped))


if __name__ == "__main__":
    main()
