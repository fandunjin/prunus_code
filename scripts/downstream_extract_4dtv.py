#!/usr/bin/env python3
"""Extract fourfold-degenerate sites and pairwise 4DTV from codon alignments."""

import argparse
import csv
from pathlib import Path


FOURFOLD_PREFIXES = {
    "GC",  # Ala
    "GG",  # Gly
    "CC",  # Pro
    "AC",  # Thr
    "GT",  # Val
    "TC",  # Ser
    "CG",  # Arg
    "CT",  # Leu
}
VALID_BASES = set("ACGT")


def read_fasta(path):
    records = {}
    order = []
    current = None
    chunks = []
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if current is not None:
                    records[current] = "".join(chunks).upper()
                current = line[1:].strip().split()[0]
                order.append(current)
                chunks = []
            else:
                chunks.append(line.strip())
        if current is not None:
            records[current] = "".join(chunks).upper()
    return order, records


def write_fasta(path, order, records):
    with Path(path).open("w", encoding="utf-8") as out:
        for seq_id in order:
            seq = records.get(seq_id, "")
            out.write(">%s\n" % seq_id)
            for idx in range(0, len(seq), 80):
                out.write(seq[idx:idx + 80] + "\n")


def is_fourfold_column(codons):
    for codon in codons:
        if len(codon) != 3:
            return False
        if "-" in codon or "N" in codon:
            return False
        if any(base not in VALID_BASES for base in codon):
            return False
        if codon[:2] not in FOURFOLD_PREFIXES:
            return False
    return True


def extract_4d_sites(records):
    lengths = {len(seq) for seq in records.values()}
    if len(lengths) != 1:
        raise ValueError("alignment has unequal sequence lengths")
    length = lengths.pop()
    if length % 3 != 0:
        raise ValueError("alignment length is not divisible by three")
    out = {seq_id: [] for seq_id in records}
    n_codons = 0
    n_4d = 0
    for pos in range(0, length, 3):
        codons = [seq[pos:pos + 3] for seq in records.values()]
        n_codons += 1
        if not is_fourfold_column(codons):
            continue
        n_4d += 1
        for seq_id, seq in records.items():
            out[seq_id].append(seq[pos + 2])
    return {seq_id: "".join(bases) for seq_id, bases in out.items()}, n_codons, n_4d


def pairwise_distance(seq_a, seq_b):
    valid = 0
    diff = 0
    for a, b in zip(seq_a, seq_b):
        if a in VALID_BASES and b in VALID_BASES:
            valid += 1
            if a != b:
                diff += 1
    if valid == 0:
        return "", 0, 0
    return diff / valid, valid, diff


def write_distance_matrix(path, order, records):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["taxon_a", "taxon_b", "valid_4d_sites", "different_4d_sites", "p_distance_4dtv"])
        for i, taxon_a in enumerate(order):
            for taxon_b in order[i + 1:]:
                dist, valid, diff = pairwise_distance(records[taxon_a], records[taxon_b])
                writer.writerow([taxon_a, taxon_b, valid, diff, dist])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment-dir", required=True, help="Directory containing codon alignments in FASTA format.")
    parser.add_argument("--pattern", default="*.fa", help="Codon alignment glob pattern.")
    parser.add_argument("--out-prefix", required=True, help="Output prefix for concatenated 4D alignment and reports.")
    parser.add_argument("--min-4d-sites", type=int, default=12, help="Minimum 4D sites required per gene alignment.")
    args = parser.parse_args()

    aln_dir = Path(args.alignment_dir)
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    concat = {}
    order = []
    stats = []
    used = 0
    skipped = 0

    for aln in sorted(aln_dir.glob(args.pattern)):
        gene = aln.stem
        gene_order, records = read_fasta(aln)
        if not order:
            order = gene_order
            concat = {seq_id: [] for seq_id in order}
        if gene_order != order:
            skipped += 1
            stats.append({"gene": gene, "status": "skip_taxon_order_mismatch", "codons": "", "fourfold_sites": ""})
            continue
        try:
            fourd, codons, n_4d = extract_4d_sites(records)
        except ValueError as exc:
            skipped += 1
            stats.append({"gene": gene, "status": "skip_%s" % str(exc).replace(" ", "_"), "codons": "", "fourfold_sites": ""})
            continue
        if n_4d < args.min_4d_sites:
            skipped += 1
            stats.append({"gene": gene, "status": "skip_too_few_4d_sites", "codons": codons, "fourfold_sites": n_4d})
            continue
        used += 1
        for seq_id in order:
            concat[seq_id].append(fourd[seq_id])
        stats.append({"gene": gene, "status": "used", "codons": codons, "fourfold_sites": n_4d})

    concat_records = {seq_id: "".join(chunks) for seq_id, chunks in concat.items()}
    write_fasta(str(out_prefix) + ".4d_sites.fa", order, concat_records)
    write_distance_matrix(str(out_prefix) + ".pairwise_4dtv.tsv", order, concat_records)
    with Path(str(out_prefix) + ".gene_4d_stats.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["gene", "status", "codons", "fourfold_sites"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(stats)
    print("used_genes=%d skipped_genes=%d taxa=%d sites=%d" % (used, skipped, len(order), len(next(iter(concat_records.values()), ""))))


if __name__ == "__main__":
    main()
