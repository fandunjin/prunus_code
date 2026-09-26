#!/usr/bin/env python3
from pathlib import Path
import gzip
import sys


def count_fasta(path):
    opener = gzip.open if path.suffix == ".gz" else open
    n = 0
    total = 0
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                n += 1
            else:
                total += len(line.strip())
    return n, total


def find_one(directory, patterns, exclude=()):
    for pattern in patterns:
        hits = sorted(
            p for p in directory.glob(pattern)
            if not any(token in p.name for token in exclude)
        )
        if hits:
            return hits[0]
    return None


def main():
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data")
    labels = sorted(
        p.name for p in root.iterdir()
        if p.is_dir() and not p.name.startswith("_")
    )
    out = root / "downloaded_genome_audit.tsv"
    with out.open("w", encoding="utf-8", newline="") as handle:
        handle.write("label\tgenome_file\tgenome_sequences\tgenome_bp\tgff\tcds\tprotein\tuse_note\n")
        for label in labels:
            d = root / label
            genome = find_one(d, ["*_genomic.fna.gz", "**/*_genomic.fna"], exclude=("cds_from_genomic",))
            gff = find_one(d, ["*_genomic.gff.gz", "**/*_genomic.gff"], exclude=("cds_from_genomic",))
            cds = find_one(d, ["*_cds_from_genomic.fna.gz", "**/*_cds_from_genomic.fna"])
            protein = find_one(d, ["*_protein.faa.gz", "**/*_protein.faa"])
            if genome and genome.exists():
                try:
                    seqs, bp = count_fasta(genome)
                except Exception:
                    seqs, bp = "ERROR", "ERROR"
            else:
                seqs, bp = "", ""
            if protein and cds and gff:
                note = "comparative_genomics_ready"
            elif genome:
                note = "genome_only_for_synteny_or_self_annotation"
            else:
                note = "missing_genome"
            row = [
                label,
                str(genome) if genome else "",
                str(seqs),
                str(bp),
                str(gff) if gff else "",
                str(cds) if cds else "",
                str(protein) if protein else "",
                note,
            ]
            handle.write("\t".join(row) + "\n")
    print(out)


if __name__ == "__main__":
    main()
