#!/usr/bin/env python3
from pathlib import Path
import gzip
import shutil
import sys


def find_one(directory, patterns):
    for pattern in patterns:
        hits = sorted(directory.glob(pattern))
        if hits:
            return hits[0]
    return None


def count_records(path):
    opener = gzip.open if path.suffix == ".gz" else open
    n = 0
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                n += 1
    return n


def copy_as_plain_fasta(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.suffix == ".gz":
        with gzip.open(src, "rb") as inp, dst.open("wb") as out:
            shutil.copyfileobj(inp, out)
    else:
        shutil.copyfile(src, dst)


def main():
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
    out_dir = root / "comparative_genomics" / "orthofinder_input" / "pep"
    report = root / "comparative_genomics" / "orthofinder_input" / "pep_manifest.tsv"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []

    self_pep = root / "final_genome" / "pzi_copy1.pep"
    if self_pep.exists():
        dst = out_dir / "pzi_copy1.pep.fa"
        copy_as_plain_fasta(self_pep, dst)
        rows.append(("pzi_copy1", "Prunus zippeliana", "genome_helixer", str(self_pep), str(dst), count_records(dst)))

    data_root = root / "data"
    for d in sorted(data_root.iterdir() if data_root.exists() else []):
        if not d.is_dir() or d.name.startswith("_"):
            continue
        protein = find_one(d, ["*_protein.faa.gz", "**/*_protein.faa.gz", "*protein*.fa.gz", "*protein*.faa"])
        if not protein:
            continue
        dst = out_dir / f"{d.name}.pep.fa"
        copy_as_plain_fasta(protein, dst)
        rows.append((d.name, "", "genome_annotation", str(protein), str(dst), count_records(dst)))

    tx_dir = root / "comparative_genomics" / "transcriptome_derived_pep"
    if tx_dir.exists():
        for protein in sorted(tx_dir.glob("*.transdecoder.pep")):
            label = protein.name.replace(".transdecoder.pep", "")
            dst = out_dir / f"{label}.pep.fa"
            copy_as_plain_fasta(protein, dst)
            rows.append((label, "", "transcriptome_transdecoder", str(protein), str(dst), count_records(dst)))

    report.parent.mkdir(parents=True, exist_ok=True)
    with report.open("w", encoding="utf-8", newline="") as handle:
        handle.write("label\tspecies\tdata_type\tsource_file\torthofinder_file\tprotein_count\n")
        for row in rows:
            handle.write("\t".join(map(str, row)) + "\n")

    print(report)
    print(f"protein_fastas={len(rows)}")


if __name__ == "__main__":
    main()
