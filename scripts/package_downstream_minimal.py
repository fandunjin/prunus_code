#!/usr/bin/env python3
import argparse
import csv
import hashlib
import shutil
from pathlib import Path


def ensure_dir(path):
    path.mkdir(parents=True, exist_ok=True)


def copy_file(src, dst):
    if not src.is_file() or src.stat().st_size == 0:
        raise SystemExit(f"ERROR: missing or empty file: {src}")
    ensure_dir(dst.parent)
    shutil.copy2(src, dst)


def md5(path, chunk_size=1024 * 1024):
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def count_gff(path):
    counts = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                counts[fields[2]] = counts.get(fields[2], 0) + 1
    return counts


def parse_lai(path):
    result = {"lai": "", "raw_lai": "", "identity": ""}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.strip().split()
            for i, val in enumerate(fields):
                key = val.lower().rstrip(":")
                if key == "lai" and i + 1 < len(fields):
                    result["lai"] = fields[i + 1]
                elif key in {"raw_lai", "rawlai"} and i + 1 < len(fields):
                    result["raw_lai"] = fields[i + 1]
                elif key in {"mean_identity", "identity"} and i + 1 < len(fields):
                    result["identity"] = fields[i + 1].rstrip("%")
    return result


def busco_short(path):
    if not path.is_file():
        return ""
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if line.startswith("C:"):
            return line
    return ""


def write_table(path, header, rows):
    ensure_dir(path.parent)
    with path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_manifest(root):
    rows = []
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        if path.name in {"MANIFEST.tsv", "CHECKSUMS.md5"}:
            continue
        rows.append({
            "relative_path": path.relative_to(root).as_posix(),
            "bytes": path.stat().st_size,
            "md5": md5(path),
        })
    write_table(root / "MANIFEST.tsv", ["relative_path", "bytes", "md5"], rows)
    with (root / "CHECKSUMS.md5").open("w", encoding="utf-8") as out:
        for row in rows:
            out.write(f"{row['md5']}  {row['relative_path']}\n")


def write_readme(root, prefix, report_rows):
    text = f"""# Downstream minimal package

This folder contains the minimal files needed for comparative genomics.

## Contents

- genome FASTA
- gene annotation GFF3
- CDS FASTA
- protein FASTA

## Notes

- Source assembly: full 2n=32 genome package
- Gene annotation: full EDTA-masked Helixer annotation
- This package intentionally excludes large QC evidence files.
"""
    (root / "README_DOWNSTREAM_MINIMAL.md").write_text(text, encoding="utf-8")


def main():
    args = parse_args()
    source = Path(args.source_dir)
    out = Path(args.out_dir)
    prefix = args.prefix
    ensure_dir(out)

    genome_src = source / "final" / "haplotypes" / f"{prefix}.copy1.hap_provenance.chr16.fa"
    gff_src = source / "annotation" / "genes" / "copy1_helixer_edta_masked" / "pzi_copy1.edta_masked.helixer.gff3"
    cds_src = source / "annotation" / "genes" / "copy1_helixer_edta_masked" / "postprocess" / "pzi_copy1.edta_masked.helixer.cds.fa"
    pep_src = source / "annotation" / "genes" / "copy1_helixer_edta_masked" / "postprocess" / "pzi_copy1.edta_masked.helixer.protein.fa"

    genome_dst = out / f"pzi_copy1.fa"
    gff_dst = out / f"pzi_copy1.gff"
    cds_dst = out / f"pzi_copy1.cds"
    pep_dst = out / f"pzi_copy1.pep"

    copy_file(genome_src, genome_dst)
    copy_file(gff_src, gff_dst)
    copy_file(cds_src, cds_dst)
    copy_file(pep_src, pep_dst)

    counts = count_gff(gff_dst)
    lai_file = source / "annotation" / "repeats" / "copy1_lai" / "pzi_copy1.all_ltr.out.LAI"
    busco_file = source / "annotation" / "genes" / "copy1_helixer_edta_masked" / "busco" / "pzi_copy1_edta_masked_helixer_protein_eudicots" / "run_eudicotyledons_odb12" / "short_summary.txt"
    lai = parse_lai(lai_file) if lai_file.is_file() else {"lai": "", "raw_lai": "", "identity": ""}
    busco = busco_short(busco_file)

    report_rows = [
        {"item": "genome_fasta", "path": str(genome_dst), "notes": "copy1 haplotype FASTA for downstream comparative genomics"},
        {"item": "gene_gff3", "path": str(gff_dst), "notes": f"copy1 Helixer annotation; genes={counts.get('gene', 0)}; cds={counts.get('CDS', 0)}; exon={counts.get('exon', 0)}"},
        {"item": "cds_fasta", "path": str(cds_dst), "notes": f"records={sum(1 for _ in cds_dst.open()) if cds_dst.is_file() else ''}"},
        {"item": "protein_fasta", "path": str(pep_dst), "notes": f"records={sum(1 for _ in pep_dst.open()) if pep_dst.is_file() else ''}"},
        {"item": "lai", "path": str(lai_file) if lai_file.is_file() else "", "notes": f"LAI={lai.get('lai','')}; raw_LAI={lai.get('raw_lai','')}; identity={lai.get('identity','')}"},
        {"item": "busco", "path": str(busco_file) if busco_file.is_file() else "", "notes": busco},
    ]
    write_table(out / f"{prefix}.downstream_minimal_report.tsv", ["item", "path", "notes"], report_rows)
    write_readme(out, prefix, report_rows)
    write_manifest(out)
    print(f"Created minimal package at: {out}")


def parse_args():
    p = argparse.ArgumentParser(description="Build a minimal downstream package with only FASTA/GFF/CDS/PEP.")
    p.add_argument("--source-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--prefix", required=True)
    return p.parse_args()


if __name__ == "__main__":
    main()
