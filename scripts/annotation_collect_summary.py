#!/usr/bin/env python3
import argparse
import csv
import hashlib
from pathlib import Path


def ensure_dir(path):
    path.mkdir(parents=True, exist_ok=True)


def find_first(base_dir, patterns):
    for pattern in patterns:
        for path in sorted(base_dir.rglob(pattern)):
            if path.is_file() and path.stat().st_size > 0:
                return path
    return None


def md5(path, chunk_size=1024 * 1024):
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def manifest(root):
    rows = []
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        if path.name in {"MANIFEST.tsv", "CHECKSUMS.md5"}:
            continue
        rel = path.relative_to(root).as_posix()
        rows.append((rel, path.stat().st_size, md5(path)))
    with (root / "MANIFEST.tsv").open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["relative_path", "bytes", "md5"])
        for row in rows:
            writer.writerow(row)
    with (root / "CHECKSUMS.md5").open("w") as out:
        for rel, size, checksum in rows:
            out.write("{}  {}\n".format(checksum, rel))


def count_gff_features(gff):
    counts = {}
    if not gff or not gff.is_file():
        return counts
    with gff.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                counts[fields[2]] = counts.get(fields[2], 0) + 1
    return counts


def parse_lai(lai_file):
    result = {"lai_file": "", "lai": "", "raw_lai": ""}
    if not lai_file or not lai_file.is_file():
        return result
    result["lai_file"] = str(lai_file)
    with lai_file.open(errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.strip().split()
            for i, value in enumerate(fields):
                key = value.lower().rstrip(":")
                if key == "lai" and i + 1 < len(fields):
                    result["lai"] = fields[i + 1]
                if key in {"raw_lai", "rawlai"} and i + 1 < len(fields):
                    result["raw_lai"] = fields[i + 1]
    return result


def parse_edta_sum(sum_file):
    result = {"edta_sum_file": "", "summary_lines": "0"}
    if not sum_file or not sum_file.is_file():
        return result
    result["edta_sum_file"] = str(sum_file)
    with sum_file.open(errors="replace") as handle:
        lines = [line for line in handle if line.strip()]
    result["summary_lines"] = str(len(lines))
    return result


def write_readme(annotation_dir):
    text = """# Annotation package

This package contains repeat annotation, LAI quality estimates, and Helixer gene prediction prepared from the frozen genome package.

Recommended downstream comparative-genomics inputs are listed in `summary/downstream_inputs.tsv`.

Interpretation:

- `copy1` is the default n=16 haplotype-provenance genome for cross-species comparative genomics.
- `chr32` EDTA/LAI results describe the frozen 2n=32 chromosome-level assembly and should be used as assembly/repeat-quality evidence.
- `copy2` Helixer output is optional and should be used for haplotype comparison, not merged into the cross-species primary gene set.
"""
    (annotation_dir / "README_ANNOTATION.md").write_text(text)


def parse_args():
    parser = argparse.ArgumentParser(description="Collect annotation outputs and write summary tables.")
    parser.add_argument("--root-dir", required=True)
    parser.add_argument("--annotation-dir", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    root = Path(args.root_dir)
    annotation = Path(args.annotation_dir)
    summary = annotation / "summary"
    ensure_dir(summary)

    copy1_edta = annotation / "repeats" / "copy1_edta"
    chr32_edta = annotation / "repeats" / "chr32_edta"
    copy1_lai = annotation / "repeats" / "copy1_lai"
    chr32_lai = annotation / "repeats" / "chr32_lai"
    copy1_helixer = annotation / "genes" / "copy1_helixer_edta_masked"
    copy2_helixer = annotation / "genes" / "copy2_helixer_edta_masked_optional"
    copy1_helixer_unmasked = annotation / "genes" / "copy1_helixer"
    copy2_helixer_unmasked = annotation / "genes" / "copy2_helixer_optional"

    repeat_rows = []
    for label, edta_dir, lai_dir in [
        ("copy1_n16", copy1_edta, copy1_lai),
        ("chr32_2n", chr32_edta, chr32_lai),
    ]:
        sum_file = find_first(edta_dir, ["*.EDTA.TEanno.sum", "*.TEanno.sum"])
        gff = find_first(edta_dir, ["*.EDTA.TEanno.gff3"])
        lai_file = find_first(lai_dir, ["*.out.LAI"])
        lai = parse_lai(lai_file)
        edta = parse_edta_sum(sum_file)
        repeat_rows.append({
            "dataset": label,
            "te_gff3": "" if gff is None else str(gff),
            "edta_sum": edta["edta_sum_file"],
            "lai_file": lai["lai_file"],
            "lai": lai["lai"],
            "raw_lai": lai["raw_lai"],
        })
    with (summary / "repeat_lai_summary.tsv").open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=["dataset", "te_gff3", "edta_sum", "lai_file", "lai", "raw_lai"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(repeat_rows)

    gene_rows = []
    for label, gene_dir in [
        ("copy1_n16_edta_masked_main", copy1_helixer),
        ("copy2_n16_edta_masked_optional", copy2_helixer),
        ("copy1_n16_unmasked_previous", copy1_helixer_unmasked),
        ("copy2_n16_unmasked_previous", copy2_helixer_unmasked),
    ]:
        gff = find_first(gene_dir, ["*.gff3"])
        counts = count_gff_features(gff)
        gene_rows.append({
            "dataset": label,
            "gff3": "" if gff is None else str(gff),
            "gene_count": counts.get("gene", 0),
            "mrna_count": counts.get("mRNA", counts.get("transcript", 0)),
            "cds_count": counts.get("CDS", 0),
            "exon_count": counts.get("exon", 0),
        })
    with (summary / "annotation_summary.tsv").open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=["dataset", "gff3", "gene_count", "mrna_count", "cds_count", "exon_count"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(gene_rows)

    downstream_rows = [
        {
            "role": "comparative_genomics_main_genome",
            "path": str(root / "final_genome" / "pzi_copy1.fa"),
            "notes": "default n=16 genome for cross-species analyses",
        },
        {
            "role": "comparative_genomics_main_gene_gff3",
            "path": str(copy1_helixer / "pzi_copy1.edta_masked.helixer.gff3"),
            "notes": "primary EDTA-masked Helixer gene annotation",
        },
        {
            "role": "comparative_genomics_main_repeat_gff3",
            "path": "" if repeat_rows[0]["te_gff3"] == "" else repeat_rows[0]["te_gff3"],
            "notes": "copy1 repeat annotation for TE comparison",
        },
        {
            "role": "assembly_quality_repeat_evidence",
            "path": "" if repeat_rows[1]["te_gff3"] == "" else repeat_rows[1]["te_gff3"],
            "notes": "2n=32 chromosome repeat annotation and LAI evidence",
        },
    ]
    with (summary / "downstream_inputs.tsv").open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=["role", "path", "notes"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(downstream_rows)

    write_readme(annotation)
    manifest(annotation)


if __name__ == "__main__":
    main()
