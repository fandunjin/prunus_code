#!/usr/bin/env python3
"""Package finalized PZI genome inputs for downstream comparative analyses."""

import argparse
import csv
import hashlib
import shutil
from pathlib import Path


CORE_PRIMARY_LABELS = ["A1B1", "A2B2"]
HAP_LABELS = ["hapA1", "hapA2", "hapB1", "hapB2"]
OPTIONAL_COMBO_LABELS = ["A1B2", "A2B1", "A1A2", "B1B2", "A1A2B1B2"]


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def require_file(path, label):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit("ERROR: missing %s: %s" % (label, path))
    return path


def copy_file(src, dst):
    src = Path(src)
    dst = Path(dst)
    require_file(src, "source file")
    ensure_dir(dst.parent)
    shutil.copy2(src, dst)
    return dst


def optional_copy(src, dst):
    src = Path(src)
    if src.is_file() and src.stat().st_size > 0:
        return copy_file(src, dst)
    return None


def md5sum(path):
    h = hashlib.md5()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def count_fasta(path):
    records = 0
    bp = 0
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                records += 1
            else:
                bp += len(line.strip())
    return records, bp


def count_gff(path):
    counts = {}
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                counts[fields[2]] = counts.get(fields[2], 0) + 1
    return counts


def write_tsv(path, rows, fieldnames):
    ensure_dir(Path(path).parent)
    with Path(path).open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def add_manifest(manifest, dst, src, file_type, analysis_use, relative_path=None):
    dst = Path(dst)
    rel = Path(relative_path) if relative_path is not None else dst
    manifest.append({
        "file_type": file_type,
        "relative_path": rel.as_posix(),
        "source": str(src),
        "bytes": str(dst.stat().st_size),
        "md5": md5sum(dst),
        "analysis_use": analysis_use,
    })


def copy_primary_set(primary_dir, out_dir, label, subdir, manifest, required=True):
    prefix = "pzi_%s.primary" % label
    genome_name = "pzi_%s.fa" % label if label.startswith("hap") else "%s.fa" % prefix
    files = {
        "genome_fasta": (primary_dir / genome_name, genome_name),
        "gff3": (primary_dir / ("%s.gff3" % prefix), "%s.gff3" % prefix),
        "cds_fasta": (primary_dir / ("%s.cds.fa" % prefix), "%s.cds.fa" % prefix),
        "protein_fasta": (primary_dir / ("%s.pep.fa" % prefix), "%s.pep.fa" % prefix),
    }
    copied = {}
    for file_type, (src, name) in files.items():
        if not src.is_file() or src.stat().st_size == 0:
            if required:
                require_file(src, "%s %s" % (label, file_type))
            return {}
        dst = copy_file(src, out_dir / subdir / name)
        copied[file_type] = dst
        add_manifest(
            manifest,
            dst,
            src,
            file_type,
            analysis_use_for_label(label, file_type),
            dst.relative_to(out_dir),
        )
    return copied


def analysis_use_for_label(label, file_type):
    if label in CORE_PRIMARY_LABELS:
        return "Species-level comparative genomics primary candidate; use primary protein for OrthoFinder/CAFE and GFF/CDS for synteny/KaKs."
    if label in HAP_LABELS:
        return "Subgenome/haplotype-specific analysis, A/B comparison, gene retention, and dominance tests."
    if label == "A1A2B1B2":
        return "Full 32-chromosome primary set for completeness checks only; not recommended as OrthoFinder main input."
    return "Alternative primary combination for sensitivity checks."


def summarize_copied(label, copied):
    if not copied:
        return None
    fa_records, fa_bp = count_fasta(copied["genome_fasta"])
    cds_records, cds_bp = count_fasta(copied["cds_fasta"])
    pep_records, pep_aa = count_fasta(copied["protein_fasta"])
    gff_counts = count_gff(copied["gff3"])
    return {
        "label": label,
        "genome_sequences": str(fa_records),
        "genome_bp": str(fa_bp),
        "gff_gene": str(gff_counts.get("gene", 0)),
        "gff_mRNA": str(gff_counts.get("mRNA", 0)),
        "gff_CDS": str(gff_counts.get("CDS", 0)),
        "cds_records": str(cds_records),
        "cds_bp": str(cds_bp),
        "pep_records": str(pep_records),
        "pep_aa": str(pep_aa),
    }


def main():
    args = parse_args()
    root = Path(args.root_dir)
    out_dir = Path(args.out_dir)
    primary_dir = Path(args.primary_dir) if args.primary_dir else root / "annotation" / "genes" / "primary_isoform_busco" / "primary"
    primary_summary_dir = primary_dir.parent / "summary"
    internal_summary_dir = root / "annotation" / "genes" / "internal_compare_eviann" / "summary"
    final_hap_dir = root / "final_haplotype_genomes_eviann"
    ensure_dir(out_dir)

    manifest = []
    summary = []

    for label in CORE_PRIMARY_LABELS:
        copied = copy_primary_set(primary_dir, out_dir, label, "species_representatives", manifest, required=True)
        row = summarize_copied(label, copied)
        if row:
            summary.append(row)

    for label in HAP_LABELS:
        copied = copy_primary_set(primary_dir, out_dir, label, "subgenomes", manifest, required=True)
        row = summarize_copied(label, copied)
        if row:
            summary.append(row)

    for label in OPTIONAL_COMBO_LABELS:
        copied = copy_primary_set(primary_dir, out_dir, label, "sensitivity_sets", manifest, required=False)
        row = summarize_copied(label, copied)
        if row:
            summary.append(row)

    metadata_sources = [
        (primary_summary_dir / "primary_summary.tsv", "metadata/primary_summary.tsv", "primary_summary", "Primary isoform file statistics."),
        (primary_summary_dir / "primary_busco_summary.tsv", "metadata/primary_busco_summary.tsv", "busco_summary", "BUSCO protein-mode comparison for single haplotypes and all combinations."),
        (internal_summary_dir / "hap_gene_count_summary.tsv", "metadata/internal_hap_gene_count_summary.tsv", "internal_summary", "Full EviAnn haplotype gene counts."),
        (internal_summary_dir / "AB_gene_number_difference.tsv", "metadata/internal_AB_gene_number_difference.tsv", "internal_summary", "A/B genome and gene-count differences."),
        (internal_summary_dir / "rbh_gene_retention_summary.tsv", "metadata/internal_rbh_gene_retention_summary.tsv", "internal_summary", "Protein RBH retention summary."),
        (internal_summary_dir / "synteny_summary.tsv", "metadata/internal_synteny_summary.tsv", "internal_summary", "Genome synteny summary."),
        (internal_summary_dir / "downstream_representative_recommendation.tsv", "metadata/internal_downstream_representative_recommendation.tsv", "internal_summary", "Earlier A+B representative recommendation."),
        (final_hap_dir / "FINAL_HAPLOTYPE_OVERVIEW.tsv", "metadata/FINAL_HAPLOTYPE_OVERVIEW.tsv", "phasing", "Final x=8 quartet and haplotype assignment table."),
        (final_hap_dir / "FINAL_VALIDATION.tsv", "metadata/FINAL_VALIDATION.tsv", "phasing", "Final haplotype validation summary."),
        (final_hap_dir / "FINAL_REVIEW_NOTES.tsv", "metadata/FINAL_REVIEW_NOTES.tsv", "phasing", "Final haplotype review notes."),
        (final_hap_dir / "pzi_A1A2B1B2.rename_map.tsv", "metadata/pzi_A1A2B1B2.rename_map.tsv", "phasing", "Old-to-new chromosome rename map."),
    ]
    for src, rel, file_type, use in metadata_sources:
        dst = optional_copy(src, out_dir / rel)
        if dst:
            add_manifest(manifest, dst, src, file_type, use, dst.relative_to(out_dir))

    write_tsv(out_dir / "SUMMARY.tsv", summary, [
        "label", "genome_sequences", "genome_bp", "gff_gene", "gff_mRNA",
        "gff_CDS", "cds_records", "cds_bp", "pep_records", "pep_aa",
    ])
    write_tsv(out_dir / "MANIFEST.tsv", manifest, [
        "file_type", "relative_path", "source", "bytes", "md5", "analysis_use",
    ])
    with (out_dir / "CHECKSUMS.md5").open("w", encoding="utf-8") as out:
        for row in manifest:
            out.write("%s  %s\n" % (row["md5"], row["relative_path"]))

    (out_dir / "README_DOWNSTREAM_INPUTS.md").write_text(
        "\n".join([
            "PZI downstream analysis input package",
            "=====================================",
            "",
            "Recommended species-level inputs:",
            "- species_representatives/pzi_A1B1.primary.*",
            "- species_representatives/pzi_A2B2.primary.*",
            "",
            "Use one of A1B1 or A2B2 as the main species representative and keep the other",
            "as a sensitivity/validation representative. Do not use a single hapA or hapB",
            "alone for species-level comparative genomics.",
            "",
            "Subgenome-specific inputs are in subgenomes/ and should be used for A/B",
            "retention, gene loss, synteny, dominance, and haplotype-specific analyses.",
            "",
            "The sensitivity_sets/ directory contains same-group, cross-group, and full",
            "A1+A2+B1+B2 primary sets for checks; the full 32-chromosome set is not",
            "recommended as the main OrthoFinder/CAFE input because BUSCO duplication is",
            "expected to be very high.",
            "",
        ]),
        encoding="utf-8",
    )
    print("Created downstream input package:", out_dir)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--primary-dir", default="")
    return parser.parse_args()


if __name__ == "__main__":
    main()
