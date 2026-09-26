#!/usr/bin/env python3
import argparse
import csv
import hashlib
import shutil
from pathlib import Path


HAPS = ["hapA1", "hapA2", "hapB1", "hapB2"]


def require_file(path, label):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit("ERROR: missing %s: %s" % (label, path))
    return path


def count_fasta(path):
    count = 0
    bp = 0
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
            else:
                bp += len(line.strip())
    return count, bp


def count_gff_features(path):
    counts = {}
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                counts[fields[2]] = counts.get(fields[2], 0) + 1
    return counts


def md5sum(path):
    h = hashlib.md5()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_tsv(path, rows, fieldnames):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def copy_file(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(str(src), str(dst))
    return dst


def optional_copy(src, dst, manifest, copied, file_type, description):
    src = Path(src)
    if not src.is_file() or src.stat().st_size == 0:
        return
    dst = copy_file(src, dst)
    copied.append(dst)
    manifest.append({
        "hap_set": "all",
        "file_type": file_type,
        "path": str(dst),
        "source": str(src),
        "bytes": str(dst.stat().st_size),
        "description": description,
    })


def parse_args():
    parser = argparse.ArgumentParser(description="Collect final phased haplotype FASTA/GFF/CDS/PEP files.")
    parser.add_argument("--root-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    root = Path(args.root_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest = []
    summary = []
    copied = []
    for hap in HAPS:
        fasta = require_file(root / "haploid" / "final_haplotypes" / ("pzi_%s.fa" % hap), "%s FASTA" % hap)
        post_dir = root / "annotation" / "genes" / ("final_%s_helixer" % hap) / "postprocess"
        prefix = "pzi_%s.helixer" % hap
        gff = require_file(post_dir / ("%s.genes.gff3" % prefix), "%s GFF3" % hap)
        cds = require_file(post_dir / ("%s.cds.fa" % prefix), "%s CDS FASTA" % hap)
        pep = require_file(post_dir / ("%s.protein.fa" % prefix), "%s protein FASTA" % hap)

        targets = {
            "fa": copy_file(fasta, out_dir / ("pzi_%s.fa" % hap)),
            "gff": copy_file(gff, out_dir / ("pzi_%s.gff3" % hap)),
            "cds": copy_file(cds, out_dir / ("pzi_%s.cds.fa" % hap)),
            "pep": copy_file(pep, out_dir / ("pzi_%s.pep.fa" % hap)),
        }
        for file_type, path in targets.items():
            copied.append(path)
            manifest.append({
                "hap_set": hap,
                "file_type": file_type,
                "path": str(path),
                "source": str({"fa": fasta, "gff": gff, "cds": cds, "pep": pep}[file_type]),
                "bytes": str(path.stat().st_size),
                "description": "Final %s file for %s." % (file_type, hap),
            })
        fa_count, fa_bp = count_fasta(targets["fa"])
        cds_count, cds_bp = count_fasta(targets["cds"])
        pep_count, pep_aa = count_fasta(targets["pep"])
        gff_counts = count_gff_features(targets["gff"])
        summary.append({
            "hap_set": hap,
            "genome_sequences": str(fa_count),
            "genome_bp": str(fa_bp),
            "gff_gene": str(gff_counts.get("gene", 0)),
            "gff_mRNA": str(gff_counts.get("mRNA", 0)),
            "gff_CDS": str(gff_counts.get("CDS", 0)),
            "cds_records": str(cds_count),
            "cds_bp": str(cds_bp),
            "pep_records": str(pep_count),
            "pep_aa": str(pep_aa),
            "status": "pass" if fa_count == 8 and gff_counts.get("gene", 0) > 0 and cds_count > 0 and pep_count > 0 else "review",
        })

    phase_dir = root / "haploid" / "final_haplotypes"
    optional_copy(
        phase_dir / "FINAL_HAPLOTYPE_OVERVIEW.tsv",
        out_dir / "FINAL_HAPLOTYPE_OVERVIEW.tsv",
        manifest,
        copied,
        "phasing_overview",
        "Main x=8 quartet and hapA1/A2/B1/B2 assignment table.",
    )
    optional_copy(
        phase_dir / "FINAL_REVIEW_NOTES.tsv",
        out_dir / "FINAL_REVIEW_NOTES.tsv",
        manifest,
        copied,
        "phasing_review",
        "Review notes for ambiguous or validation-relevant phasing evidence.",
    )
    optional_copy(
        phase_dir / "FINAL_VALIDATION.tsv",
        out_dir / "FINAL_VALIDATION.tsv",
        manifest,
        copied,
        "phasing_validation",
        "Validation table proving 8 quartets, 32 chromosomes, and 8 chromosomes per hap set.",
    )
    optional_copy(
        phase_dir / "pzi_A1A2B1B2.rename_map.tsv",
        out_dir / "pzi_A1A2B1B2.rename_map.tsv",
        manifest,
        copied,
        "rename_map",
        "Old-to-new chromosome name map for chrXX_A1/A2/B1/B2.",
    )

    write_tsv(out_dir / "MANIFEST.tsv", manifest, ["hap_set", "file_type", "path", "source", "bytes", "description"])
    write_tsv(out_dir / "SUMMARY.tsv", summary, [
        "hap_set", "genome_sequences", "genome_bp", "gff_gene", "gff_mRNA", "gff_CDS",
        "cds_records", "cds_bp", "pep_records", "pep_aa", "status",
    ])
    with (out_dir / "CHECKSUMS.md5").open("w", encoding="utf-8") as out:
        for path in sorted(copied, key=lambda p: p.name):
            out.write("%s  %s\n" % (md5sum(path), path.name))
    (out_dir / "README_FINAL_GENOMES.md").write_text(
        "\n".join([
            "Final phased haplotype genome files",
            "====================================",
            "",
            "This directory contains the minimal downstream set for four x=8 haplotypes.",
            "Each haplotype has FASTA, GFF3, CDS FASTA, and protein FASTA.",
            "",
            "Files are named pzi_hapA1.*, pzi_hapA2.*, pzi_hapB1.*, and pzi_hapB2.*.",
            "The package also includes FINAL_HAPLOTYPE_OVERVIEW.tsv, FINAL_REVIEW_NOTES.tsv,",
            "FINAL_VALIDATION.tsv, and pzi_A1A2B1B2.rename_map.tsv so downstream analyses",
            "retain the phasing logic and chromosome renaming context.",
            "See MANIFEST.tsv for source paths and SUMMARY.tsv for basic counts.",
            "",
        ]),
        encoding="utf-8",
    )

    bad = [row["hap_set"] for row in summary if row["status"] != "pass"]
    if bad:
        raise SystemExit("ERROR: final genome collection has review/failing hap sets: %s" % ",".join(bad))
    print("Collected final haplotype genomes:", out_dir)


if __name__ == "__main__":
    main()
