#!/usr/bin/env python3
import argparse
import csv
import hashlib
import os
import shutil
from pathlib import Path


def require_file(path, label):
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit("ERROR: {} does not exist or is empty: {}".format(label, path))


def ensure_dir(path):
    path.mkdir(parents=True, exist_ok=True)


def copy_file(src, dst):
    if not src.is_file():
        return False
    ensure_dir(dst.parent)
    shutil.copy2(str(src), str(dst))
    return True


def copy_first(candidates, dst, required=False, label=None):
    for src in candidates:
        if copy_file(src, dst):
            return str(src)
    if required:
        raise SystemExit("ERROR: missing required file for {} -> {}".format(label or dst.name, dst))
    return ""


def read_tsv(path):
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fieldnames, rows):
    ensure_dir(path.parent)
    with path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def file_md5(path, chunk_size=1024 * 1024):
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def md5_manifest(final_dir):
    rows = []
    for path in sorted(p for p in final_dir.rglob("*") if p.is_file()):
        if path.name in {"CHECKSUMS.md5", "MANIFEST.tsv"}:
            continue
        rel = path.relative_to(final_dir).as_posix()
        rows.append((file_md5(path), rel, path.stat().st_size))
    with (final_dir / "CHECKSUMS.md5").open("w") as out:
        for checksum, rel, size in rows:
            out.write("{}  {}\n".format(checksum, rel))
    with (final_dir / "MANIFEST.tsv").open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["relative_path", "bytes", "md5"])
        for checksum, rel, size in rows:
            writer.writerow([rel, size, checksum])


def summarize_busco(path):
    text = path.read_text(errors="replace") if path.is_file() else ""
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("C:"):
            return line
    return "not_available"


def summarize_quast(path):
    if not path.is_file() or path.stat().st_size == 0:
        return {}
    result = {}
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header or len(header) < 2:
            return result
        for fields in reader:
            if len(fields) >= 2:
                result[fields[0]] = fields[1]
    return result


def parse_syri_summary(path):
    result = {}
    if not path.is_file():
        return result
    section = None
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line == "#Structural annotations":
                section = "structural"
                continue
            if line == "#Sequence annotations":
                section = "sequence"
                continue
            if line.startswith("#") or section is None:
                continue
            parts = line.split("\t")
            if len(parts) >= 4:
                key = "{}_{}".format(section, parts[0].lower().replace(" ", "_").replace("(", "").replace(")", ""))
                result[key + "_count"] = parts[1]
                result[key + "_length_ref"] = parts[2]
                result[key + "_length_qry"] = parts[3]
    return result


def confidence_label(cis_percent, telomere_ends, provenance_status, dominance):
    notes = []
    level = "high"
    if cis_percent is None:
        level = "review"
        notes.append("missing Hi-C contact summary")
    elif cis_percent < 55.0:
        level = "low_hic_review"
        notes.append("low Hi-C cis_percent")
    elif cis_percent < 60.0:
        level = "moderate"
        notes.append("moderate Hi-C cis_percent")

    if telomere_ends == 0:
        if level == "high":
            level = "moderate"
        notes.append("no terminal TTTAGGG signal in 100 kb windows")
    elif telomere_ends == 1:
        if level == "high":
            level = "moderate"
        notes.append("single terminal TTTAGGG signal")

    if provenance_status and provenance_status != "confident":
        level = "review"
        notes.append("hap provenance not confident")
    if dominance is not None and dominance < 0.90:
        level = "review"
        notes.append("hap provenance dominance below 0.90")

    if not notes:
        notes.append("supported by Hi-C, telomere, and hap provenance evidence")
    return level, "; ".join(notes)


def build_chromosome_confidence(source_dir, final_dir, prefix):
    seq_rows = read_tsv(source_dir / "{}.sequence_stats.tsv".format(prefix))
    hic_rows = read_tsv(source_dir / "hic_qc" / "chr_contact_support.tsv")
    tidk_summary_path = source_dir / "tidk_chr_end_summary.tsv"
    if not tidk_summary_path.is_file():
        build_tidk_end_summary(source_dir, prefix, tidk_summary_path)
    tidk_rows = read_tsv(tidk_summary_path)
    provenance_rows = read_tsv(source_dir / "haplotype_provenance" / "{}.chr_hap_provenance.tsv".format(prefix))
    if not provenance_rows:
        provenance_rows = read_tsv(source_dir / "{}.chr_hap_provenance.tsv".format(prefix))

    hic = {row.get("chr"): row for row in hic_rows}
    tidk = {row.get("chr"): row for row in tidk_rows}
    prov = {row.get("chr"): row for row in provenance_rows}

    out_rows = []
    for row in seq_rows:
        seq_id = row.get("seq_id", "")
        if not seq_id.startswith("chr"):
            continue
        cis = None
        if seq_id in hic and hic[seq_id].get("cis_percent"):
            cis = float(hic[seq_id]["cis_percent"])
        telomere_ends = None
        if seq_id in tidk and tidk[seq_id].get("telomere_ends"):
            telomere_ends = int(tidk[seq_id]["telomere_ends"])
        dominance = None
        if seq_id in prov and prov[seq_id].get("dominance"):
            dominance = float(prov[seq_id]["dominance"])
        status = prov.get(seq_id, {}).get("status", "")
        level, note = confidence_label(cis, telomere_ends, status, dominance)
        out_rows.append({
            "chr": seq_id,
            "length": row.get("length", ""),
            "n_bases": row.get("n_bases", ""),
            "gap_count": row.get("gap_count", ""),
            "hic_cis_percent": "" if cis is None else "{:.4f}".format(cis),
            "hic_total_contacts": hic.get(seq_id, {}).get("total_contacts", ""),
            "telomere_ends_100kb": "" if telomere_ends is None else str(telomere_ends),
            "hap_dominant_source": prov.get(seq_id, {}).get("dominant_hap", ""),
            "hap_dominance": "" if dominance is None else "{:.6f}".format(dominance),
            "hap_status": status,
            "final_confidence": level,
            "notes": note,
        })
    fields = [
        "chr", "length", "n_bases", "gap_count", "hic_cis_percent", "hic_total_contacts",
        "telomere_ends_100kb", "hap_dominant_source", "hap_dominance", "hap_status",
        "final_confidence", "notes",
    ]
    write_tsv(final_dir / "tables" / "{}.chromosome_confidence.tsv".format(prefix), fields, out_rows)
    return out_rows


def build_tidk_end_summary(source_dir, prefix, out_path):
    seq_rows = read_tsv(source_dir / "{}.sequence_stats.tsv".format(prefix))
    chr_lengths = {}
    for row in seq_rows:
        seq_id = row.get("seq_id", "")
        if seq_id.startswith("chr") and row.get("length"):
            chr_lengths[seq_id] = int(row["length"])

    candidates = [
        source_dir / "homolog_syri_tidk" / "tidk_pzi_2n32_final_telomeric_repeat_windows.tsv",
        source_dir / "tidk_pzi_2n32_final_telomeric_repeat_windows.tsv",
    ]
    tidk_path = None
    for candidate in candidates:
        if candidate.is_file() and candidate.stat().st_size > 0:
            tidk_path = candidate
            break
    if tidk_path is None:
        return

    sums = {chrom: {"left": 0, "right": 0} for chrom in chr_lengths}
    for row in read_tsv(tidk_path):
        chrom = row.get("id", "")
        if chrom not in chr_lengths:
            continue
        window = int(row.get("window", "0") or 0)
        repeats = int(row.get("forward_repeat_number", "0") or 0) + int(row.get("reverse_repeat_number", "0") or 0)
        if window <= 100000:
            sums[chrom]["left"] += repeats
        if window >= chr_lengths[chrom] - 100000:
            sums[chrom]["right"] += repeats

    rows = []
    for chrom in sorted(chr_lengths):
        left = sums[chrom]["left"]
        right = sums[chrom]["right"]
        rows.append({
            "chr": chrom,
            "left100kb": str(left),
            "right100kb": str(right),
            "telomere_ends": str(int(left > 20) + int(right > 20)),
        })
    write_tsv(out_path, ["chr", "left100kb", "right100kb", "telomere_ends"], rows)


def write_final_summary(source_dir, final_dir, prefix, confidence_rows):
    quast = summarize_quast(source_dir / "qc" / "quast" / "report.tsv")
    busco = summarize_busco(source_dir / "qc" / "short_summary.specific.eudicotyledons_odb12.busco_pzi_2n32_final.txt")
    low_hic = [r["chr"] for r in confidence_rows if r["final_confidence"] == "low_hic_review"]
    moderate = [r["chr"] for r in confidence_rows if r["final_confidence"] == "moderate"]

    rows = [
        {"metric": "final_version", "value": prefix},
        {"metric": "main_assembly", "value": "2n=32 chromosome-level genome plus unplaced sequences"},
        {"metric": "total_length_bp", "value": quast.get("Total length", "not_available")},
        {"metric": "num_sequences", "value": quast.get("# contigs", "not_available")},
        {"metric": "largest_sequence_bp", "value": quast.get("Largest contig", "not_available")},
        {"metric": "n50_bp", "value": quast.get("N50", "not_available")},
        {"metric": "gc_percent", "value": quast.get("GC (%)", "not_available")},
        {"metric": "busco_eudicotyledons_odb12", "value": busco},
        {"metric": "low_hic_review_chromosomes", "value": ",".join(low_hic) if low_hic else "none"},
        {"metric": "moderate_review_chromosomes", "value": ",".join(moderate) if moderate else "none"},
        {"metric": "annotation_status", "value": "not included"},
    ]
    write_tsv(final_dir / "tables" / "{}.final_genome_summary.tsv".format(prefix), ["metric", "value"], rows)


def write_readme(final_dir, prefix):
    text = """# Final genome package: {prefix}

This directory is the frozen genome package for downstream analyses.
The primary assembly is the 2n=32 chromosome-level FASTA in `genome/`.
The `haplotypes/` directory contains provenance-corrected n=16 copy FASTAs for comparative analyses only.

## Directory layout

- `genome/`: final 2n=32 genome FASTA, chromosome-only FASTA, unplaced FASTA, indexes, AGP, sizes, and checksums.
- `haplotypes/`: corrected homolog-pair table and hap1/hap2-provenance n=16 FASTAs.
- `tables/`: sequence statistics, gap table, rename map, chromosome confidence table, and final summary.
- `qc/`: copied QC evidence from QUAST, BUSCO, FCS/GX, Hi-C, TIDK, SyRI, and plotsr.
- `MANIFEST.tsv` and `CHECKSUMS.md5`: file inventory and checksums for the frozen package.

## Interpretation notes

- BUSCO duplication is expected because this package keeps both homologous chromosome copies.
- `chr27`, `chr28`, `chr29`, and `chr30` have lower Hi-C contact support than most chromosomes.
  They are retained because telomere evidence, contamination checks, and haplotype provenance support them.
- The n=16 copy FASTAs are helper datasets for comparison and plotting, not replacements for the primary 2n=32 assembly.
- Gene annotation and functional annotation are intentionally not included in this package.
""".format(prefix=prefix)
    (final_dir / "README_FINAL_GENOME.md").write_text(text)


def copy_package_files(args):
    source = Path(args.source_dir)
    final = Path(args.final_dir)
    prefix = args.prefix
    require_file(source / "{}.final.fa".format(prefix), "final FASTA")

    genome_dir = final / "genome"
    hap_dir = final / "haplotypes"
    tables_dir = final / "tables"
    qc_dir = final / "qc"
    for path in [genome_dir, hap_dir, tables_dir, qc_dir]:
        ensure_dir(path)

    genome_files = [
        "{}.final.fa", "{}.final.fa.fai", "{}.chromosomes.fa", "{}.chromosomes.fa.fai",
        "{}.unplaced.fa", "{}.unplaced.fa.fai", "{}.n16_representative_copy1.fa",
        "{}.n16_representative_copy1.fa.fai", "{}.chrom.sizes", "{}.chromosomes.sizes",
        "{}.agp", "{}.md5",
    ]
    for pattern in genome_files:
        name = pattern.format(prefix)
        copy_file(source / name, genome_dir / name)

    table_files = [
        "{}.rename_map.tsv", "{}.homolog_pairs.tsv", "{}.sequence_stats.tsv",
        "{}.gap_stats.tsv",
    ]
    for pattern in table_files:
        name = pattern.format(prefix)
        copy_file(source / name, tables_dir / name)

    hap_sources = [
        "{}.chr_hap_provenance.tsv", "{}.corrected_homolog_pairs.tsv",
        "{}.copy1.hap_provenance.chr16.fa", "{}.copy1.hap_provenance.chr16.fa.fai",
        "{}.copy2.hap_provenance.chr16.fa", "{}.copy2.hap_provenance.chr16.fa.fai",
    ]
    for pattern in hap_sources:
        name = pattern.format(prefix)
        copy_first([source / "haplotype_provenance" / name, source / name], hap_dir / name)

    copy_first([source / "homolog_syri_tidk" / "pzi.copy_fastas.map.tsv", source / "pzi.copy_fastas.map.tsv"], hap_dir / "pzi.copy_fastas.map.tsv")
    copy_first(
        [
            source / "homolog_syri_tidk" / "pzi.copy1_vs_copy2.same_chr_alignment.tsv",
            source / "pzi.copy1_vs_copy2.same_chr_alignment.tsv",
        ],
        hap_dir / "pzi.copy1_vs_copy2.same_chr_alignment.tsv",
    )
    copy_first(
        [
            source / "homolog_syri_tidk" / "review_checklist.txt",
            source / "review_checklist.txt",
        ],
        qc_dir / "syri_plotsr" / "review_checklist.txt",
    )

    copy_tree_files(source / "hic_qc", qc_dir / "hic", [".tsv", ".png"])
    copy_tree_files(source / "qc" / "quast", qc_dir / "quast", [".txt", ".tsv", ".pdf", ".html", ".log"])
    copy_tree_files(source / "qc", qc_dir / "fcs_busco", [".txt", ".rpt", ".log"], max_depth=1)

    syri_dirs = [source / "homolog_syri_tidk", source]
    for base in syri_dirs:
        if not base.is_dir():
            continue
        for path in base.iterdir():
            if not path.is_file():
                continue
            name = path.name
            if name.startswith("tidk_") and path.suffix in {".svg", ".tsv"}:
                copy_file(path, qc_dir / "tidk" / name)
            if name.startswith("ppa_vs_pzi_") and (
                name.endswith(".png") or name.endswith(".summary") or name.endswith(".orientation.tsv") or name.endswith(".syri.log")
            ):
                copy_file(path, qc_dir / "syri_plotsr" / name)
            if args.include_large_evidence and name.startswith("ppa_vs_pzi_") and (
                name.endswith(".syri.out") or name.endswith(".vcf") or name.endswith(".paf") or name.endswith(".bam") or name.endswith(".bam.bai")
            ):
                copy_file(path, qc_dir / "syri_plotsr" / "large_evidence" / name)

    confidence = build_chromosome_confidence(source, final, prefix)
    copy_file(source / "tidk_chr_end_summary.tsv", tables_dir / "tidk_chr_end_summary.tsv")
    write_final_summary(source, final, prefix, confidence)
    write_readme(final, prefix)
    md5_manifest(final)


def copy_tree_files(src_dir, dst_dir, suffixes, max_depth=None):
    if not src_dir.is_dir():
        return
    src_dir = src_dir.resolve()
    for path in src_dir.rglob("*"):
        if not path.is_file():
            continue
        rel = path.relative_to(src_dir)
        if max_depth is not None and len(rel.parts) > max_depth:
            continue
        if path.suffix in suffixes:
            copy_file(path, dst_dir / rel)


def parse_args():
    parser = argparse.ArgumentParser(description="Create a frozen final genome package for downstream analyses.")
    parser.add_argument("--source-dir", required=True)
    parser.add_argument("--final-dir", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--include-large-evidence", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    copy_package_files(args)


if __name__ == "__main__":
    main()
