#!/usr/bin/env python3
"""Collect assembly, annotation, ncRNA, functional, EDTA and LAI metrics.

This is a read-only report step.  It never reruns an annotation program.  The
report distinguishes true haplotype-level evidence from a complete-species
file that is merely shared by all haplotype rows.
"""

import argparse
import csv
import gzip
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable, Iterator, Optional


HAPLOTYPES = {
    "pzi": ("hapA1", "hapA2", "hapB1", "hapB2"),
    "psp": ("hapA1", "hapA2", "hapB1", "hapB2"),
    "ppa": ("hapA", "hapB"),
}


def open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def nonempty(path):
    return bool(path and path.is_file() and path.stat().st_size > 0)


def first_existing(paths):
    for path in paths:
        if nonempty(path):
            return path
    return None


def first_recursive(directories, patterns):
    seen = set()
    for directory in directories:
        if not directory.is_dir():
            continue
        for pattern in patterns:
            for path in sorted(directory.rglob(pattern)):
                resolved = path.resolve()
                if resolved in seen:
                    continue
                seen.add(resolved)
                if nonempty(path):
                    return path
    return None


def nx(lengths, fraction):
    if not lengths:
        return 0
    target = sum(lengths) * fraction
    cumulative = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= target:
            return length
    return 0


def fasta_stats(path):
    result = {
        "genome_bp": "NA", "scaffold_count": "NA", "contig_count": "NA",
        "scaffold_n50": "NA", "scaffold_n90": "NA", "contig_n50": "NA",
        "contig_n90": "NA", "gc_pct": "NA", "gap_bp": "NA",
    }
    if not nonempty(path):
        return result

    scaffolds = []
    contigs = []
    scaffold_length = 0
    contig_length = 0
    gc = 0
    acgt = 0
    gap_bp = 0
    in_record = False

    def finish_record() -> None:
        nonlocal scaffold_length, contig_length
        if not in_record:
            return
        if contig_length:
            contigs.append(contig_length)
        scaffolds.append(scaffold_length)
        scaffold_length = 0
        contig_length = 0

    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish_record()
                in_record = True
                continue
            if not in_record:
                continue
            sequence = line.upper()
            scaffold_length += len(sequence)
            for base in sequence:
                if base == "N":
                    gap_bp += 1
                    if contig_length:
                        contigs.append(contig_length)
                        contig_length = 0
                else:
                    contig_length += 1
                    if base in "GC":
                        gc += 1
                    if base in "ACGT":
                        acgt += 1
    finish_record()

    result.update({
        "genome_bp": sum(scaffolds),
        "scaffold_count": len(scaffolds),
        "contig_count": len(contigs),
        "scaffold_n50": nx(scaffolds, 0.50),
        "scaffold_n90": nx(scaffolds, 0.90),
        "contig_n50": nx(contigs, 0.50),
        "contig_n90": nx(contigs, 0.90),
        "gc_pct": round(100.0 * gc / acgt, 3) if acgt else "NA",
        "gap_bp": gap_bp,
    })
    return result


def fasta_record_count(path):
    if not nonempty(path):
        return "NA", "NA"
    count = 0
    total = 0
    with open_text(path) as handle:
        for raw in handle:
            if raw.startswith(">"):
                count += 1
            elif count:
                total += len(raw.strip())
    return count, total


def canonical_feature(feature: str) -> str:
    key = feature.lower()
    if key in {"mrna", "transcript", "rna"}:
        return "mRNA"
    if key == "gene":
        return "gene"
    if key == "cds":
        return "CDS"
    if key == "exon":
        return "exon"
    return feature


def gff_stats(path):
    result = {
        "gene_count": "NA", "mrna_count": "NA", "cds_count": "NA",
        "exon_count": "NA", "cds_bp": "NA", "exon_bp": "NA",
        "gff_feature_count": "NA", "gff_bad_rows": "NA",
    }
    if not nonempty(path):
        return result
    counts = Counter()
    feature_bp = Counter()
    bad_rows = 0
    with open_text(path) as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                bad_rows += 1
                continue
            feature = canonical_feature(fields[2])
            counts[feature] += 1
            try:
                start, end = int(fields[3]), int(fields[4])
                if start >= 1 and end >= start:
                    feature_bp[feature] += end - start + 1
            except ValueError:
                bad_rows += 1
    result.update({
        "gene_count": counts.get("gene", 0),
        "mrna_count": counts.get("mRNA", 0),
        "cds_count": counts.get("CDS", 0),
        "exon_count": counts.get("exon", 0),
        "cds_bp": feature_bp.get("CDS", 0),
        "exon_bp": feature_bp.get("exon", 0),
        "gff_feature_count": sum(counts.values()),
        "gff_bad_rows": bad_rows,
    })
    return result


def parse_edta_sum(path):
    result = {
        "edta_genome_bp": "NA", "edta_repeat_bp": "NA",
        "edta_repeat_pct": "NA", "edta_summary_lines": "NA",
    }
    if not nonempty(path):
        return result
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    result["edta_summary_lines"] = sum(bool(line.strip()) for line in lines)
    for line in lines:
        match = re.search(r"Total Length:\s*([\d,]+)\s*bp", line, re.I)
        if match:
            result["edta_genome_bp"] = int(match.group(1).replace(",", ""))
        match = re.search(
            r"(?:total interspersed|^\s*Total)\s+\d+\s+([\d,]+)\s+([\d.]+)%",
            line,
            re.I,
        )
        if match:
            result["edta_repeat_bp"] = int(match.group(1).replace(",", ""))
            result["edta_repeat_pct"] = float(match.group(2))
    return result


def parse_lai(path):
    if not nonempty(path):
        return "NA"
    with open_text(path) as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if fields and fields[0].lower() == "whole_genome":
                if len(fields) > 6:
                    return fields[6]
                break
    return "NA"


def parse_busco(path):
    result = {
        "busco_complete_pct": "NA", "busco_single_copy_pct": "NA",
        "busco_duplicated_pct": "NA", "busco_fragmented_pct": "NA",
        "busco_missing_pct": "NA",
    }
    if not nonempty(path):
        return result
    text = path.read_text(encoding="utf-8", errors="replace")
    match = re.search(
        r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%",
        text,
    )
    if match:
        result.update({
            "busco_complete_pct": float(match.group(1)),
            "busco_single_copy_pct": float(match.group(2)),
            "busco_duplicated_pct": float(match.group(3)),
            "busco_fragmented_pct": float(match.group(4)),
            "busco_missing_pct": float(match.group(5)),
        })
    return result


def parse_key_value_tsv(path):
    values = {}
    if not nonempty(path):
        return values
    with open_text(path) as handle:
        for index, raw in enumerate(handle):
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if index == 0 and len(fields) >= 2 and fields[0].lower() in {"metric", "key"}:
                continue
            if len(fields) >= 2:
                values[fields[0]] = fields[1]
    return values


def count_annotation_table(path, go_table=False):
    result = {"rows": "NA", "unique_proteins": "NA", "go_proteins": "NA"}
    if not nonempty(path):
        return result
    rows = 0
    identifiers = set()
    go_ids = set()
    with open_text(path) as handle:
        header = None
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if header is None and fields[0].lower() in {"source", "protein_id", "qseqid"}:
                header = [item.lower() for item in fields]
                continue
            rows += 1
            id_index = 0
            if header and "qseqid" in header:
                id_index = header.index("qseqid")
            elif header and "protein_id" in header:
                id_index = header.index("protein_id")
            if id_index < len(fields) and fields[id_index] not in {"", "-"}:
                identifiers.add(fields[id_index])
            if go_table and len(fields) > 1 and fields[1] not in {"", "-"}:
                go_ids.add(fields[0])
    result["rows"] = rows
    result["unique_proteins"] = len(identifiers)
    result["go_proteins"] = len(go_ids)
    return result


def path_scope(directory, hap):
    """Mark a file as haplotype-level only when its path names this hap."""
    parts = {part.lower() for part in directory.parts}
    return "haplotype" if hap.lower() in parts else "species_complete"


def parse_ncrna(directory, hap):
    result = {
        "source": "", "scope": "missing", "total_features": "NA",
        "trna_count": "NA", "rfam_count": "NA",
    }
    if not directory or not directory.is_dir():
        return result
    combined = first_existing([
        directory / "ncrna.combined.gff3",
        directory / "combined.gff3",
    ]) or first_recursive([directory], ["*.ncrna.combined.gff3", "*combined.gff3"])
    trna = first_existing([directory / "tRNAscan-SE.gff3"])
    trna = trna or first_recursive([directory], ["*.tRNAscan-SE.gff3"])
    rfam = first_existing([directory / "rfam.cmscan.gff3"])
    rfam = rfam or first_recursive([directory], ["*.rfam.cmscan.gff3"])
    summary = first_existing([directory / "ncrna_summary.tsv"])
    summary = summary or first_recursive([directory], ["*.ncrna_summary.tsv"])
    if combined:
        files = [(combined, "combined")]
    else:
        files = [(item, label) for item, label in ((trna, "trna"), (rfam, "rfam")) if item]
    if files:
        total = 0
        trna_count = 0
        rfam_count = 0
        for path, label in files:
            with open_text(path) as handle:
                for raw in handle:
                    if not raw.strip() or raw.startswith("#"):
                        continue
                    fields = raw.rstrip("\n").split("\t")
                    if len(fields) != 9:
                        continue
                    total += 1
                    feature = fields[2].lower()
                    source = fields[1].lower()
                    if "trna" in feature or "trna" in source or label == "trna":
                        trna_count += 1
                    if label == "rfam" or any(token in source for token in ("rfam", "infernal", "cmscan")):
                        rfam_count += 1
        result.update({"source": str(files[0][0]), "total_features": total,
                       "trna_count": trna_count, "rfam_count": rfam_count})
    elif summary:
        counts = {}
        with open_text(summary) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                counts[row.get("class", "")] = row.get("count", "NA")
        result.update({
            "source": str(summary),
            "total_features": counts.get("total", "NA"),
            "trna_count": counts.get("tRNA", "NA"),
            "rfam_count": counts.get("Rfam_cm_hit", "NA"),
        })
        if result["total_features"] == "NA":
            try:
                result["total_features"] = sum(int(value) for key, value in counts.items()
                                               if key in {"tRNA", "Rfam_cm_hit"})
            except (TypeError, ValueError):
                pass
    else:
        return result
    result["scope"] = path_scope(directory, hap)
    return result


def locate_ncrna_dir(root, species, hap, ncrna_root):
    candidates = [
        ncrna_root / species / hap,
        root / "annotation" / "repartition" / "ncrna" / species / hap,
        root / "annotation" / "ncrna" / species / hap,
        root / "annotation" / "ncrna" / species,
        root / "data" / species / "haploid" / "annotation" / "ncrna" / hap,
    ]
    for directory in candidates:
        if directory.is_dir() and first_recursive([directory], ["ncrna.combined.gff3", "ncrna_summary.tsv", "*.gff3"]):
            return directory
    return None


def parse_functional(directory, hap):
    result = {
        "source": "", "scope": "missing", "hit_rows": "NA",
        "annotated_proteins": "NA", "go_proteins": "NA",
    }
    if not directory or not directory.is_dir():
        return result
    protein_table = first_existing([
        directory / "protein_annotation.tsv",
        directory / "swissprot.diamond.tsv",
    ]) or first_recursive([directory], ["*.protein_annotation.tsv", "*.swissprot.diamond.tsv"])
    go_table = first_existing([directory / "go_annotation.tsv"])
    go_table = go_table or first_recursive([directory], ["*.go_annotation.tsv"])
    if not protein_table and not go_table:
        return result
    protein = count_annotation_table(protein_table)
    go = count_annotation_table(go_table, go_table=True)
    result.update({
        "source": str(protein_table or go_table),
        "hit_rows": protein["rows"],
        "annotated_proteins": protein["unique_proteins"],
        "go_proteins": go["go_proteins"] if nonempty(go_table) else "NA",
    })
    result["scope"] = path_scope(directory, hap)
    return result


def locate_functional_dir(root, species, hap, functional_root):
    candidates = [
        functional_root / species / hap,
        root / "annotation" / "repartition" / "functional" / species / hap,
        root / "annotation" / "functional" / species / hap,
        root / "annotation" / "functional" / species,
        root / "data" / species / "haploid" / "annotation" / "functional" / hap,
    ]
    for directory in candidates:
        if directory.is_dir() and first_recursive([directory], [
            "protein_annotation.tsv", "swissprot.diamond.tsv", "go_annotation.tsv",
            "*.protein_annotation.tsv", "*.swissprot.diamond.tsv", "*.go_annotation.tsv",
        ]):
            return directory
    return None


def parse_postprocess_summary(path):
    values = parse_key_value_tsv(path)
    return {
        "low_confidence_models": values.get("low_confidence_models", "NA"),
        "pass_models": values.get("flag_PASS", "NA"),
    }


def locate_files(root, input_root, gene_root, lai_root, species, hap):
    prefix = f"{species}_{hap}"
    old_hap_root = root / "data" / species / "haploid"
    genome = first_existing([
        input_root / species / "genome" / f"{prefix}.fa",
        input_root / species / "genome" / f"{prefix}.fasta",
        old_hap_root / "genome" / f"{prefix}.fa",
        root / "final_haplotype_genomes_eviann_old" / f"{prefix}.fa",
        root / "final_haplotype_genomes_old" / f"{prefix}.fa",
        root / "data" / species / "haploid_rename" / "genome" / f"{prefix}.fa",
    ])

    edta_dirs = [
        input_root / species / "edta" / hap,
        old_hap_root / "annotation" / "repeats" / f"{hap}_edta",
        root / "data" / species / "haploid_rename" / "edta" / hap,
    ]
    te_gff = first_recursive(edta_dirs, [f"{prefix}.fa.mod.EDTA.TEanno.gff3", "*.EDTA.TEanno.gff3"])
    te_sum = first_recursive(edta_dirs, [f"{prefix}.fa.mod.EDTA.TEanno.sum", "*.EDTA.TEanno.sum"])

    eviann_dir = gene_root / species / f"{hap}_eviann" / "postprocess"
    helixer_dir = gene_root / species / f"{hap}_helixer"
    old_gene_dir = old_hap_root / "annotation" / "genes" / f"{hap}_helixer_edta_masked"
    annotation = first_existing([
        eviann_dir / f"{prefix}.eviann.genes.gff3",
        input_root / species / "gff" / f"{prefix}.gff3",
        old_gene_dir / "postprocess" / f"{prefix}.edta_masked.helixer.genes.gff3",
        root / "final_haplotype_genomes_eviann_old" / f"{prefix}.gff3",
        root / "final_haplotype_genomes_old" / f"{prefix}.gff3",
        root / "data" / species / "haploid_rename" / "gff" / f"{prefix}.gff3",
    ])
    helixer = first_existing([
        helixer_dir / f"{prefix}.edta_masked.helixer.gff3",
        old_gene_dir / f"{prefix}.edta_masked.helixer.gff3",
    ])

    cds = pep = post_summary = None
    if annotation:
        name = annotation.name
        if name.endswith(".genes.gff3"):
            stem = name[:-len(".genes.gff3")]
            cds_candidates = [annotation.parent / f"{stem}.cds.fa", annotation.parent / f"{stem}.cds.fasta"]
            pep_candidates = [annotation.parent / f"{stem}.protein.fa", annotation.parent / f"{stem}.pep.fa"]
            post_summary = annotation.parent / f"{stem}.postprocess_summary.tsv"
        else:
            stem = name.rsplit(".gff3", 1)[0]
            cds_candidates = [annotation.parent / f"{stem}.cds.fa"]
            pep_candidates = [annotation.parent / f"{stem}.pep.fa", annotation.parent / f"{stem}.protein.fa"]
        cds = first_existing(cds_candidates + [input_root / species / "cds" / f"{prefix}.cds.fa",
                                               root / "final_haplotype_genomes_eviann_old" / f"{prefix}.cds.fa",
                                               root / "final_haplotype_genomes_old" / f"{prefix}.cds.fa",
                                               root / "data" / species / "haploid_rename" / "cds" / f"{prefix}.cds.fa"])
        pep = first_existing(pep_candidates + [input_root / species / "pep" / f"{prefix}.pep.fa",
                                               root / "final_haplotype_genomes_eviann_old" / f"{prefix}.pep.fa",
                                               root / "final_haplotype_genomes_old" / f"{prefix}.pep.fa",
                                               root / "data" / species / "haploid_rename" / "pep" / f"{prefix}.pep.fa"])

    lai_dirs = [
        lai_root / species / f"{hap}_lai",
        root / "annotation" / "repartition" / "lai" / species / f"{hap}_lai",
        old_hap_root / "annotation" / "repeats" / f"{hap}_lai",
        root / "annotation" / "repeats" / "final_haplotypes" / f"{hap}_lai",
    ]
    lai = first_recursive(lai_dirs, ["*.LAI"])

    busco_dirs = [
        root / "haploid_repartition" / "busco",
        root / "annotation" / "repartition" / "busco",
        old_gene_dir / "busco",
        old_hap_root / "annotation" / "genes" / f"{hap}_helixer_edta_masked" / "busco",
    ]
    busco = first_recursive(busco_dirs, [f"*{prefix}*short_summary*.txt", "short_summary*.txt"])
    return {
        "genome": genome, "te_gff": te_gff, "te_sum": te_sum,
        "annotation": annotation, "helixer": helixer, "cds": cds,
        "pep": pep, "post_summary": post_summary, "lai": lai,
        "busco": busco,
    }


def annotation_source(path):
    if not path:
        return "missing"
    lower = str(path).lower()
    if "eviann" in lower:
        return "eviann_or_split_eviann"
    if "helixer" in lower:
        return "helixer"
    return "split_annotation"


def collect_row(root, input_root, gene_root, lai_root, ncrna_root, functional_root, species, hap):
    files = locate_files(root, input_root, gene_root, lai_root, species, hap)
    assembly = fasta_stats(files["genome"])
    annotation = gff_stats(files["annotation"])
    helixer = gff_stats(files["helixer"])
    edta = parse_edta_sum(files["te_sum"])
    ncrna = parse_ncrna(locate_ncrna_dir(root, species, hap, ncrna_root), hap)
    functional = parse_functional(locate_functional_dir(root, species, hap, functional_root), hap)
    busco = parse_busco(files["busco"])
    postprocess = parse_postprocess_summary(files["post_summary"])
    cds_records, cds_bp = fasta_record_count(files["cds"])
    pep_records, pep_aa = fasta_record_count(files["pep"])

    core_ok = nonempty(files["genome"]) and nonempty(files["annotation"])
    row = {
        "species": species,
        "haplotype": hap,
        "core_status": "pass" if core_ok else "missing_genome_or_annotation",
        "genome_fasta": str(files["genome"] or ""),
        "annotation_gff3": str(files["annotation"] or ""),
        "annotation_source": annotation_source(files["annotation"]),
        "annotation_status": "pass" if nonempty(files["annotation"]) else "missing",
        "helixer_gff3": str(files["helixer"] or ""),
        "helixer_status": "pass" if nonempty(files["helixer"]) else "missing",
        "edta_te_gff3": str(files["te_gff"] or ""),
        "edta_summary": str(files["te_sum"] or ""),
        "edta_status": "pass" if nonempty(files["te_sum"]) else "missing",
        "lai_result": str(files["lai"] or ""),
        "lai_whole_genome": parse_lai(files["lai"]),
        "lai_status": "pass" if nonempty(files["lai"]) else "missing",
        "ncrna_source": ncrna["source"],
        "ncrna_scope": ncrna["scope"],
        "ncrna_total_features": ncrna["total_features"],
        "trna_count": ncrna["trna_count"],
        "rfam_count": ncrna["rfam_count"],
        "ncrna_status": "pass" if ncrna["scope"] != "missing" else "missing",
        "functional_source": functional["source"],
        "functional_scope": functional["scope"],
        "functional_hit_rows": functional["hit_rows"],
        "functional_annotated_proteins": functional["annotated_proteins"],
        "go_annotated_proteins": functional["go_proteins"],
        "functional_status": "pass" if functional["scope"] != "missing" else "missing",
        "busco_summary": str(files["busco"] or ""),
        "busco_status": "pass" if nonempty(files["busco"]) else "missing",
        "postprocess_summary": str(files["post_summary"] or ""),
        "cds_fasta": str(files["cds"] or ""),
        "protein_fasta": str(files["pep"] or ""),
        "cds_records": cds_records,
        "cds_fasta_bp": cds_bp,
        "protein_records": pep_records,
        "protein_fasta_aa": pep_aa,
    }
    row.update(assembly)
    row.update(annotation)
    row.update({
        "helixer_gene_count": helixer["gene_count"],
        "helixer_mrna_count": helixer["mrna_count"],
        "helixer_cds_count": helixer["cds_count"],
        "helixer_exon_count": helixer["exon_count"],
    })
    row.update(edta)
    row.update(busco)
    row.update(postprocess)
    return row


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--input-root", type=Path)
    parser.add_argument("--gene-root", type=Path)
    parser.add_argument("--lai-root", type=Path)
    parser.add_argument("--ncrna-root", type=Path)
    parser.add_argument("--functional-root", type=Path)
    parser.add_argument("--species", nargs="+", default=list(HAPLOTYPES))
    parser.add_argument("--strict", action="store_true",
                        help="fail if a genome or primary annotation is missing")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.root
    input_root = args.input_root or root / "haploid_repartition" / "inputs"
    gene_root = args.gene_root or root / "annotation" / "repartition" / "genes"
    lai_root = args.lai_root or root / "annotation" / "repartition" / "lai"
    ncrna_root = args.ncrna_root or root / "annotation" / "repartition" / "ncrna"
    functional_root = args.functional_root or root / "annotation" / "repartition" / "functional"
    unknown = [species for species in args.species if species not in HAPLOTYPES]
    if unknown:
        print("ERROR: unsupported species: {}".format(", ".join(unknown)), file=sys.stderr)
        return 2

    rows = []
    for species in args.species:
        for hap in HAPLOTYPES[species]:
            rows.append(collect_row(root, input_root, gene_root, lai_root,
                                    ncrna_root, functional_root, species, hap))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    columns = [
        "species", "haplotype", "core_status",
        "genome_fasta", "genome_bp", "scaffold_count", "contig_count",
        "scaffold_n50", "scaffold_n90", "contig_n50", "contig_n90",
        "gc_pct", "gap_bp",
        "annotation_gff3", "annotation_source", "annotation_status",
        "gene_count", "mrna_count", "cds_count", "exon_count", "cds_bp", "exon_bp",
        "cds_fasta", "cds_records", "cds_fasta_bp", "protein_fasta", "protein_records", "protein_fasta_aa",
        "helixer_gff3", "helixer_status", "helixer_gene_count", "helixer_mrna_count",
        "helixer_cds_count", "helixer_exon_count",
        "edta_te_gff3", "edta_summary", "edta_status", "edta_genome_bp",
        "edta_repeat_bp", "edta_repeat_pct", "edta_summary_lines",
        "lai_result", "lai_whole_genome", "lai_status",
        "ncrna_source", "ncrna_scope", "ncrna_total_features", "trna_count", "rfam_count", "ncrna_status",
        "functional_source", "functional_scope", "functional_hit_rows",
        "functional_annotated_proteins", "go_annotated_proteins", "functional_status",
        "busco_summary", "busco_status", "busco_complete_pct", "busco_single_copy_pct",
        "busco_duplicated_pct", "busco_fragmented_pct", "busco_missing_pct",
        "postprocess_summary", "low_confidence_models", "pass_models",
        "gff_feature_count", "gff_bad_rows",
    ]
    output = args.output_dir / "haplotype_annotation_summary.tsv"
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    manifest = args.output_dir / "haplotype_annotation_summary_manifest.txt"
    missing_core = [f"{row['species']}/{row['haplotype']}" for row in rows if row["core_status"] != "pass"]
    species_scoped_ncrna = sorted({row["species"] for row in rows if row["ncrna_scope"] == "species_complete"})
    species_scoped_functional = sorted({row["species"] for row in rows if row["functional_scope"] == "species_complete"})
    with manifest.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("Read-only haplotype annotation summary\n")
        handle.write("Input root: {}\n".format(input_root))
        handle.write("Gene root: {}\n".format(gene_root))
        handle.write("LAI root: {}\n".format(lai_root))
        handle.write("Rows: {}\n".format(len(rows)))
        handle.write("Missing genome/primary annotation: {}\n".format(", ".join(missing_core) or "none"))
        handle.write("ncRNA rows backed only by complete-species files: {}\n".format(", ".join(species_scoped_ncrna) or "none"))
        handle.write("Functional rows backed only by complete-species files: {}\n".format(", ".join(species_scoped_functional) or "none"))
        handle.write("No annotation, ncRNA, functional or LAI program was rerun.\n")

    print("[OK] Summary: {}".format(output))
    print("[OK] Manifest: {}".format(manifest))
    if missing_core:
        print("[WARN] Missing core rows: {}".format(", ".join(missing_core)), file=sys.stderr)
        if args.strict:
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
