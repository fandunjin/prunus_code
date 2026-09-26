#!/usr/bin/env python3
"""Summarize PSP A1/A2/B1/B2 assembly, repeat, LAI and gene statistics."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


HAPLOTYPES = ("A1", "A2", "B1", "B2")


def nx(lengths: list[int], fraction: float) -> int:
    if not lengths:
        return 0
    target = sum(lengths) * fraction
    cumulative = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= target:
            return length
    return 0


def fasta_stats(path: Path) -> dict[str, float | int]:
    scaffold_lengths: list[int] = []
    contig_lengths: list[int] = []
    current_length = 0
    current_contig = 0
    gc = 0
    acgt = 0
    gaps = 0

    with path.open(encoding="ascii") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_length:
                    scaffold_lengths.append(current_length)
                if current_contig:
                    contig_lengths.append(current_contig)
                current_length = 0
                current_contig = 0
                continue

            seq = line.upper()
            current_length += len(seq)
            gc += seq.count("G") + seq.count("C")
            acgt += sum(seq.count(base) for base in "ACGT")
            gaps += seq.count("N")
            for segment in re.split(r"(N+)", seq):
                if not segment:
                    continue
                if segment.startswith("N"):
                    if current_contig:
                        contig_lengths.append(current_contig)
                        current_contig = 0
                else:
                    current_contig += len(segment)

    if current_length:
        scaffold_lengths.append(current_length)
    if current_contig:
        contig_lengths.append(current_contig)

    return {
        "Genome size (bp)": sum(scaffold_lengths),
        "Number of contigs": len(contig_lengths),
        "Contig N50 (bp)": nx(contig_lengths, 0.50),
        "Contig N75 (bp)": nx(contig_lengths, 0.75),
        "Number of scaffolds": len(scaffold_lengths),
        "Scaffold N50 (bp)": nx(scaffold_lengths, 0.50),
        "GC (%)": round(100 * gc / acgt, 2) if acgt else 0.0,
        "Gap bases": gaps,
    }


def parse_lai(path: Path) -> float:
    with path.open(encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields and fields[0] == "whole_genome":
                return float(fields[6])
    raise ValueError(f"No whole_genome LAI row in {path}")


def parse_edta(path: Path) -> dict[str, float | int]:
    text = path.read_text(encoding="utf-8")
    top = text.split("Repeat Stats", 1)[0]
    total_length_match = re.search(r"^Total Length:\s+(\d+)\s+bp$", top, re.MULTILINE)
    if not total_length_match:
        raise ValueError(f"No Total Length in {path}")

    category = None
    subclasses: dict[str, list[int]] = {}
    total_bp = None
    total_percent = None
    for line in top.splitlines():
        category_match = re.match(r"^(LINE|SINE|LTR|TIR|nonLTR|nonTIR|rDNA)\s+--", line)
        if category_match:
            category = category_match.group(1)
            subclasses.setdefault(category, [])
            continue

        subclass_match = re.match(r"^\s{2,}(\S+)\s+\d+\s+(\d+)\s+([\d.]+)%", line)
        if subclass_match and category:
            subclasses.setdefault(category, []).append(int(subclass_match.group(2)))
            continue

        total_match = re.match(r"^Total\s+\d+\s+(\d+)\s+([\d.]+)%", line)
        if total_match:
            total_bp = int(total_match.group(1))
            total_percent = float(total_match.group(2))

    def subclass_bp(parent: str, name: str) -> int:
        pattern = rf"^\s+{re.escape(name)}\s+\d+\s+(\d+)\s+[\d.]+%"
        section = top.split(f"{parent}", 1)[1] if parent in top else ""
        match = re.search(pattern, section, re.MULTILINE | re.IGNORECASE)
        return int(match.group(1)) if match else 0

    if total_bp is None or total_percent is None:
        raise ValueError(f"No final Total row in {path}")

    return {
        "EDTA genome size": int(total_length_match.group(1)),
        "LINE (bp)": sum(subclasses.get("LINE", [])),
        "SINE (bp)": sum(subclasses.get("SINE", [])),
        "LTR (bp)": sum(subclasses.get("LTR", [])),
        "LTR-Gypsy (bp)": subclass_bp("LTR", "Gypsy"),
        "LTR-Copia (bp)": subclass_bp("LTR", "Copia"),
        "LTR-Unknown (bp)": subclass_bp("LTR", "unknown"),
        "Unclassified (bp)": 0,
        "Total repeat (bp)": total_bp,
        "Total repeat (%)": total_percent,
    }


def parse_postprocess(path: Path) -> dict[str, int]:
    metrics: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            key, value = line.rstrip("\n").split("\t", 1)
            metrics[key] = value
    return {
        "Protein-coding gene": int(metrics["gff3_feature_gene"]),
        "mRNA": int(metrics["gff3_feature_mRNA"]),
        "CDS": int(metrics["gff3_feature_CDS"]),
        "Exon": int(metrics["gff3_feature_exon"]),
    }


def collect(root: Path, hap: str) -> dict[str, float | int]:
    tag = f"hap{hap}"
    stats = fasta_stats(root / "genome" / f"psp_{tag}.fa")
    stats["LAI"] = parse_lai(
        root / "annotation" / "repeats" / f"{tag}_lai" / f"psp_{tag}.all_ltr.out.LAI"
    )
    edta = parse_edta(
        root
        / "annotation"
        / "repeats"
        / f"{tag}_edta"
        / f"psp_{tag}.fa.mod.EDTA.TEanno.sum"
    )
    edta_genome_size = edta.pop("EDTA genome size")
    if stats["Genome size (bp)"] - stats["Gap bases"] != edta_genome_size:
        raise ValueError(
            f"Non-N genome size differs between FASTA and EDTA for {hap}: "
            f"{stats['Genome size (bp)'] - stats['Gap bases']} versus {edta_genome_size}"
        )
    stats.pop("Gap bases")
    stats.update(edta)
    stats.update(
        parse_postprocess(
            root
            / "annotation"
            / "genes"
            / f"{tag}_helixer_edta_masked"
            / "postprocess"
            / f"psp_{tag}.edta_masked.helixer.postprocess_summary.tsv"
        )
    )
    return stats


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    data = {hap: collect(args.root, hap) for hap in HAPLOTYPES}
    rows = [
        "Genome size (bp)",
        "Number of contigs",
        "Contig N50 (bp)",
        "Contig N75 (bp)",
        "Number of scaffolds",
        "Scaffold N50 (bp)",
        "GC (%)",
        "LAI",
        "Repeat elements",
        "LINE (bp)",
        "SINE (bp)",
        "LTR (bp)",
        "LTR-Gypsy (bp)",
        "LTR-Copia (bp)",
        "LTR-Unknown (bp)",
        "Unclassified (bp)",
        "Total (bp)",
        "Total (%)",
        "Protein-coding gene",
        "Gene number",
        "mRNA",
        "CDS",
        "Exon",
    ]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("Haplotype\t" + "\t".join(HAPLOTYPES) + "\n")
        for row in rows:
            if row in {"Repeat elements", "Protein-coding gene"}:
                handle.write(row + "\t\t\t\t\n")
                continue
            source_row = {
                "Gene number": "Protein-coding gene",
                "Total (bp)": "Total repeat (bp)",
                "Total (%)": "Total repeat (%)",
            }.get(row, row)
            values = [data[hap][source_row] for hap in HAPLOTYPES]
            handle.write(row + "\t" + "\t".join(str(value) for value in values) + "\n")


if __name__ == "__main__":
    main()
