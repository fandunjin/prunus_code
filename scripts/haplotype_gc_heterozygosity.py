#!/usr/bin/env python3
"""Prepare phased references and summarize GC and residual heterozygosity."""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path


HAPS = ("A1", "A2", "B1", "B2")
GROUPS = {
    "A1": ("A1",),
    "A2": ("A2",),
    "B1": ("B1",),
    "B2": ("B2",),
    "copy1": ("A1", "B1"),
    "copy2": ("A2", "B2"),
}


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="ascii", errors="replace")
    return path.open("r", encoding="ascii", errors="replace")


def parse_assignments(items: list[str]) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for item in items:
        if "=" not in item:
            raise SystemExit(f"Expected LABEL=FASTA, received: {item}")
        label, value = item.split("=", 1)
        if label not in HAPS:
            raise SystemExit(f"Unsupported haplotype label: {label}")
        path = Path(value)
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Missing or empty FASTA: {path}")
        result[label] = path
    missing = [label for label in HAPS if label not in result]
    if missing:
        raise SystemExit(f"Missing FASTA assignments: {', '.join(missing)}")
    return result


def fasta_records(path: Path):
    name = None
    chunks: list[str] = []
    with open_text(path) as handle:
        for raw in handle:
            if raw.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = raw[1:].strip().split()[0]
                if not name:
                    raise SystemExit(f"Empty FASTA identifier in {path}")
                chunks = []
            else:
                chunks.append(raw.strip())
    if name is not None:
        yield name, "".join(chunks)


def write_reference(fastas: dict[str, Path], output: Path, mapping: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    seen: set[str] = set()
    with output.open("w", encoding="ascii") as fasta_out, mapping.open(
        "w", encoding="utf-8", newline=""
    ) as map_out:
        writer = csv.writer(map_out, delimiter="\t", lineterminator="\n")
        writer.writerow(("seqid", "haplotype", "length"))
        for label in HAPS:
            for seqid, sequence in fasta_records(fastas[label]):
                if seqid in seen:
                    raise SystemExit(f"Duplicate sequence ID across haplotypes: {seqid}")
                seen.add(seqid)
                writer.writerow((seqid, label, len(sequence)))
                fasta_out.write(f">{seqid}\n")
                for start in range(0, len(sequence), 80):
                    fasta_out.write(sequence[start : start + 80] + "\n")


def sequence_stats(fastas: dict[str, Path]) -> dict[str, dict[str, int]]:
    stats: dict[str, dict[str, int]] = {}
    for label in HAPS:
        counts = defaultdict(int)
        for _, sequence in fasta_records(fastas[label]):
            counts["sequences"] += 1
            upper = sequence.upper()
            counts["total_bp"] += len(upper)
            counts["a_bp"] += upper.count("A")
            counts["c_bp"] += upper.count("C")
            counts["g_bp"] += upper.count("G")
            counts["t_bp"] += upper.count("T")
            counts["n_bp"] += upper.count("N")
        counts["acgt_bp"] = sum(counts[f"{base}_bp"] for base in "acgt")
        counts["gc_bp"] = counts["g_bp"] + counts["c_bp"]
        counts["other_bp"] = counts["total_bp"] - counts["acgt_bp"] - counts["n_bp"]
        stats[label] = dict(counts)
    return stats


def write_gc(fastas: dict[str, Path], output: Path) -> None:
    hap_stats = sequence_stats(fastas)
    fields = (
        "label",
        "source_haplotypes",
        "sequences",
        "total_bp",
        "acgt_bp",
        "gc_bp",
        "n_bp",
        "other_bp",
        "gc_percent_acgt",
    )
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for label, members in GROUPS.items():
            row = {key: sum(hap_stats[x].get(key, 0) for x in members) for key in (
                "sequences", "total_bp", "acgt_bp", "gc_bp", "n_bp", "other_bp"
            )}
            row["label"] = label
            row["source_haplotypes"] = "+".join(members)
            row["gc_percent_acgt"] = (
                f"{100.0 * row['gc_bp'] / row['acgt_bp']:.6f}" if row["acgt_bp"] else "NA"
            )
            writer.writerow(row)


def read_map(path: Path) -> tuple[dict[str, str], dict[str, int]]:
    labels: dict[str, str] = {}
    lengths: dict[str, int] = {}
    with path.open("r", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            labels[row["seqid"]] = row["haplotype"]
            lengths[row["seqid"]] = int(row["length"])
    return labels, lengths


def read_gc(path: Path) -> dict[str, dict[str, str]]:
    with path.open("r", encoding="utf-8") as handle:
        return {row["label"]: row for row in csv.DictReader(handle, delimiter="\t")}


def variant_counts(
    path: Path,
    contig_labels: dict[str, str],
    min_qual: float,
    min_depth: int,
    max_depth: int,
    min_alt_fraction: float,
    max_alt_fraction: float,
) -> tuple[dict[str, int], dict[str, int]]:
    passing = defaultdict(int)
    all_het = defaultdict(int)
    with open_text(path) as handle:
        for raw in handle:
            if not raw or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom, _, _, ref, alt, qual_text = fields[:6]
            label = contig_labels.get(chrom)
            if label is None or len(ref) != 1 or len(alt) != 1 or "," in alt:
                continue
            format_keys = fields[8].split(":")
            sample_values = fields[9].split(":")
            sample = dict(zip(format_keys, sample_values))
            gt = sample.get("GT", "").replace("|", "/")
            if gt not in {"0/1", "1/0"}:
                continue
            all_het[label] += 1
            try:
                qual = float(qual_text)
                depth = int(sample.get("DP", "0"))
                ad = [int(x) for x in sample.get("AD", "").split(",")[:2]]
            except (ValueError, TypeError):
                continue
            if len(ad) != 2 or sum(ad) == 0:
                continue
            alt_fraction = ad[1] / sum(ad)
            if (
                qual >= min_qual
                and min_depth <= depth <= max_depth
                and min_alt_fraction <= alt_fraction <= max_alt_fraction
            ):
                passing[label] += 1
    return dict(passing), dict(all_het)


def coverage_stats(
    path: Path, contig_labels: dict[str, str], contig_lengths: dict[str, int]
) -> dict[str, dict[str, float]]:
    totals: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    with path.open("r", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            seqid = row.get("#rname") or row.get("rname")
            label = contig_labels.get(seqid or "")
            if label is None:
                continue
            span = contig_lengths[seqid]
            totals[label]["span"] += span
            totals[label]["covered_bp"] += int(row["covbases"])
            totals[label]["depth_sum"] += float(row["meandepth"]) * span
    return {label: dict(values) for label, values in totals.items()}


def callable_stats(path: Path | None, contig_labels: dict[str, str]) -> dict[str, int]:
    if path is None:
        return {}
    totals = defaultdict(int)
    with open_text(path) as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            label = contig_labels.get(fields[0])
            if label is None:
                continue
            try:
                at_min_depth = int(fields[-2])
                above_max_depth = int(fields[-1])
            except ValueError:
                continue
            totals[label] += max(0, at_min_depth - above_max_depth)
    return dict(totals)


def write_summary(a: argparse.Namespace) -> None:
    gc = read_gc(a.gc)
    labels, lengths = read_map(a.contig_map)
    passing, all_het = variant_counts(
        a.vcf,
        labels,
        a.min_qual,
        a.min_depth,
        a.max_depth,
        a.min_alt_fraction,
        a.max_alt_fraction,
    )
    coverage = coverage_stats(a.coverage, labels, lengths)
    callable_by_hap = callable_stats(a.callable, labels)
    fields = (
        "label", "source_haplotypes", "total_bp", "acgt_bp", "gc_percent_acgt",
        "covered_bp", "covered_percent", "callable_bp", "denominator_type", "mean_depth", "raw_het_snps",
        "filtered_het_snps", "heterozygosity_rate_per_bp", "heterozygosity_percent",
        "het_snps_per_kb", "min_qual", "min_depth", "max_depth", "alt_fraction_range",
    )
    with a.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for label, members in GROUPS.items():
            acgt_bp = int(gc[label]["acgt_bp"])
            covered_bp = int(sum(coverage.get(x, {}).get("covered_bp", 0) for x in members))
            span = sum(coverage.get(x, {}).get("span", 0) for x in members)
            depth_sum = sum(coverage.get(x, {}).get("depth_sum", 0) for x in members)
            callable_bp = sum(callable_by_hap.get(x, 0) for x in members)
            filtered = sum(passing.get(x, 0) for x in members)
            raw = sum(all_het.get(x, 0) for x in members)
            denominator = callable_bp if a.callable is not None else covered_bp
            denominator_type = (
                f"callable_depth_{a.min_depth}_{a.max_depth}"
                if a.callable is not None
                else "covered_bp_fallback"
            )
            rate = filtered / denominator if denominator else 0.0
            writer.writerow({
                "label": label,
                "source_haplotypes": "+".join(members),
                "total_bp": gc[label]["total_bp"],
                "acgt_bp": acgt_bp,
                "gc_percent_acgt": gc[label]["gc_percent_acgt"],
                "covered_bp": covered_bp,
                "covered_percent": f"{100.0 * covered_bp / span:.6f}" if span else "NA",
                "callable_bp": callable_bp if a.callable is not None else "NA",
                "denominator_type": denominator_type,
                "mean_depth": f"{depth_sum / span:.6f}" if span else "NA",
                "raw_het_snps": raw,
                "filtered_het_snps": filtered,
                "heterozygosity_rate_per_bp": f"{rate:.10g}",
                "heterozygosity_percent": f"{100.0 * rate:.8f}",
                "het_snps_per_kb": f"{1000.0 * rate:.8f}",
                "min_qual": a.min_qual,
                "min_depth": a.min_depth,
                "max_depth": a.max_depth,
                "alt_fraction_range": f"{a.min_alt_fraction}-{a.max_alt_fraction}",
            })


def last_merqury_row(path: Path) -> list[str]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    rows: list[list[str]] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if fields and fields[0].lower() not in {"assembly", "asm"}:
                rows.append(fields)
    return rows[-1] if rows else []


def write_merqury_summary(a: argparse.Namespace) -> None:
    gc = read_gc(a.gc)
    runs: dict[str, Path] = {}
    for item in a.run:
        if "=" not in item:
            raise SystemExit(f"Expected LABEL=OUTPUT_PREFIX, received: {item}")
        label, prefix = item.split("=", 1)
        runs[label] = Path(prefix)
    missing = [label for label in GROUPS if label not in runs]
    if missing:
        raise SystemExit(f"Missing Merqury output prefixes: {', '.join(missing)}")

    fields = (
        "label", "source_haplotypes", "sequences", "total_bp", "acgt_bp",
        "gc_percent_acgt", "merqury_qv", "merqury_error_rate",
        "merqury_completeness_percent", "qv_file", "completeness_file",
        "spectra_cn_hist", "interpretation",
    )
    with a.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for label, members in GROUPS.items():
            prefix = runs[label]
            qv_file = Path(str(prefix) + ".qv")
            completeness_file = Path(str(prefix) + ".completeness.stats")
            spectra_hist = Path(str(prefix) + ".spectra-cn.hist")
            qv_row = last_merqury_row(qv_file)
            completeness_row = last_merqury_row(completeness_file)
            writer.writerow({
                "label": label,
                "source_haplotypes": "+".join(members),
                "sequences": gc[label]["sequences"],
                "total_bp": gc[label]["total_bp"],
                "acgt_bp": gc[label]["acgt_bp"],
                "gc_percent_acgt": gc[label]["gc_percent_acgt"],
                "merqury_qv": qv_row[-2] if len(qv_row) >= 2 else "NA",
                "merqury_error_rate": qv_row[-1] if qv_row else "NA",
                "merqury_completeness_percent": completeness_row[-1] if completeness_row else "NA",
                "qv_file": str(qv_file),
                "completeness_file": str(completeness_file),
                "spectra_cn_hist": str(spectra_hist),
                "interpretation": "Merqury error rate is not biological heterozygosity; inspect spectra-cn/asm for residual duplication and haplotype separation.",
            })


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    reference = sub.add_parser("reference")
    reference.add_argument("--fasta", action="append", required=True)
    reference.add_argument("--output", required=True, type=Path)
    reference.add_argument("--contig-map", required=True, type=Path)

    gc = sub.add_parser("gc")
    gc.add_argument("--fasta", action="append", required=True)
    gc.add_argument("--output", required=True, type=Path)

    summary = sub.add_parser("summarize")
    summary.add_argument("--gc", required=True, type=Path)
    summary.add_argument("--contig-map", required=True, type=Path)
    summary.add_argument("--vcf", required=True, type=Path)
    summary.add_argument("--coverage", required=True, type=Path)
    summary.add_argument("--callable", type=Path)
    summary.add_argument("--output", required=True, type=Path)
    summary.add_argument("--min-qual", type=float, default=30.0)
    summary.add_argument("--min-depth", type=int, default=8)
    summary.add_argument("--max-depth", type=int, default=200)
    summary.add_argument("--min-alt-fraction", type=float, default=0.25)
    summary.add_argument("--max-alt-fraction", type=float, default=0.75)

    merqury = sub.add_parser("merqury-summary")
    merqury.add_argument("--gc", required=True, type=Path)
    merqury.add_argument("--run", action="append", required=True)
    merqury.add_argument("--output", required=True, type=Path)

    a = parser.parse_args()
    if a.command == "reference":
        write_reference(parse_assignments(a.fasta), a.output, a.contig_map)
    elif a.command == "gc":
        write_gc(parse_assignments(a.fasta), a.output)
    elif a.command == "summarize":
        write_summary(a)
    else:
        write_merqury_summary(a)


if __name__ == "__main__":
    main()
