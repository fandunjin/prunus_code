#!/usr/bin/env python3
"""Build AGR input tables from OrthoFinder and normalized JCVI inputs."""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import tempfile
from collections import Counter
from pathlib import Path


COMMON_SUFFIXES = (
    ".pep.fa",
    ".pep.faa",
    ".protein.fa",
    ".protein.faa",
    ".fasta",
    ".faa",
    ".fa",
    ".pep",
    ".protein",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create orthogroups.csv, agr_gff.csv and agr_chr.csv for AGR."
    )
    parser.add_argument("--orthogroups", required=True, type=Path)
    parser.add_argument("--jcvi-input-root", required=True, type=Path)
    parser.add_argument("--genome-root", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--species", nargs="+", default=None)
    parser.add_argument("--min-match-rate", type=float, default=0.90)
    return parser.parse_args()


def die(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def strip_known_suffixes(name: str) -> list[str]:
    candidates = [name]
    current = name
    changed = True
    while changed:
        changed = False
        lower = current.lower()
        for suffix in COMMON_SUFFIXES:
            if lower.endswith(suffix):
                current = current[: -len(suffix)]
                if current and current not in candidates:
                    candidates.append(current)
                changed = True
                break
    return candidates


def resolve_species_column(raw: str, available: set[str]) -> str | None:
    raw = raw.strip()
    candidates = strip_known_suffixes(raw)
    candidates.extend(
        candidate.replace("-", "_") for candidate in list(candidates)
    )
    matches = []
    for candidate in candidates:
        if candidate in available and candidate not in matches:
            matches.append(candidate)
    if len(matches) > 1:
        die(f"Ambiguous OrthoFinder species column {raw!r}: {matches}")
    return matches[0] if matches else None


def read_seqids(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8", errors="replace")
    seqids = [item.strip() for item in re.split(r"[\s,]+", text) if item.strip()]
    if not seqids:
        die(f"No sequence IDs found in {path}")
    if len(seqids) != len(set(seqids)):
        die(f"Duplicate sequence IDs found in {path}")
    return seqids


def read_bed(path: Path, allowed_seqids: set[str]) -> list[tuple[str, int, int, str, str]]:
    records = []
    seen_genes = set()
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                die(f"BED has fewer than 6 columns at {path}:{line_number}")
            seqid, start0, end, gene, _, strand = fields[:6]
            if seqid not in allowed_seqids:
                continue
            if gene in seen_genes:
                die(f"Duplicate BED GeneName {gene!r} in {path}")
            seen_genes.add(gene)
            try:
                start1 = int(start0) + 1
                stop = int(end)
            except ValueError:
                die(f"Invalid BED coordinate at {path}:{line_number}")
            if start1 < 1 or stop < start1:
                die(f"Invalid BED interval at {path}:{line_number}")
            if strand not in {"+", "-"}:
                strand = "."
            records.append((seqid, start1, stop, gene, strand))
    if not records:
        die(f"No usable BED records found in {path}")
    return records


def read_fasta_lengths(path: Path, wanted: set[str]) -> dict[str, int]:
    lengths: dict[str, int] = {}
    current = None
    current_length = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if current in wanted:
                    lengths[current] = current_length
                current = line[1:].strip().split()[0]
                current_length = 0
            elif current in wanted:
                current_length += len(line.strip())
    if current in wanted:
        lengths[current] = current_length
    missing = sorted(wanted - lengths.keys())
    if missing:
        preview = ", ".join(missing[:10])
        raise ValueError(
            f"{path} lacks {len(missing)} required sequence IDs: {preview}"
        )
    return lengths


def infer_lengths_from_bed(
    seqids: list[str], records: list[tuple[str, int, int, str, str]]
) -> dict[str, int]:
    lengths = {seqid: 0 for seqid in seqids}
    for seqid, _, stop, _, _ in records:
        lengths[seqid] = max(lengths[seqid], stop)
    missing = [seqid for seqid, length in lengths.items() if length <= 0]
    if missing:
        die(
            "Cannot infer chromosome lengths from BED for: "
            + ", ".join(missing[:10])
        )
    return lengths


def split_orthofinder_cell(cell: str) -> list[str]:
    if not cell.strip():
        return []
    return [item.strip() for item in cell.split(",") if item.strip()]


def main() -> None:
    args = parse_args()
    for path in (args.orthogroups, args.jcvi_input_root, args.genome_root):
        if not path.exists():
            die(f"Required path does not exist: {path}")
    if not 0 <= args.min_match_rate <= 1:
        die("--min-match-rate must be between 0 and 1")

    available = {
        path.name
        for path in args.jcvi_input_root.iterdir()
        if path.is_dir()
    }
    if args.species:
        selected_species = list(dict.fromkeys(args.species))
        missing = [species for species in selected_species if species not in available]
        if missing:
            die(f"Species missing from JCVI input root: {', '.join(missing)}")
    else:
        selected_species = sorted(available)

    species_data = {}
    for species in selected_species:
        input_dir = args.jcvi_input_root / species
        bed_path = input_dir / f"{species}.bed"
        seqids_path = input_dir / f"{species}.seqids"
        for path in (bed_path, seqids_path):
            if not path.is_file():
                die(f"Missing required input for {species}: {path}")
        seqids = read_seqids(seqids_path)
        records = read_bed(bed_path, set(seqids))
        genome_candidates = [
            (args.outdir / "genome" / f"{species}.genome.fa", "agr_genome_copy"),
            (args.genome_root / species / "genome.fa", "staged_genome_fasta"),
        ]
        lengths = None
        length_source = ""
        genome_errors = []
        for genome_path, source_name in genome_candidates:
            if not genome_path.is_file():
                continue
            try:
                lengths = read_fasta_lengths(genome_path, set(seqids))
                length_source = source_name
                break
            except ValueError as error:
                genome_errors.append(str(error))
        if lengths is None:
            lengths = infer_lengths_from_bed(seqids, records)
            if genome_errors:
                length_source = "BED_max_gene_end_invalid_genome_fasta"
                print(f"WARNING: no valid chromosome FASTA was found for {species}:")
                for error in genome_errors:
                    print(f"  {error}")
            else:
                length_source = "BED_max_gene_end_no_genome_fasta"
                print(f"WARNING: no chromosome FASTA was found for {species}.")
            print("  ChrSize will use the last gene end from the normalized BED.")
        gene_index = {record[3] for record in records}
        gene_counts = Counter(record[0] for record in records)
        species_data[species] = {
            "seqids": seqids,
            "records": records,
            "lengths": lengths,
            "length_source": length_source,
            "gene_index": gene_index,
            "gene_counts": gene_counts,
        }

    args.outdir.mkdir(parents=True, exist_ok=True)
    tmpdir = Path(tempfile.mkdtemp(prefix="agr_prepare_", dir=args.outdir))
    try:
        gff_tmp = tmpdir / "agr_gff.csv"
        chr_tmp = tmpdir / "agr_chr.csv"
        og_tmp = tmpdir / "orthogroups.csv"
        summary_tmp = tmpdir / "species_summary.tsv"
        unmatched_tmp = tmpdir / "unmatched_orthofinder_genes.tsv"

        with gff_tmp.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter=";", lineterminator="\n")
            writer.writerow(["Species", "GeneName", "ChrName", "Start", "Stop", "Strand"])
            for species in selected_species:
                for seqid, start, stop, gene, strand in species_data[species]["records"]:
                    writer.writerow([species, gene, seqid, start, stop, strand])

        with chr_tmp.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter=";", lineterminator="\n")
            writer.writerow(["Species", "ChrName", "ChrSize", "ChrNum", "Effectif"])
            for species in selected_species:
                data = species_data[species]
                for chr_num, seqid in enumerate(data["seqids"], 1):
                    writer.writerow(
                        [
                            species,
                            seqid,
                            data["lengths"][seqid],
                            chr_num,
                            data["gene_counts"].get(seqid, 0),
                        ]
                    )

        total_by_species = Counter()
        matched_by_species = Counter()
        unmatched_rows = []
        with args.orthogroups.open("r", encoding="utf-8-sig", errors="replace", newline="") as source, og_tmp.open(
            "w", encoding="utf-8", newline=""
        ) as destination:
            reader = csv.reader(source, delimiter="\t")
            writer = csv.writer(destination, delimiter=";", lineterminator="\n")
            header = next(reader, None)
            if not header or header[0] != "Orthogroup":
                die(f"Unexpected OrthoFinder header in {args.orthogroups}")

            selected_columns = {}
            for column_index, raw_name in enumerate(header[1:], 1):
                resolved = resolve_species_column(raw_name, set(selected_species))
                if resolved:
                    if resolved in selected_columns.values():
                        die(f"Multiple OrthoFinder columns resolve to {resolved}")
                    selected_columns[column_index] = resolved
            unresolved = [species for species in selected_species if species not in selected_columns.values()]
            if unresolved:
                die(
                    "Selected species missing from OrthoFinder header after name normalization: "
                    + ", ".join(unresolved)
                )

            writer.writerow(["Orthogroup", "Species", "GeneName"])
            for row in reader:
                if not row:
                    continue
                orthogroup = row[0].strip()
                if not orthogroup:
                    continue
                for column_index, species in selected_columns.items():
                    cell = row[column_index] if column_index < len(row) else ""
                    for gene in split_orthofinder_cell(cell):
                        total_by_species[species] += 1
                        if gene in species_data[species]["gene_index"]:
                            writer.writerow([orthogroup, species, gene])
                            matched_by_species[species] += 1
                        else:
                            unmatched_rows.append((species, orthogroup, gene))

        failed = []
        with summary_tmp.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "Species",
                    "Chromosomes_or_sequences",
                    "ChrSize_source",
                    "BED_genes",
                    "OrthoFinder_genes",
                    "Matched_genes",
                    "Match_rate",
                ]
            )
            for species in selected_species:
                total = total_by_species[species]
                matched = matched_by_species[species]
                rate = matched / total if total else 0.0
                writer.writerow(
                    [
                        species,
                        len(species_data[species]["seqids"]),
                        species_data[species]["length_source"],
                        len(species_data[species]["records"]),
                        total,
                        matched,
                        f"{rate:.6f}",
                    ]
                )
                if rate < args.min_match_rate:
                    failed.append(f"{species}={rate:.3%}")

        with unmatched_tmp.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["Species", "Orthogroup", "GeneName"])
            writer.writerows(unmatched_rows)

        if failed:
            os.replace(summary_tmp, args.outdir / "species_summary.tsv")
            os.replace(
                unmatched_tmp,
                args.outdir / "unmatched_orthofinder_genes.tsv",
            )
            die(
                "OrthoFinder-to-BED match rate is below the required threshold for: "
                + ", ".join(failed)
                + f". Inspect {args.outdir / 'species_summary.tsv'} and "
                + f"{args.outdir / 'unmatched_orthofinder_genes.tsv'}."
            )

        readme_tmp = tmpdir / "README.txt"
        readme_tmp.write_text(
            "AGR input tables generated from:\n"
            f"  OrthoFinder: {args.orthogroups}\n"
            f"  JCVI inputs: {args.jcvi_input_root}\n"
            f"  Genomes: {args.genome_root}\n"
            f"  Species: {' '.join(selected_species)}\n\n"
            "Files:\n"
            "  orthogroups.csv: Orthogroup;Species;GeneName\n"
            "  agr_gff.csv: Species;GeneName;ChrName;Start;Stop;Strand\n"
            "  agr_chr.csv: Species;ChrName;ChrSize;ChrNum;Effectif\n"
            "  species_summary.tsv: identifier matching and chromosome-count QC\n"
            "  unmatched_orthofinder_genes.tsv: excluded OrthoFinder IDs\n",
            encoding="utf-8",
        )

        for name in (
            "orthogroups.csv",
            "agr_gff.csv",
            "agr_chr.csv",
            "species_summary.tsv",
            "unmatched_orthofinder_genes.tsv",
            "README.txt",
        ):
            os.replace(tmpdir / name, args.outdir / name)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    print(f"AGR inputs written to: {args.outdir}")
    print(f"Species: {' '.join(selected_species)}")


if __name__ == "__main__":
    main()
