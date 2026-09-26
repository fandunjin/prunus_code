#!/usr/bin/env python3
"""Generate one native WGDI dot plot per focal ancestor-unit/chromosome pair.

The plot inputs are restricted to one reference unit and one breakpoint-
associated target chromosome.  No whole-genome concatenation is performed.
Existing 0903 WGDI BLAST tables are filtered by the same GFF gene IDs, then
WGDI's native ``-d`` renderer writes the pairwise SVG/PNG outputs.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import time
from pathlib import Path


TARGETS = (
    ("pzi_hapA1", "PZI", "A"),
    ("psp_hapA1", "PSP", "A"),
    ("ppa_hapA", "PPA", "A"),
    ("pzi_hapB1", "PZI", "B"),
    ("psp_hapB1", "PSP", "B"),
    ("ppa_hapB", "PPA", "B"),
    ("rch", "RCH", "outgroup"),
)

# Retain every extant chromosome that carries one of the six focal
# rearrangement-associated ancestral units.  AK1 (C2B) occurs on A-side
# C5A/C6A/C8A and B-side C6B/C7B; AK3 (GTR03) is represented by A-side C8A
# and B-side C4B.  The remaining chromosomes are omitted from this focused
# breakpoint view.
FOCUS_CHROMS = {
    "A": ("C1A", "C2A", "C3A", "C5A", "C6A", "C7A", "C8A"),
    # C3B is retained because AK1-associated transitions are present on
    # PSP-B C3B in the 0903 junction table; omitting it would hide a
    # supported extant combination from the validation set.
    "B": ("C1B", "C3B", "C4B", "C6B", "C7B", "C8B"),
    "outgroup": ("NC_037088.1", "NC_037089.1", "NC_037090.1", "NC_037091.1", "NC_037092.1", "NC_037093.1", "NC_037094.1"),
}

# AK1=C2B, AK2=C5B and AK3--AK9=GTR03--GTR09.  The present focal display
# retains only the units used by the selected breakpoint pairs.  AK1, AK3,
# AK4, AK6 and AK7 are omitted from this figure.
UNIT_TO_AK = {
    "C2B": "AK1",
    "C5B": "AK2",
    "GTR03": "AK3",
    "GTR05": "AK5",
    "GTR08": "AK8",
    "GTR09": "AK9",
    "GTR04": "AK4",
    "GTR06": "AK6",
    "GTR07": "AK7",
}


def wsl_path(path: Path) -> str:
    """Convert a Windows path under a local drive to the WSL mount path."""
    resolved = path.resolve()
    drive = resolved.drive.rstrip(":").lower()
    return f"/mnt/{drive}{resolved.as_posix()[2:]}"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def gff_rows(path: Path, chromosome: str) -> tuple[list[str], set[str]]:
    rows: list[str] = []
    gene_ids: set[str] = set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 6 or fields[0] != chromosome:
                continue
            rows.append(line if line.endswith("\n") else line + "\n")
            gene_ids.add(fields[1])
    return rows, gene_ids


def lens_row(path: Path, chromosome: str) -> str:
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if fields and fields[0] == chromosome:
                return line if line.endswith("\n") else line + "\n"
    raise ValueError(f"chromosome {chromosome} is missing from {path}")


def filter_blast(path: Path, query_ids: set[str], subject_ids: set[str], destination: Path) -> int:
    """Retain BLAST rows for the selected source and target chromosomes."""
    count = 0
    with path.open(encoding="utf-8", errors="replace") as source, destination.open(
        "w", encoding="utf-8", newline=""
    ) as target:
        for line in source:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 2:
                continue
            if fields[0] in query_ids and fields[1] in subject_ids:
                target.write(line if line.endswith("\n") else line + "\n")
                count += 1
    return count


def write_pair_inputs(
    source_gff: Path,
    source_lens: Path,
    target_gff: Path,
    target_lens: Path,
    blast: Path,
    source_unit: str,
    target_chrom: str,
    pair_dir: Path,
) -> tuple[Path, Path, Path, int]:
    pair_dir.mkdir(parents=True, exist_ok=True)
    source_rows, source_ids = gff_rows(source_gff, source_unit)
    target_rows, target_ids = gff_rows(target_gff, target_chrom)
    if not source_rows or not target_rows:
        raise ValueError(f"missing GFF rows for {source_unit} or {target_chrom}")

    source_gff_out = pair_dir / "source.gff"
    target_gff_out = pair_dir / "target.gff"
    source_lens_out = pair_dir / "source.lens"
    target_lens_out = pair_dir / "target.lens"
    blast_out = pair_dir / "pair.blast.tsv"
    source_gff_out.write_text("".join(source_rows), encoding="utf-8")
    target_gff_out.write_text("".join(target_rows), encoding="utf-8")
    source_lens_out.write_text(lens_row(source_lens, source_unit), encoding="utf-8")
    target_lens_out.write_text(lens_row(target_lens, target_chrom), encoding="utf-8")
    hits = filter_blast(blast, source_ids, target_ids, blast_out)
    return source_gff_out, source_lens_out, target_gff_out, hits


def write_conf(
    path: Path,
    source_gff: Path,
    source_lens: Path,
    target_gff: Path,
    target_lens: Path,
    blast: Path,
    output: Path,
    source_label: str,
    target_label: str,
) -> None:
    values = {
        "gff1": wsl_path(source_gff),
        "gff2": wsl_path(target_gff),
        "lens1": wsl_path(source_lens),
        "lens2": wsl_path(target_lens),
        "blast": wsl_path(blast),
        "genome1_name": source_label,
        "genome2_name": target_label,
        "multiple": 1,
        "score": 100,
        "evalue": "1e-5",
        "repeat_number": 10,
        "position": "order",
        "blast_reverse": "false",
        "ancestor_left": "none",
        "ancestor_top": "none",
        "markersize": 0.5,
        "figsize": "10,10",
        "savefig": wsl_path(output),
    }
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("[dotplot]\n")
        for key, value in values.items():
            handle.write(f"{key} = {value}\n")


def run_wgdi(config: Path, attempts: int = 3) -> None:
    command = [
        # Keep Matplotlib non-interactive.  Passing empty DISPLAY/Qt variables
        # can make WSL fail with E_ACCESSDENIED before WGDI starts.
        "wsl", "--exec", "env", "MPLBACKEND=Agg",
        "/opt/plant-wakr-venv/bin/wgdi", "-d", wsl_path(config),
    ]
    last_error: subprocess.CalledProcessError | None = None
    for attempt in range(1, attempts + 1):
        try:
            subprocess.run(command, check=True)
            return
        except subprocess.CalledProcessError as error:
            last_error = error
            if attempt < attempts:
                time.sleep(2 * attempt)
    assert last_error is not None
    raise last_error


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("0903"))
    parser.add_argument("--output-dir", type=Path, default=Path("0903/figures/fig4_breakpoint_native_pairs"))
    parser.add_argument("--prepare-only", action="store_true",
                        help="write filtered inputs/configs without invoking WSL/WGDI")
    parser.add_argument("--targets", default="all",
                        help="comma-separated target labels to process, or all")
    parser.add_argument("--units", default="all",
                        help="comma-separated reference units (for example C2B,C5B,GTR03), or all")
    args = parser.parse_args()

    root = args.root.resolve()
    out = args.output_dir.resolve()
    pairs_dir = out / "pairs"
    out.mkdir(parents=True, exist_ok=True)
    pairs_dir.mkdir(parents=True, exist_ok=True)

    reference_gff = root / "karyotype/fig4_breakpoint_candidates/reference/prunus9_reference.gff"
    reference_lens = root / "karyotype/fig4_breakpoint_candidates/reference/prunus9_reference.lens"
    blocks_path = root / "karyotype/fig4_breakpoint_candidates/selected_primary_blocks.tsv"
    mapping_root = root / "karyotype/fig4_breakpoint_candidates/mapping"
    prepared_root = root / "karyotype/prepared"

    blocks = read_tsv(blocks_path)
    requested_targets = {item.strip() for item in args.targets.split(",") if item.strip()}
    target_records = TARGETS if args.targets.strip().lower() == "all" else tuple(item for item in TARGETS if item[0] in requested_targets)
    if not target_records:
        raise SystemExit(f"No target labels matched --targets={args.targets}")
    requested_units = set(UNIT_TO_AK) if args.units.strip().lower() == "all" else {item.strip() for item in args.units.split(",") if item.strip()}
    unknown_units = requested_units - set(UNIT_TO_AK)
    if unknown_units:
        raise SystemExit(f"Unknown reference units: {', '.join(sorted(unknown_units))}")
    manifest_rows: list[dict[str, object]] = []
    for target, species, subgenome in target_records:
        target_gff = prepared_root / target / f"{target}.gff"
        target_lens = prepared_root / target / f"{target}.lens"
        blast = mapping_root / target / f"prunus9_vs_{target}.blast.tsv"
        target_focus = set(FOCUS_CHROMS[subgenome])
        selected = [
            row for row in blocks
            if row.get("target_label") == target
            and row.get("source_unit") in requested_units
            and row.get("target_chromosome") in target_focus
        ]
        pair_keys = sorted({(row["source_unit"], row["target_chromosome"]) for row in selected})
        for source_unit, target_chrom in pair_keys:
            pair_name = safe_name(f"{target}__{UNIT_TO_AK[source_unit]}_{target_chrom}")
            pair_dir = pairs_dir / pair_name
            try:
                source_gff, source_lens, target_gff_one, hits = write_pair_inputs(
                    reference_gff,
                    reference_lens,
                    target_gff,
                    target_lens,
                    blast,
                    source_unit,
                    target_chrom,
                    pair_dir,
                )
            except ValueError:
                continue
            if hits == 0:
                continue
            source_label = UNIT_TO_AK[source_unit]
            target_label = target_chrom
            outputs: dict[str, Path] = {}
            for extension in ("svg", "png"):
                output = out / f"{pair_name}.{extension}"
                config = pair_dir / f"{pair_name}.{extension}.conf"
                write_conf(config, source_gff, source_lens, target_gff_one,
                           pair_dir / "target.lens", pair_dir / "pair.blast.tsv",
                           output, source_label, target_label)
                if not args.prepare_only:
                    run_wgdi(config)
                outputs[extension] = output
            event_context = (
                "A-associated" if subgenome == "A"
                else "B-associated" if subgenome == "B"
                else "outgroup"
            )
            manifest_rows.append({
                "pair": pair_name,
                "species": species,
                "target": target,
                "subgenome": subgenome,
                "event_context": event_context,
                "ancestor_unit": UNIT_TO_AK[source_unit],
                "reference_unit": source_unit,
                "target_chromosome": target_chrom,
                "selected_primary_blocks": sum(
                    row["source_unit"] == source_unit and row["target_chromosome"] == target_chrom
                    for row in selected
                ),
                "blast_hits": hits,
                "svg": str(outputs["svg"]),
                "png": str(outputs["png"]),
            })

    manifest = out / "pair_manifest.tsv"
    fields = list(manifest_rows[0]) if manifest_rows else [
        "pair", "species", "target", "subgenome", "event_context", "ancestor_unit",
        "reference_unit", "target_chromosome", "selected_primary_blocks", "blast_hits", "svg", "png",
    ]
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest_rows)
    print(f"Generated {len(manifest_rows)} native WGDI pair plots")
    print(manifest)


if __name__ == "__main__":
    main()
