#!/usr/bin/env python3
"""Prepare compact native WGDI 2x2 dot plots for selected chromosome pairs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


COMPARISONS = {
    "pzi_c1a_c5a_vs_gtr_chr02_chr09": {
        "left": ("pzi_hapA1", ("C1A", "C5A")),
        "right": ("gtr", ("Chr02", "Chr09")),
        "blast": "pzi_hapA1_vs_gtr",
    },
    "pzi_c1a_c5a_vs_rch_089_092": {
        "left": ("pzi_hapA1", ("C1A", "C5A")),
        "right": ("rch", ("NC_037089.1", "NC_037092.1")),
        "blast": "pzi_hapA1_vs_rch",
    },
    "pcam_283_285_vs_gtr_chr02_chr09": {
        "left": ("pcam", ("CM114283.1", "CM114285.1")),
        "right": ("gtr", ("Chr02", "Chr09")),
        "blast": "pcam_vs_gtr",
    },
    "pcam_283_285_vs_rch_089_092": {
        "left": ("pcam", ("CM114283.1", "CM114285.1")),
        "right": ("rch", ("NC_037089.1", "NC_037092.1")),
        "blast": "pcam_vs_rch",
    },
    "gtr_chr02_chr09_vs_rch_089_092": {
        "left": ("gtr", ("Chr02", "Chr09")),
        "right": ("rch", ("NC_037089.1", "NC_037092.1")),
        "blast": "gtr_vs_rch",
    },
    "pzi_c1b_vs_gtr_chr04_chr08": {
        "left": ("pzi_hapB1", ("C1B",)),
        "right": ("gtr", ("Chr04", "Chr08")),
        "blast": "pzi_hapB1_vs_gtr",
    },
    "pcam_cm114280_vs_gtr_chr04_chr08": {
        "left": ("pcam", ("CM114280.1",)),
        "right": ("gtr", ("Chr04", "Chr08")),
        "blast": "pcam_vs_gtr",
    },
    "ppa_c4a_vs_gtr_chr04_chr08": {
        "left": ("ppa_hapA", ("C4A",)),
        "right": ("gtr", ("Chr04", "Chr08")),
        "blast": "ppa_hapA_vs_gtr",
    },
    "pzi_c1b_vs_rch_chr091_chr093_chr094": {
        "left": ("pzi_hapB1", ("C1B",)),
        "right": ("rch", ("NC_037091.1", "NC_037093.1", "NC_037094.1")),
        "blast": "pzi_hapB1_vs_rch",
    },
    "pzi_c7a_c4a_vs_c1b": {
        "left": ("pzi_hapA1", ("C7A", "C4A")),
        "right": ("pzi_hapB1", ("C1B",)),
        "blast": "pzi_hapA1_vs_pzi_hapB1",
        "blast_root": "pzi",
    },
    "pcam_cm114280_vs_pzi_c1b": {
        "left": ("pcam", ("CM114280.1",)),
        "right": ("pzi_hapB1", ("C1B",)),
        "blast": "pzi_hapB1_vs_pcam",
        "blast_reverse": "true",
    },
}


def wsl_path(path: Path) -> str:
    resolved = path.resolve()
    if not resolved.drive:
        return resolved.as_posix()
    drive = resolved.drive.rstrip(":").lower()
    return f"/mnt/{drive}{resolved.as_posix()[2:]}"


def read_gff_rows(path: Path, chromosomes: tuple[str, ...]) -> tuple[list[str], set[str]]:
    wanted = set(chromosomes)
    rows: list[str] = []
    gene_ids: set[str] = set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) >= 6 and fields[0] in wanted:
                rows.append(line if line.endswith("\n") else line + "\n")
                gene_ids.add(fields[1])
    if set(chromosomes) - {line.split("\t", 1)[0] for line in rows}:
        missing = sorted(set(chromosomes) - {line.split("\t", 1)[0] for line in rows})
        raise ValueError(f"missing chromosomes in {path}: {', '.join(missing)}")
    return rows, gene_ids


def read_lens_rows(path: Path, chromosomes: tuple[str, ...]) -> list[str]:
    wanted = set(chromosomes)
    by_name: dict[str, str] = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if fields and fields[0] in wanted:
                by_name[fields[0]] = line if line.endswith("\n") else line + "\n"
    missing = [chrom for chrom in chromosomes if chrom not in by_name]
    if missing:
        raise ValueError(f"missing chromosomes in {path}: {', '.join(missing)}")
    return [by_name[chrom] for chrom in chromosomes]


def filter_blast(path: Path, left_ids: set[str], right_ids: set[str], destination: Path) -> int:
    count = 0
    with path.open(encoding="utf-8", errors="replace") as source, destination.open(
        "w", encoding="utf-8", newline=""
    ) as target:
        for line in source:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) >= 2 and fields[0] in left_ids and fields[1] in right_ids:
                target.write(line if line.endswith("\n") else line + "\n")
                count += 1
    return count


def write_config(
    path: Path,
    left_gff: Path,
    left_lens: Path,
    right_gff: Path,
    right_lens: Path,
    blast: Path,
    output: Path,
    left_name: str,
    right_name: str,
    repeat_number: int,
    markersize: float,
    score: float,
    blast_reverse: str,
) -> None:
    values = {
        "gff1": wsl_path(left_gff),
        "gff2": wsl_path(right_gff),
        "lens1": wsl_path(left_lens),
        "lens2": wsl_path(right_lens),
        "blast": wsl_path(blast),
        "genome1_name": left_name,
        "genome2_name": right_name,
        "multiple": 1,
        "score": score,
        "evalue": "1e-5",
        "repeat_number": repeat_number,
        "position": "order",
        "blast_reverse": blast_reverse,
        "ancestor_left": "none",
        "ancestor_top": "none",
        "markersize": markersize,
        "figsize": "8,8",
        "savefig": wsl_path(output),
    }
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("[dotplot]\n")
        for key, value in values.items():
            handle.write(f"{key} = {value}\n")


def read_score_map(path: Path | None) -> dict[str, float]:
    if path is None:
        return {}
    scores: dict[str, float] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            fields = text.split()
            if len(fields) != 2:
                raise ValueError(f"invalid score map line {line_number}: {line.rstrip()}")
            number = float(fields[1])
            scores[fields[0]] = int(number) if number.is_integer() else number
    return scores


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("0903"))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("0903/figures/wgdi_native_2x2_breakpoint_pairs"),
    )
    parser.add_argument("--repeat-number", type=int, default=1)
    parser.add_argument("--markersize", type=float, default=0.35)
    parser.add_argument("--score", type=int, default=300)
    parser.add_argument(
        "--score-map",
        type=Path,
        default=None,
        help="optional two-column figure score map; entries override --score",
    )
    args = parser.parse_args()
    if args.repeat_number < 1:
        raise SystemExit("--repeat-number must be >= 1")

    root = args.root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    score_map = read_score_map(args.score_map)
    prepared = root / "karyotype/prepared"
    blast_root = root / "karyotype/wgdi_native/cross"
    manifest: list[dict[str, object]] = []

    for name, spec in COMPARISONS.items():
        left_label, left_chroms = spec["left"]
        right_label, right_chroms = spec["right"]
        left_dir = prepared / left_label
        right_dir = prepared / right_label
        work = output_dir / "inputs" / name
        work.mkdir(parents=True, exist_ok=True)
        left_gff_rows, left_ids = read_gff_rows(left_dir / f"{left_label}.gff", left_chroms)
        right_gff_rows, right_ids = read_gff_rows(right_dir / f"{right_label}.gff", right_chroms)
        left_gff = work / "left.gff"
        right_gff = work / "right.gff"
        left_lens = work / "left.lens"
        right_lens = work / "right.lens"
        blast = work / "selected.blast.tsv"
        left_gff.write_text("".join(left_gff_rows), encoding="utf-8")
        right_gff.write_text("".join(right_gff_rows), encoding="utf-8")
        left_lens.write_text("".join(read_lens_rows(left_dir / f"{left_label}.lens", left_chroms)), encoding="utf-8")
        right_lens.write_text("".join(read_lens_rows(right_dir / f"{right_label}.lens", right_chroms)), encoding="utf-8")
        source_root = root / "karyotype/wgdi_native" / spec.get("blast_root", "cross")
        source_blast = source_root / spec["blast"] / f"{spec['blast']}.blast.tsv"
        blast_reverse = spec.get("blast_reverse", "false")
        source_left_ids, source_right_ids = left_ids, right_ids
        if blast_reverse == "true":
            source_left_ids, source_right_ids = right_ids, left_ids
        hits = filter_blast(source_blast, source_left_ids, source_right_ids, blast)
        if hits == 0:
            raise SystemExit(f"No selected BLAST hits for {name}")
        output = output_dir / f"{name}.svg"
        config = work / f"{name}.svg.conf"
        figure_score = score_map.get(name, args.score)
        write_config(
            config,
            left_gff,
            left_lens,
            right_gff,
            right_lens,
            blast,
            output,
            left_label,
            right_label,
            args.repeat_number,
            args.markersize,
            figure_score,
            blast_reverse,
        )
        manifest.append(
            {
                "figure": name,
                "left_genome": left_label,
                "left_chromosomes": ",".join(left_chroms),
                "right_genome": right_label,
                "right_chromosomes": ",".join(right_chroms),
                "blast_hits": hits,
                "score": figure_score,
                "repeat_number": args.repeat_number,
                "markersize": args.markersize,
                "blast_reverse": blast_reverse,
                "config": str(config),
                "svg": str(output),
            }
        )

    manifest_path = output_dir / "manifest.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest)
    print(f"Prepared {len(manifest)} native WGDI 2x2 SVG configurations")
    print(manifest_path)


if __name__ == "__main__":
    main()
