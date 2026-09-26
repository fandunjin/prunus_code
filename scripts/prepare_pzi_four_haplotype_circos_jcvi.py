#!/usr/bin/env python3
"""Prepare native Circos tracks and JCVI block links for four PZI haplotypes."""

from __future__ import annotations

import argparse
import gzip
import math
import os
import re
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Iterable, TextIO


HAPS = ("A1", "A2", "B1", "B2")
# Restore the original all-A2/all-B1 display convention outside C01/C02,
# while preserving the requested C01/C02 directions: C01_B1 forward,
# C01_B2 reverse, C02_B1 reverse and C02_B2 forward.  Relative to that layout,
# C08_B1 is toggled back to forward and C08_B2 is toggled to reverse.  The
# operation does not alter FASTA, GFF3, BED or JCVI coordinates.
REVERSED_DISPLAY_CHROMS = frozenset(
    [f"C{i:02d}_A2" for i in range(1, 9)]
    + [f"C{i:02d}_B1" for i in range(2, 8)]
    + ["C01_B2", "C08_B2"]
)
HAP_COLORS = {
    "A1": "hapa1",
    "A2": "hapa2",
    "B1": "hapb1",
    "B2": "hapb2",
}
# Source chromosome IDs are ``chr01_A1`` in the original package and
# ``C01_A1`` in the renamed downstream package.  Both encode the same
# quartet number and haplotype.
CHR_RE = re.compile(r"^(?:chr|C)(\d{2})_(A1|A2|B1|B2)$")
DISPLAY_RE = re.compile(r"^C(\d{2})_(A1|A2|B1|B2)$")
QUARTET_LABELS = {
    1: "C01",
    2: "C02",
    3: "C03",
    4: "C04",
    5: "C05",
    6: "C06",
    7: "C07",
    8: "C08",
}
PAIR_LABELS = tuple(QUARTET_LABELS[number] for number in range(1, 9))
PAIR_COLORS_A1B1 = {label: f"{label.lower()}_ab1" for label in PAIR_LABELS}
PAIR_COLORS_A2B2 = {label: f"{label.lower()}_ab2" for label in PAIR_LABELS}
# Manual alias correction requested after inspection of the A-B self-JCVI
# pattern.  The corresponding A1 and A2 identities are exchanged together;
# all source sequence IDs and coordinates remain unchanged.
DISPLAY_ALIAS_OVERRIDES = {
    "C01_A1": "C06_A1",
    "C06_A1": "C01_A1",
    "C01_A2": "C06_A2",
    "C06_A2": "C01_A2",
}
# Okabe-Ito-inspired qualitative colors provide red-green color-vision
# accessibility and restrained publication styling.  Matched A1-B1/A2-B2
# variants have comparable luminance with only a small hue shift.
PAIR_RGB_A1B1 = {
    "C01": "0,114,178",    # blue
    "C02": "230,159,0",    # orange
    "C03": "0,158,115",    # bluish green
    "C04": "213,94,0",     # vermilion
    "C05": "204,121,167",  # reddish purple
    "C06": "86,180,233",   # sky blue
    "C07": "190,168,0",    # mustard yellow
    "C08": "80,80,80",     # charcoal
}
PAIR_RGB_A2B2 = {
    "C01": "38,105,176",   # blue variant
    "C02": "215,150,18",   # orange variant
    "C03": "28,148,116",   # green variant
    "C04": "201,96,21",    # vermilion variant
    "C05": "190,122,173",  # purple variant
    "C06": "96,170,221",   # sky-blue variant
    "C07": "177,156,18",   # mustard variant
    "C08": "96,96,96",     # charcoal variant
}
# Circos only accepts integer _aN suffixes.  Use 100 alpha steps so the
# requested native-scale a5.1 can be represented without the decimal suffix.
# After visual review, use 65% transparency so the ribbons are stronger.  This
# keeps PNG and SVG/PDF rendering consistent; a
# four-component RGBA definition is not propagated correctly by Circos' SVG
# writer.
AUTO_ALPHA_STEPS = 100
HIGHLIGHT_ALPHA_INDEX = 65
BACKGROUND_ALPHA_INDEX = 50
HIGHLIGHT_ALPHA = f"a{HIGHLIGHT_ALPHA_INDEX}"
BACKGROUND_ALPHA = f"a{BACKGROUND_ALPHA_INDEX}"
HIGHLIGHT_TRANSPARENCY = HIGHLIGHT_ALPHA_INDEX / AUTO_ALPHA_STEPS
COPY1_PAIR_ROWS = (
    ("C01", "06", "15", "06", "15", "manual_A_C01_C06_swap"),
    ("C02", "02", "11", "02", "11", "rearranged_primary"),
    ("C03", "03", "12", "03", "12", "rearranged_primary"),
    ("C04", "04", "07", "04", "07", "fixed_backbone"),
    ("C05", "05", "16", "05", "16", "fixed_backbone"),
    ("C06", "01", "13", "01", "13", "manual_A_C01_C06_swap"),
    ("C07", "08", "10", "08", "10", "fixed_backbone"),
    ("C08", "09", "14", "09", "14", "fixed_backbone"),
)
LINK_SPECS = (
    ("A1", "A2", "background"),
    ("B1", "B2", "background"),
    ("A1", "B1", "highlight"),
    ("A2", "B2", "highlight"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument(
        "--hap-rename-dir",
        type=Path,
        help="Renamed package root containing genome/gff/cds/pep subdirectories",
    )
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--jcvi-dir", type=Path)
    parser.add_argument("--window", type=int, default=1_000_000)
    parser.add_argument(
        "--stage-jcvi-inputs",
        action="store_true",
        help="Create hap*.bed and hap*.pep inputs, then exit before track generation",
    )
    return parser.parse_args()


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def require_file(path: Path, label: str) -> Path:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty {label}: {path}")
    return path


def first_existing(paths: Iterable[Path], label: str) -> Path:
    candidates = list(paths)
    for path in candidates:
        if path.is_file() and path.stat().st_size > 0:
            return path
    joined = "\n  ".join(str(path) for path in candidates)
    raise FileNotFoundError(f"No usable {label} found. Checked:\n  {joined}")


def renamed_package_path(root: Path, package_dir: Path | None = None) -> Path:
    return (package_dir or root / "haploid_rename").resolve()


def renamed_package_available(root: Path, package_dir: Path | None = None) -> bool:
    """Return true only when the complete renamed four-haplotype package exists."""
    package = renamed_package_path(root, package_dir)
    return all(
        (package / subdir / f"pzi_hap{hap}.{suffix}").is_file()
        and (package / subdir / f"pzi_hap{hap}.{suffix}").stat().st_size > 0
        for hap in HAPS
        for subdir, suffix in (
            ("genome", "fa"),
            ("gff", "gff3"),
            ("cds", "cds.fa"),
            ("pep", "pep.fa"),
        )
    )


def renamed_seqid_aliases(
    root: Path, package_dir: Path | None = None
) -> dict[str, str]:
    """Map pre-package compact chromosome IDs to final Cxx_hap aliases."""
    package_map = renamed_package_path(root, package_dir) / "chromosome_rename_map.tsv"
    if not package_map.is_file() or package_map.stat().st_size == 0:
        return {}
    aliases: dict[str, str] = {}
    with open_text(package_map) as handle:
        header_line = handle.readline().rstrip("\n")
        header = header_line.split("\t")
        required = {"source_seqid", "display_seqid"}
        if not required.issubset(header):
            raise ValueError(
                f"Renamed package map lacks columns {sorted(required - set(header))}: {package_map}"
            )
        index = {name: i for i, name in enumerate(header)}
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < len(header):
                continue
            source = fields[index["source_seqid"]]
            display = fields[index["display_seqid"]]
            if source in aliases and aliases[source] != display:
                raise ValueError(f"Conflicting renamed package aliases for {source}")
            aliases[source] = display
    return aliases


def genome_path(root: Path, hap: str, package_dir: Path | None = None) -> Path:
    renamed_dir = renamed_package_path(root, package_dir)
    return first_existing(
        (
            renamed_dir / "genome" / f"pzi_hap{hap}.fa",
            root / "haploid" / "final_haplotypes" / f"pzi_hap{hap}.fa",
        ),
        f"hap{hap} genome",
    )


def annotation_paths(
    root: Path, hap: str, package_dir: Path | None = None
) -> tuple[Path, Path]:
    label = f"pzi_hap{hap}"
    renamed_dir = renamed_package_path(root, package_dir)
    canonical_dir = root / "downstream_analysis" / "work" / "canonical" / label
    subgenome_prefix = root / "downstream_inputs" / "subgenomes" / label
    package_dir = root / "final_haplotype_genomes_eviann"
    primary_dir = root / "annotation" / "genes" / "primary_isoform_busco" / "primary"
    post_dir = root / "annotation" / "genes" / "final_haplotypes_eviann" / f"hap{hap}" / "postprocess"
    candidates = (
        (
            renamed_dir / "gff" / f"{label}.gff3",
            renamed_dir / "pep" / f"{label}.pep.fa",
            renamed_dir / "cds" / f"{label}.cds.fa",
            "renamed haploid package",
        ),
        (
            canonical_dir / "primary.gff3",
            canonical_dir / "primary.pep.fa",
            canonical_dir / "primary.cds.fa",
            "downstream canonical (same as previous JCVI)",
        ),
        (
            Path(f"{subgenome_prefix}.primary.gff3"),
            Path(f"{subgenome_prefix}.primary.pep.fa"),
            Path(f"{subgenome_prefix}.primary.cds.fa"),
            "downstream subgenomes (same as previous JCVI fallback)",
        ),
        (
            package_dir / f"{label}.gff3",
            package_dir / f"{label}.pep.fa",
            package_dir / f"{label}.cds.fa",
            "collected EviAnn package (same as previous JCVI fallback)",
        ),
        (
            primary_dir / f"{label}.primary.gff3",
            primary_dir / f"{label}.primary.pep.fa",
            primary_dir / f"{label}.primary.cds.fa",
            "primary_isoform_busco",
        ),
        (
            post_dir / f"{label}.eviann.genes.gff3",
            post_dir / f"{label}.eviann.protein.fa",
            post_dir / f"{label}.eviann.cds.fa",
            "EviAnn postprocess",
        ),
    )
    for gff, pep, cds, source in candidates:
        canonical_requires_cds = source.startswith("downstream canonical")
        cds_ok = cds.is_file() and cds.stat().st_size > 0
        if (
            gff.is_file()
            and gff.stat().st_size > 0
            and pep.is_file()
            and pep.stat().st_size > 0
            and (cds_ok if canonical_requires_cds else True)
        ):
            print(f"[INFO] hap{hap} annotation source: {source}")
            return gff, pep
    checked = "\n".join(f"  {source}: {gff} | {pep} | {cds}" for gff, pep, cds, source in candidates)
    raise FileNotFoundError(
        f"No matched GFF3/protein pair found for hap{hap}. Checked:\n{checked}"
    )


def gene_gff_path(root: Path, hap: str, package_dir: Path | None = None) -> Path:
    return annotation_paths(root, hap, package_dir)[0]


def peptide_path(root: Path, hap: str, package_dir: Path | None = None) -> Path:
    return annotation_paths(root, hap, package_dir)[1]


def te_gff_path(root: Path, hap: str) -> Path:
    filename = f"pzi_hap{hap}.fa.mod.EDTA.TEanno.gff3"
    return first_existing(
        (
            root / "annotation" / "repeats" / "final_haplotypes" / "panEDTA_uniform" / filename,
            root / "annotation" / "repeats" / "final_haplotypes" / "uniform" / filename,
            root / "annotation" / "repeats" / "final_haplotypes" / f"hap{hap}_edta" / filename,
        ),
        f"hap{hap} EDTA GFF3",
    )


def parse_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attributes[key] = value
    return attributes


def load_quartet_assignment(
    root: Path,
    package_dir: Path | None = None,
) -> tuple[dict[str, str], dict[str, str], dict[int, tuple[str, str]]]:
    """Load the authoritative quartet assignment and build Circos-only aliases.

    Legacy inputs retain source IDs such as ``chr02_A1`` and receive Circos
    aliases at render time.  Inputs from ``haploid_rename`` already use final
    ``C01_A1`` aliases, so those IDs are kept unchanged and their provenance is
    recovered from ``chromosome_rename_map.tsv``.
    """
    path = root / "haploid" / "final_haplotypes" / "details" / "final_quartet_assignment.tsv"
    require_file(path, "final quartet assignment table")
    source_to_display: dict[str, str] = {}
    display_to_source: dict[str, str] = {}
    quartet_meta: dict[int, tuple[str, str]] = {}
    renamed_dir = renamed_package_path(root, package_dir)
    package_mode = renamed_package_available(root, renamed_dir)
    package_aliases: dict[tuple[str, str], str] = {}
    if package_mode:
        package_map = renamed_dir / "chromosome_rename_map.tsv"
        require_file(package_map, "renamed package chromosome map")
        with open_text(package_map) as package_handle:
            package_header_line = package_handle.readline().rstrip("\n")
            package_header = package_header_line.split("\t")
            package_required = {"original_chr32", "display_seqid", "haplotype"}
            if not package_required.issubset(package_header):
                raise ValueError(
                    f"Renamed package map lacks columns {sorted(package_required - set(package_header))}: {package_map}"
                )
            package_index = {name: i for i, name in enumerate(package_header)}
            for raw in package_handle:
                if not raw.strip() or raw.startswith("#"):
                    continue
                fields = raw.rstrip("\n").split("\t")
                if len(fields) < len(package_header):
                    continue
                hap = fields[package_index["haplotype"]]
                if hap in HAPS:
                    package_aliases[(fields[package_index["original_chr32"]], hap)] = fields[
                        package_index["display_seqid"]
                    ]
    # The assignment table records the pre-compaction chromosome IDs
    # (chr01..chr32).  Final haplotype FASTAs use compact chr01..chr08 IDs;
    # resolve that translation when the authoritative rename map is present.
    rename_path = root / "haploid" / "final_haplotypes" / "pzi_A1A2B1B2.rename_map.tsv"
    old_to_new: dict[tuple[str, str], str] = {}
    if rename_path.is_file() and rename_path.stat().st_size > 0:
        with open_text(rename_path) as rename_handle:
            rename_header_line = rename_handle.readline().rstrip("\n")
            rename_header = rename_header_line.split("\t")
            rename_required = {"old_chr", "new_chr", "hap_set"}
            if rename_required.issubset(rename_header):
                rename_index = {name: i for i, name in enumerate(rename_header)}
                for raw in rename_handle:
                    if not raw.strip() or raw.startswith("#"):
                        continue
                    fields = raw.rstrip("\n").split("\t")
                    if len(fields) < len(rename_header):
                        continue
                    old_chr = fields[rename_index["old_chr"]]
                    new_chr = fields[rename_index["new_chr"]]
                    hap_set = fields[rename_index["hap_set"]]
                    if hap_set in HAPS:
                        old_to_new[(old_chr, hap_set)] = new_chr
    required = {"quartet_id", "hapA1_chr", "hapA2_chr", "hapB1_chr", "hapB2_chr", "confidence"}
    with open_text(path) as handle:
        header_line = handle.readline().rstrip("\n")
        if not header_line:
            raise ValueError(f"Empty quartet assignment table: {path}")
        header = header_line.split("\t")
        missing_columns = required - set(header)
        if missing_columns:
            raise ValueError(
                f"Quartet assignment table lacks columns {sorted(missing_columns)}: {path}"
            )
        index = {name: i for i, name in enumerate(header)}
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < len(header):
                raise ValueError(f"Malformed quartet assignment row in {path}: {raw.rstrip()}")
            quartet_id = fields[index["quartet_id"]]
            match = re.search(r"(?:final_)?quartet(\d+)$", quartet_id)
            if match is None:
                raise ValueError(f"Cannot parse quartet number from {quartet_id!r} in {path}")
            number = int(match.group(1))
            if number in quartet_meta:
                raise ValueError(f"Duplicate quartet number q{number:02d} in {path}")
            confidence = fields[index["confidence"]]
            quartet_meta[number] = (quartet_id, confidence)
            for hap in HAPS:
                assignment_chr = fields[index[f"hap{hap}_chr"]]
                nominal_display = f"{QUARTET_LABELS[number]}_{hap}"
                if package_mode:
                    # FASTA/GFF records in haploid_rename already carry their
                    # final Cxx_hap aliases.  Map the original chr32 assignment
                    # to that alias, and keep the rendered alias unchanged.
                    source = package_aliases.get((assignment_chr, hap))
                    if source is None:
                        raise ValueError(
                            f"Renamed package map has no entry for {assignment_chr}/{hap}"
                        )
                    display = source
                else:
                    source = old_to_new.get((assignment_chr, hap))
                    if source is None:
                        # Final files are normally compacted in quartet order.
                        source = f"chr{number:02d}_{hap}"
                    display = DISPLAY_ALIAS_OVERRIDES.get(nominal_display, nominal_display)
                if source in source_to_display:
                    raise ValueError(f"Source chromosome assigned more than once: {source}")
                if display in display_to_source:
                    raise ValueError(f"Display chromosome assigned more than once: {display}")
                source_to_display[source] = display
                display_to_source[display] = source

    expected_prefix = "C" if package_mode else "chr"
    expected = {f"{expected_prefix}{i:02d}_{hap}" for hap in HAPS for i in range(1, 9)}
    missing = expected - set(source_to_display)
    extra = set(source_to_display) - expected
    if missing or extra:
        raise ValueError(
            "Quartet assignment does not cover the expected chromosomes; "
            f"missing={sorted(missing)}, extra={sorted(extra)}"
        )
    expected_display = {f"{QUARTET_LABELS[i]}_{hap}" for i in range(1, 9) for hap in HAPS}
    missing_display = expected_display - set(display_to_source)
    if missing_display:
        raise ValueError(f"Quartet assignment is missing aliases: {sorted(missing_display)}")
    return source_to_display, display_to_source, quartet_meta


def stage_jcvi_inputs(
    root: Path, out_dir: Path, package_dir: Path | None = None
) -> None:
    input_dir = out_dir / "work" / "jcvi_inputs"
    input_dir.mkdir(parents=True, exist_ok=True)
    manifest = input_dir / "input_manifest.tsv"
    with manifest.open("w", encoding="utf-8") as report:
        report.write("haplotype\tbed\tproteins\tgene_count\n")
        for hap in HAPS:
            gff, pep = annotation_paths(root, hap, package_dir)
            records: list[tuple[str, int, int, str, str]] = []
            seen: set[str] = set()
            with open_text(gff) as handle:
                for raw in handle:
                    if raw.startswith("#"):
                        continue
                    fields = raw.rstrip("\n").split("\t")
                    if len(fields) < 9 or fields[2] not in {"mRNA", "transcript"}:
                        continue
                    attrs = parse_attributes(fields[8])
                    if attrs.get("primary_isoform", "true").lower() not in {"true", "1", "yes"}:
                        continue
                    feature_id = attrs.get("ID")
                    if not feature_id or feature_id in seen:
                        continue
                    seen.add(feature_id)
                    records.append(
                        (fields[0], int(fields[3]) - 1, int(fields[4]), feature_id, fields[6])
                    )
            records.sort(key=lambda row: (row[0], row[1], row[2], row[3]))
            if not records:
                raise ValueError(f"No primary transcripts parsed from {gff}")
            bed = input_dir / f"hap{hap}.bed"
            with bed.open("w", encoding="utf-8") as handle:
                for chrom, start, end, feature_id, strand in records:
                    handle.write(f"{chrom}\t{start}\t{end}\t{feature_id}\t0\t{strand}\n")
            pep_link = input_dir / f"hap{hap}.pep"
            if pep_link.is_symlink() or pep_link.exists():
                pep_link.unlink()
            try:
                os.symlink(pep.resolve(), pep_link)
            except OSError:
                # Windows commonly disallows symlink creation without elevated
                # privileges. HPC/Linux runs use symlinks; this fallback keeps
                # local validation possible.
                shutil.copy2(pep, pep_link)
            report.write(f"hap{hap}\t{bed}\t{pep}\t{len(records)}\n")
            print(f"[INFO] staged hap{hap}: {len(records)} genes")
    (input_dir / ".staged.ok").write_text("JCVI inputs staged\n", encoding="ascii")


def fasta_windows(path: Path, window: int) -> tuple[dict[str, int], dict[str, list[tuple[int, int, float]]]]:
    sizes: dict[str, int] = {}
    windows: dict[str, list[tuple[int, int, float]]] = defaultdict(list)
    name: str | None = None
    position = 0
    win_start = 0
    win_gc = 0
    win_called = 0
    win_bases = 0

    def flush_window() -> None:
        nonlocal win_start, win_gc, win_called, win_bases
        if name is not None and win_bases:
            gc = win_gc / win_called if win_called else 0.0
            windows[name].append((win_start, win_start + win_bases, gc))
            win_start += win_bases
            win_gc = win_called = win_bases = 0

    def flush_sequence() -> None:
        nonlocal position
        if name is not None:
            flush_window()
            sizes[name] = position
        position = 0

    with open_text(path) as handle:
        for raw in handle:
            if raw.startswith(">"):
                flush_sequence()
                name = raw[1:].strip().split()[0]
                if name in sizes or name in windows:
                    raise ValueError(f"Duplicate FASTA sequence ID in {path}: {name}")
                win_start = 0
                continue
            if name is None:
                continue
            seq = raw.strip().upper()
            offset = 0
            while offset < len(seq):
                take = min(window - win_bases, len(seq) - offset)
                chunk = seq[offset : offset + take]
                win_gc += chunk.count("G") + chunk.count("C")
                win_called += sum(chunk.count(base) for base in "ACGT")
                win_bases += take
                position += take
                offset += take
                if win_bases == window:
                    flush_window()
    flush_sequence()
    return sizes, windows


def validate_chromosomes(sizes: dict[str, int], hap: str, renamed: bool = False) -> None:
    prefix = "C" if renamed else "chr"
    expected = {f"{prefix}{i:02d}_{hap}" for i in range(1, 9)}
    missing = expected - set(sizes)
    extra = set(sizes) - expected
    if missing or extra:
        raise ValueError(
            f"Unexpected hap{hap} chromosomes; missing={sorted(missing)}, extra={sorted(extra)}"
        )


def parse_gene_counts(path: Path, sizes: dict[str, int], window: int) -> dict[str, list[int]]:
    counts = {chrom: [0] * math.ceil(length / window) for chrom, length in sizes.items()}
    with open_text(path) as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene" or fields[0] not in counts:
                continue
            midpoint = ((int(fields[3]) - 1) + int(fields[4])) // 2
            index = min(midpoint // window, len(counts[fields[0]]) - 1)
            counts[fields[0]][index] += 1
    return counts


def parse_te_intervals(
    path: Path,
    sizes: dict[str, int],
    seqid_aliases: dict[str, str] | None = None,
) -> dict[str, dict[str, list[tuple[int, int]]]]:
    classes = {
        "all": defaultdict(list),
        "gypsy": defaultdict(list),
        "copia": defaultdict(list),
    }
    with open_text(path) as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom = (seqid_aliases or {}).get(fields[0], fields[0])
            if chrom not in sizes:
                continue
            start = max(0, int(fields[3]) - 1)
            end = min(sizes[chrom], int(fields[4]))
            if end <= start:
                continue
            classes["all"][chrom].append((start, end))
            attrs = parse_attributes(fields[8])
            classification = attrs.get("classification", "").lower()
            feature_type = fields[2].lower()
            if "gypsy" in classification or "gypsy" in feature_type:
                classes["gypsy"][chrom].append((start, end))
            if "copia" in classification or "copia" in feature_type:
                classes["copia"][chrom].append((start, end))
    return classes


def window_coverage(
    intervals: dict[str, list[tuple[int, int]]], sizes: dict[str, int], window: int
) -> dict[str, list[int]]:
    coverage = {chrom: [0] * math.ceil(length / window) for chrom, length in sizes.items()}
    for chrom, chrom_intervals in intervals.items():
        merged: list[list[int]] = []
        for start, end in sorted(chrom_intervals):
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        for start, end in merged:
            for index in range(start // window, (end - 1) // window + 1):
                win_start = index * window
                win_end = min((index + 1) * window, sizes[chrom])
                coverage[chrom][index] += max(0, min(end, win_end) - max(start, win_start))
    return coverage


def percentile(values: list[float], fraction: float) -> float:
    if not values:
        return 1.0
    ordered = sorted(values)
    index = min(len(ordered) - 1, max(0, round((len(ordered) - 1) * fraction)))
    return ordered[index]


def coverage_plot_max(values: list[float]) -> float:
    """Choose a readable coverage scale while limiting only extreme outliers."""
    if not values:
        return 0.1
    observed = max(values)
    robust = percentile(values, 0.995) * 1.10
    # Keep a small headroom above the visible range; values above this scale
    # are rare windows and are clipped by Circos rather than flattening the plot.
    return min(1.0, max(0.02, robust, observed * 0.98))


def load_bed(path: Path) -> dict[str, tuple[str, int, int]]:
    records: dict[str, tuple[str, int, int]] = {}
    with open_text(path) as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            records[fields[3]] = (fields[0], int(fields[1]), int(fields[2]))
    if not records:
        raise ValueError(f"No BED records parsed from {path}")
    return records


def block_interval(
    first_id: str,
    last_id: str,
    bed: dict[str, tuple[str, int, int]],
) -> tuple[str, int, int]:
    if first_id not in bed or last_id not in bed:
        missing = first_id if first_id not in bed else last_id
        raise KeyError(f"JCVI anchor gene is absent from BED: {missing}")
    first = bed[first_id]
    last = bed[last_id]
    if first[0] != last[0]:
        raise ValueError(f"JCVI simple block crosses chromosomes: {first_id}, {last_id}")
    return first[0], min(first[1], last[1]), max(first[2], last[2])


def parse_jcvi_blocks(
    simple_path: Path,
    left_bed: dict[str, tuple[str, int, int]],
    right_bed: dict[str, tuple[str, int, int]],
) -> list[tuple[str, int, int, str, int, int, str, str]]:
    blocks: list[tuple[str, int, int, str, int, int, str, str]] = []
    with open_text(simple_path) as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.split()
            if len(fields) != 6:
                raise ValueError(f"Expected six columns in {simple_path}: {raw.rstrip()}")
            left = block_interval(fields[0], fields[1], left_bed)
            right = block_interval(fields[2], fields[3], right_bed)
            blocks.append((*left, *right, fields[4], fields[5]))
    if not blocks:
        raise ValueError(f"No JCVI blocks parsed from {simple_path}")
    return blocks


def chromosome_layout(
    sizes: dict[str, int],
    source_to_display: dict[str, str],
    display_to_source: dict[str, str],
) -> tuple[list[str], dict[str, float]]:
    # Clockwise from 12 o'clock: B1 upper-right, B2 lower-right,
    # A2 lower-left, A1 upper-left.
    order = (
        [f"{label}_B1" for label in PAIR_LABELS]
        + [f"{label}_B2" for label in reversed(PAIR_LABELS)]
        + [f"{label}_A2" for label in PAIR_LABELS]
        + [f"{label}_A1" for label in reversed(PAIR_LABELS)]
    )
    if set(order) != set(display_to_source):
        raise ValueError("Chromosome layout aliases do not match the quartet assignment table")
    totals = {hap: 0 for hap in HAPS}
    for source, display in source_to_display.items():
        if source not in sizes:
            raise ValueError(f"Quartet assignment chromosome is absent from FASTA: {source}")
        match = DISPLAY_RE.fullmatch(display)
        if match is None:
            raise ValueError(f"Invalid display chromosome alias: {display}")
        totals[match.group(2)] += sizes[source]
    target = sum(totals.values()) / len(totals)
    factors = {hap: target / totals[hap] for hap in HAPS}
    return order, factors


def display_haplotype(chromosome: str) -> str:
    """Return the haplotype suffix from a Circos display chromosome."""
    match = DISPLAY_RE.fullmatch(chromosome)
    if match is None:
        raise ValueError(f"Invalid display chromosome alias: {chromosome}")
    return match.group(2)


def link_display_orientation(
    raw_orientation: str,
    left_display: str,
    right_display: str,
) -> str:
    """Return the orientation used by a displayed Circos link.

    The endpoint reversal must be XORed with JCVI's native orientation for
    gray links.  Reversal is determined from the explicit chromosome-level
    display set.  This changes only the rendered direction; source coordinates
    and true local inversions remain unchanged.
    """
    if raw_orientation not in {"+", "-"}:
        raise ValueError(f"Unexpected JCVI block orientation: {raw_orientation!r}")
    inverted = raw_orientation == "-"
    inverted ^= left_display in REVERSED_DISPLAY_CHROMS
    inverted ^= right_display in REVERSED_DISPLAY_CHROMS
    return "-" if inverted else "+"


def write_circos_config(
    conf_dir: Path,
    chromosome_order: list[str],
    scale_factors: dict[str, float],
    gene_max: float,
    gc_min: float,
    gc_max: float,
    gypsy_max: float,
    copia_max: float,
) -> None:
    chromosomes = ";".join(chromosome_order)
    scale_items = [
        f"{chrom}={scale_factors[DISPLAY_RE.fullmatch(chrom).group(2)]:.8f}"
        for chrom in chromosome_order
    ]
    scales = ",".join(scale_items)
    reversed_chromosomes = ",".join(
        chrom for chrom in chromosome_order if chrom in REVERSED_DISPLAY_CHROMS
    )
    pair_color_defs = "\n".join(
        f"{PAIR_COLORS_A1B1[label]} = {PAIR_RGB_A1B1[label]}\n"
        f"{PAIR_COLORS_A2B2[label]} = {PAIR_RGB_A2B2[label]}"
        for label in PAIR_LABELS
    )
    config = f"""<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>

<colors>
hapa1        = 55,95,145
hapa2        = 120,168,205
hapb1        = 186,85,40
hapb2        = 222,151,55
gene_track   = 85,90,95
te_track     = 112,135,52
gypsy_track  = 196,91,82
copia_track  = 76,120,168
gc_track     = 42,88,108
within_gray  = 175,175,175
{pair_color_defs}
</colors>

karyotype = data/karyotype.txt
chromosomes_units = 1000000
chromosomes_display_default = no
chromosomes = {chromosomes}
chromosomes_scale = {scales}
# Reverse only the explicitly requested complete ideograms.  The underlying
# FASTA, GFF3, BED and JCVI block coordinates remain unchanged.
chromosomes_reverse = {reversed_chromosomes}

<image>
background = white
dir   = figures
file  = pzi_four_haplotypes_jcvi_circos.png
png   = yes
svg   = yes
radius = 1800p
# Shift by half the top A/B gap so the approximately 10-degree opening is
# centered at 12 o'clock.
angle_offset = -85.0
auto_alpha_colors = yes
auto_alpha_steps  = {AUTO_ALPHA_STEPS}
24bit             = yes
</image>

<ideogram>
<spacing>
default = 0.0025r
<pairwise C01_A1 C01_B1>
# Pairwise r values are multiples of default spacing; with one opening,
# 12.4r gives approximately 10 degrees.  Every other boundary uses default.
spacing = 12.4r
</pairwise>
</spacing>
# Move only the chromosome ring outward.  The plot radii below are compensated
# for the change in the ideogram inner radius so the inner track geometry stays
# at approximately its previous absolute canvas position.
radius           = 0.940r
thickness        = 22p
fill             = yes
stroke_color     = white
stroke_thickness = 1p
show_label       = yes
label_font       = default
label_radius     = 1.030r
label_size       = 17p
label_parallel   = yes
show_bands       = no
</ideogram>

show_ticks       = yes
show_tick_labels = yes
<ticks>
radius           = 1r
color            = grey
thickness        = 1p
multiplier       = 1e-6
format           = %d
<tick>
spacing          = 10u
size             = 7p
show_label       = yes
label_size       = 10p
label_offset     = 3p
</tick>
</ticks>

<plots>
<plot>
type=histogram
file=data/gene_density.txt
r0=0.9348r
r1=0.9742r
min=0
max={gene_max:.6f}
fill_color=gene_track
color=gene_track
thickness=0
</plot>
<plot>
type=histogram
file=data/te_coverage.txt
r0=0.8660r
r1=0.9053r
min=0
max=1
fill_color=te_track
color=te_track
thickness=0
</plot>
<plot>
type=histogram
file=data/gypsy_coverage.txt
r0=0.7971r
r1=0.8364r
min=0
max={gypsy_max:.6f}
fill_color=gypsy_track
color=gypsy_track
thickness=0
</plot>
<plot>
type=histogram
file=data/copia_coverage.txt
r0=0.7282r
r1=0.7676r
min=0
max={copia_max:.6f}
fill_color=copia_track
color=copia_track
thickness=0
</plot>
<plot>
type=line
file=data/gc_fraction.txt
r0=0.6593r
r1=0.6987r
min={gc_min:.6f}
max={gc_max:.6f}
color=gc_track
thickness=2p
</plot>
</plots>

<links>
<link>
file=data/links_background.txt
radius=0.6298r
bezier_radius=0.0689r
ribbon=yes
thickness=1p
crest=0.45
</link>
<link>
file=data/links_highlight.txt
radius=0.6199r
bezier_radius=0.0590r
ribbon=yes
thickness=1p
crest=0.45
</link>
</links>
"""
    (conf_dir / "circos.conf").write_text(config, encoding="ascii")


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    out_dir = args.out_dir.resolve()
    renamed_dir = renamed_package_path(root, args.hap_rename_dir)
    if args.stage_jcvi_inputs:
        stage_jcvi_inputs(root, out_dir, renamed_dir)
        return

    jcvi_dir = (args.jcvi_dir or out_dir / "work" / "jcvi").resolve()
    package_mode = renamed_package_available(root, renamed_dir)
    if package_mode:
        print(f"[INFO] four-haplotype inputs: {renamed_dir}")
    te_seqid_aliases = renamed_seqid_aliases(root, renamed_dir) if package_mode else {}
    provenance_source_ids = {display: source for source, display in te_seqid_aliases.items()}
    source_to_display, display_to_source, quartet_meta = load_quartet_assignment(
        root, renamed_dir
    )
    data_dir = out_dir / "data"
    conf_dir = out_dir / "conf"
    figures_dir = out_dir / "figures"
    for directory in (data_dir, conf_dir, figures_dir):
        directory.mkdir(parents=True, exist_ok=True)

    all_sizes: dict[str, int] = {}
    all_windows: dict[str, list[tuple[int, int, float]]] = {}
    all_gene_counts: dict[str, list[int]] = {}
    all_te: dict[str, list[int]] = {}
    all_gypsy: dict[str, list[int]] = {}
    all_copia: dict[str, list[int]] = {}
    input_rows: list[tuple[str, str, str, str]] = []

    for hap in HAPS:
        genome = genome_path(root, hap, renamed_dir)
        gene_gff = gene_gff_path(root, hap, renamed_dir)
        te_gff = te_gff_path(root, hap)
        print(f"[INFO] hap{hap} genome: {genome}")
        print(f"[INFO] hap{hap} genes : {gene_gff}")
        print(f"[INFO] hap{hap} TE    : {te_gff}")
        sizes, windows = fasta_windows(genome, args.window)
        validate_chromosomes(sizes, hap, renamed=package_mode)
        gene_counts = parse_gene_counts(gene_gff, sizes, args.window)
        te_intervals = parse_te_intervals(te_gff, sizes, te_seqid_aliases)
        all_sizes.update(sizes)
        all_windows.update(windows)
        all_gene_counts.update(gene_counts)
        all_te.update(window_coverage(te_intervals["all"], sizes, args.window))
        all_gypsy.update(window_coverage(te_intervals["gypsy"], sizes, args.window))
        all_copia.update(window_coverage(te_intervals["copia"], sizes, args.window))
        input_rows.append((hap, str(genome), str(gene_gff), str(te_gff)))

    chromosome_order, scale_factors = chromosome_layout(
        all_sizes, source_to_display, display_to_source
    )
    with (data_dir / "karyotype.txt").open("w", encoding="ascii") as handle:
        for chrom in chromosome_order:
            match = DISPLAY_RE.fullmatch(chrom)
            assert match is not None
            quartet, hap = match.groups()
            source = display_to_source[chrom]
            handle.write(f"chr - {chrom} {chrom} 0 {all_sizes[source]} {HAP_COLORS[hap]}\n")

    with (data_dir / "display_orientation.tsv").open("w", encoding="ascii") as handle:
        handle.write("display_chromosome\torientation\n")
        for chrom in chromosome_order:
            orientation = "reverse" if chrom in REVERSED_DISPLAY_CHROMS else "forward"
            handle.write(f"{chrom}\t{orientation}\n")

    with (data_dir / "chromosome_rename_map.tsv").open("w", encoding="ascii") as handle:
        handle.write(
            "group\tsource_group\tquartet_id\tconfidence\tdisplay_chromosome\t"
            "original_chromosome\thaplotype\tmapping_note\n"
        )
        for number in sorted(quartet_meta):
            for hap in HAPS:
                display = f"{QUARTET_LABELS[number]}_{hap}"
                source = display_to_source[display]
                provenance_source = provenance_source_ids.get(source, source)
                source_match = CHR_RE.fullmatch(provenance_source)
                if source_match is None:
                    raise ValueError(f"Invalid source chromosome alias: {provenance_source}")
                source_number = int(source_match.group(1))
                source_group = QUARTET_LABELS[source_number]
                quartet_id, confidence = quartet_meta[source_number]
                mapping_note = (
                    "manual_A_C01_C06_swap"
                    if source_group != QUARTET_LABELS[number]
                    else "quartet_default"
                )
                handle.write(
                    f"{QUARTET_LABELS[number]}\t{source_group}\t{quartet_id}\t{confidence}\t"
                    f"{display}\t{provenance_source}\t{hap}\t{mapping_note}\n"
                )

    with (data_dir / "copy1_pair_map.tsv").open("w", encoding="ascii") as handle:
        handle.write("group\tcopy1_A1_slot\tcopy1_B1_slot\tcopy2_A2_slot\tcopy2_B2_slot\tstatus\n")
        for group, a1_slot, b1_slot, a2_slot, b2_slot, status in COPY1_PAIR_ROWS:
            handle.write(
                f"{group}\t{a1_slot}\t{b1_slot}\t{a2_slot}\t{b2_slot}\t{status}\n"
            )

    with (data_dir / "layout_scaling.tsv").open("w", encoding="ascii") as handle:
        handle.write("haplotype\ttotal_bp\tangular_scale\tquadrant\n")
        quadrants = {"A1": "upper-left", "B1": "upper-right", "A2": "lower-left", "B2": "lower-right"}
        for hap in HAPS:
            total = sum(
                all_sizes[display_to_source[f"C{i:02d}_{hap}"]]
                for i in range(1, 9)
            )
            handle.write(f"hap{hap}\t{total}\t{scale_factors[hap]:.8f}\t{quadrants[hap]}\n")

    gene_values = [float(value) for counts in all_gene_counts.values() for value in counts]
    gene_max = max(1.0, percentile(gene_values, 0.99))
    gc_values = [gc for windows in all_windows.values() for _, _, gc in windows if gc > 0]
    gc_low = percentile(gc_values, 0.01)
    gc_high = percentile(gc_values, 0.99)
    gc_padding = max(0.005, (gc_high - gc_low) * 0.08)
    gc_min = max(0.0, gc_low - gc_padding)
    gc_max = min(1.0, gc_high + gc_padding)
    gypsy_values = [
        all_gypsy[source][index] / (end - start)
        for source in all_windows
        for index, (start, end, _gc) in enumerate(all_windows[source])
    ]
    copia_values = [
        all_copia[source][index] / (end - start)
        for source in all_windows
        for index, (start, end, _gc) in enumerate(all_windows[source])
    ]
    gypsy_max = coverage_plot_max(gypsy_values)
    copia_max = coverage_plot_max(copia_values)

    with (data_dir / "window_tracks.tsv").open("w", encoding="ascii") as raw, \
         (data_dir / "gene_density.txt").open("w", encoding="ascii") as genes, \
         (data_dir / "te_coverage.txt").open("w", encoding="ascii") as tes, \
         (data_dir / "gypsy_coverage.txt").open("w", encoding="ascii") as gypsy_out, \
         (data_dir / "copia_coverage.txt").open("w", encoding="ascii") as copia_out, \
         (data_dir / "gc_fraction.txt").open("w", encoding="ascii") as gcs:
        raw.write("chrom\tstart\tend\tgene_count\tte_bp\tte_fraction\tgypsy_bp\tgypsy_fraction\tcopia_bp\tcopia_fraction\tgc_fraction\n")
        for chrom in chromosome_order:
            source = display_to_source[chrom]
            for index, (start, end, gc) in enumerate(all_windows[source]):
                width = end - start
                gene_count = all_gene_counts[source][index]
                te_fraction = all_te[source][index] / width
                gypsy_fraction = all_gypsy[source][index] / width
                copia_fraction = all_copia[source][index] / width
                raw.write(
                    f"{chrom}\t{start}\t{end}\t{gene_count}\t{all_te[source][index]}\t{te_fraction:.6f}\t"
                    f"{all_gypsy[source][index]}\t{gypsy_fraction:.6f}\t{all_copia[source][index]}\t"
                    f"{copia_fraction:.6f}\t{gc:.6f}\n"
                )
                genes.write(f"{chrom} {start} {end} {min(gene_count, gene_max):.6f}\n")
                tes.write(f"{chrom} {start} {end} {min(te_fraction, 1.0):.6f}\n")
                gypsy_out.write(f"{chrom} {start} {end} {min(gypsy_fraction, gypsy_max):.6f}\n")
                copia_out.write(f"{chrom} {start} {end} {min(copia_fraction, copia_max):.6f}\n")
                gcs.write(f"{chrom} {start} {end} {gc:.6f}\n")

    staged_beds = {hap: load_bed(out_dir / "work" / "jcvi_inputs" / f"hap{hap}.bed") for hap in HAPS}
    block_summary: list[tuple[str, str, int, int]] = []
    link_audit_rows: list[tuple[str, str, str, int, int, str, int, int, str, str, str]] = []
    endpoint_counts: dict[tuple[str, str, str], int] = defaultdict(int)
    endpoint_bp: dict[tuple[str, str, str], int] = defaultdict(int)
    background_path = data_dir / "links_background.txt"
    highlight_path = data_dir / "links_highlight.txt"
    with background_path.open("w", encoding="ascii") as background, highlight_path.open("w", encoding="ascii") as highlight:
        for left, right, layer in LINK_SPECS:
            simple = require_file(jcvi_dir / f"{left}_{right}" / "blocks.simple", f"{left}-{right} JCVI simple blocks")
            blocks = parse_jcvi_blocks(simple, staged_beds[left], staged_beds[right])
            output = background if layer == "background" else highlight
            cross_chrom = 0
            for left_chr, left_start, left_end, right_chr, right_start, right_end, score, orientation in blocks:
                left_display = source_to_display[left_chr]
                right_display = source_to_display[right_chr]
                if left_display.split("_", 1)[0] != right_display.split("_", 1)[0]:
                    cross_chrom += 1
                if layer == "background":
                    link_color = "within_gray"
                else:
                    display_match = DISPLAY_RE.fullmatch(left_display)
                    if display_match is None:
                        raise ValueError(f"Invalid display chromosome alias in JCVI link: {left_display}")
                    pair_label = f"C{display_match.group(1)}"
                    if left == "A1":
                        link_color = PAIR_COLORS_A1B1[pair_label]
                    else:
                        link_color = PAIR_COLORS_A2B2[pair_label]
                z = 5 if layer == "background" else 100
                if orientation not in {"+", "-"}:
                    raise ValueError(f"Unexpected JCVI block orientation: {orientation!r}")
                displayed_orientation = link_display_orientation(
                    orientation, left_display, right_display
                )
                if layer == "highlight":
                    # For the colored A-B layer, use JCVI's native orientation
                    # for twist detection so the explicit display reversals
                    # are not rendered twice.  Native '-' blocks remain
                    # visible and are retained in link_block_audit.tsv.
                    twist = ",twist=yes" if orientation == "-" else ",twist=no"
                else:
                    twist = (
                        ",twist=yes" if displayed_orientation == "-" else ",twist=no"
                    )
                alpha = BACKGROUND_ALPHA if layer == "background" else HIGHLIGHT_ALPHA
                color_name = f"{link_color}_{alpha}"
                output.write(
                    f"{left_display} {left_start} {left_end} {right_display} {right_start} {right_end} "
                    f"color={color_name},z={z}{twist}\n"
                )
                comparison = f"{left}-{right}"
                link_audit_rows.append(
                    (
                        comparison,
                        layer,
                        left_display,
                        left_start,
                        left_end,
                        right_display,
                        right_start,
                        right_end,
                        score,
                        orientation,
                        displayed_orientation,
                    )
                )
                for display, start, end in (
                    (left_display, left_start, left_end),
                    (right_display, right_start, right_end),
                ):
                    key = (comparison, layer, display)
                    endpoint_counts[key] += 1
                    endpoint_bp[key] += end - start
            block_summary.append((f"{left}-{right}", layer, len(blocks), cross_chrom))
            print(f"[INFO] {left}-{right}: {len(blocks)} JCVI blocks")

    with (data_dir / "jcvi_block_summary.tsv").open("w", encoding="ascii") as handle:
        handle.write("pair\tlayer\tblocks\tcross_chromosome_blocks\n")
        for row in block_summary:
            handle.write("\t".join(map(str, row)) + "\n")

    # Keep direct evidence for chromosome-level interpretation.  In
    # particular, this distinguishes a real cross-quartet rearrangement from
    # a chromosome that lost all blocks during the JCVI filtering step.
    with (data_dir / "link_block_audit.tsv").open("w", encoding="ascii") as handle:
        handle.write(
            "comparison\tlayer\tleft_chromosome\tleft_start\tleft_end\t"
            "right_chromosome\tright_start\tright_end\tscore\t"
            "native_orientation\tdisplayed_orientation\n"
        )
        for row in link_audit_rows:
            handle.write("\t".join(map(str, row)) + "\n")

    with (data_dir / "link_coverage_by_chromosome.tsv").open("w", encoding="ascii") as handle:
        handle.write("comparison\tlayer\thaplotype\tdisplay_chromosome\tblocks\tendpoint_bp\n")
        for left, right, layer in LINK_SPECS:
            comparison = f"{left}-{right}"
            for hap in (left, right):
                for display in chromosome_order:
                    if display_haplotype(display) != hap:
                        continue
                    key = (comparison, layer, display)
                    handle.write(
                        f"{comparison}\t{layer}\t{hap}\t{display}\t"
                        f"{endpoint_counts[key]}\t{endpoint_bp[key]}\n"
                    )

    with (data_dir / "input_manifest.tsv").open("w", encoding="utf-8") as handle:
        handle.write("haplotype\tgenome\tgene_gff3\tte_gff3\n")
        for row in input_rows:
            handle.write("\t".join(row) + "\n")

    with (data_dir / "color_legend.tsv").open("w", encoding="ascii") as handle:
        handle.write("element\tcolor_definition\tdescription\n")
        handle.write("hapa1\t55,95,145\thapA1 ideograms, upper-left\n")
        handle.write("hapb1\t186,85,40\thapB1 ideograms, upper-right\n")
        handle.write("hapa2\t120,168,205\thapA2 ideograms, lower-left\n")
        handle.write("hapb2\t222,151,55\thapB2 ideograms, lower-right\n")
        handle.write("gene_track\t85,90,95\tgene count per window\n")
        handle.write("te_track\t112,135,52\tnon-overlapping total TE coverage\n")
        handle.write("gypsy_track\t196,91,82\tnon-overlapping Gypsy coverage\n")
        handle.write("copia_track\t76,120,168\tnon-overlapping Copia coverage\n")
        handle.write("gc_track\t42,88,108\tGC fraction\n")
        handle.write("within_gray\t175,175,175\tA1-A2 and B1-B2 JCVI blocks\n")
        handle.write(f"gypsy_ymax\t{gypsy_max:.6f}\tGypsy histogram display maximum\n")
        handle.write(f"copia_ymax\t{copia_max:.6f}\tCopia histogram display maximum\n")
        handle.write(
            f"highlight_alpha\t{HIGHLIGHT_ALPHA}\t"
            f"tuned highlight transparency ({HIGHLIGHT_TRANSPARENCY:.2%})\n"
        )
        handle.write(
            f"background_alpha\t{BACKGROUND_ALPHA}\t"
            f"background transparency ({BACKGROUND_ALPHA_INDEX / AUTO_ALPHA_STEPS:.2%})\n"
        )
        for label in PAIR_LABELS:
            handle.write(
                f"{PAIR_COLORS_A1B1[label]}_{HIGHLIGHT_ALPHA}\t"
                f"{PAIR_RGB_A1B1[label]}; {HIGHLIGHT_ALPHA}\t"
                f"{label} A1-B1 blocks\n"
            )
            handle.write(
                f"{PAIR_COLORS_A2B2[label]}_{HIGHLIGHT_ALPHA}\t"
                f"{PAIR_RGB_A2B2[label]}; {HIGHLIGHT_ALPHA}\t"
                f"{label} A2-B2 blocks\n"
            )

    with (data_dir / "track_layout.tsv").open("w", encoding="ascii") as handle:
        handle.write("panel_label\ttrack\tr0\tr1\tmin\tmax\tnote\n")
        handle.write(f"A\tgene_density\t0.9348r\t0.9742r\t0\t{gene_max:.6f}\tgene count per window\n")
        handle.write("B\ttotal_TE\t0.8660r\t0.9053r\t0\t1\tnon-overlapping total TE coverage\n")
        handle.write(f"C\tGypsy\t0.7971r\t0.8364r\t0\t{gypsy_max:.6f}\tcompressed display y-axis\n")
        handle.write(f"D\tCopia\t0.7282r\t0.7676r\t0\t{copia_max:.6f}\tcompressed display y-axis\n")
        handle.write(f"E\tGC\t0.6593r\t0.6987r\t{gc_min:.6f}\t{gc_max:.6f}\tGC fraction\n")

    write_circos_config(
        conf_dir,
        chromosome_order,
        scale_factors,
        gene_max,
        gc_min,
        gc_max,
        gypsy_max,
        copia_max,
    )
    print(f"[DONE] JCVI-based Circos inputs: {out_dir}")


if __name__ == "__main__":
    main()
