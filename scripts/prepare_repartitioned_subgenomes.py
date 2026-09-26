#!/usr/bin/env python3
"""Prepare the revised A/B chromosome packages for PZI, PPA and PSP.

The revised partition is deliberately keyed to PZI copy1 chromosome IDs.  It
does not reuse the former C01--C08 quartet labels and it does not infer a
one-to-one homeolog map.  A chromosome can therefore be linked to more than
one chromosome in the audit table, while each output haplotype still contains
exactly eight chromosomes.

The script creates genome-only inputs for a fresh EDTA -> Helixer/EviAnn ->
BUSCO run.  Existing annotation files are not copied into the new package;
this prevents an old annotation with old chromosome IDs from being silently
attached to the new partition.

PZI partitioning can start from the already resolved copy1/copy2
16-chromosome FASTA files, or from the complete 2n=32 FASTA.  The latter is
the coordinate-compatible mode for splitting a complete-assembly annotation;
the copy mode remains available for phase-specific genome analyses.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
from collections import Counter
from pathlib import Path
from typing import Iterable, TextIO


A_COPY1 = ("chr01", "chr02", "chr05", "chr06", "chr07", "chr10", "chr12", "chr14")
B_COPY1 = ("chr03", "chr04", "chr08", "chr09", "chr11", "chr13", "chr15", "chr16")
HOMOLOGY_EDGES = (
    ("chr01", "chr13"),
    ("chr02", "chr15"),
    ("chr02", "chr11"),
    ("chr05", "chr16"),
    ("chr06", "chr03"),
    ("chr12", "chr03"),
    ("chr07", "chr04"),
    ("chr10", "chr08"),
    ("chr14", "chr09"),
)


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def fasta_records(path: Path) -> Iterable[tuple[str, str, str]]:
    name = None
    header = ""
    chunks: list[str] = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, header, "".join(chunks)
                header = line[1:]
                name = header.split(None, 1)[0]
                chunks = []
            elif line.strip():
                if name is None:
                    raise ValueError(f"Sequence data precedes FASTA header in {path}")
                chunks.append(line.strip())
    if name is not None:
        yield name, header, "".join(chunks)


def load_fasta(path: Path) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    for name, header, seq in fasta_records(path):
        if name in records:
            raise ValueError(f"Duplicate FASTA sequence ID in {path}: {name}")
        records[name] = (header, seq)
    if not records:
        raise ValueError(f"FASTA is empty: {path}")
    return records


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"TSV has no header: {path}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_fasta(path: Path, rows: list[tuple[str, str, str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for new_id, source_id, header, sequence in rows:
            handle.write(f">{new_id} source_seqid={source_id} source_header={header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def norm_chr(value: str) -> str:
    value = value.strip()
    if re.fullmatch(r"\d+", value):
        return f"chr{int(value):02d}"
    match = re.fullmatch(r"chr(\d+)", value, flags=re.IGNORECASE)
    if match:
        return f"chr{int(match.group(1)):02d}"
    return value


def natural_key(value: str) -> tuple[object, ...]:
    parts: list[object] = []
    token = ""
    for char in value:
        if char.isdigit():
            token += char
        else:
            if token:
                parts.append(int(token))
                token = ""
            parts.append(char)
    if token:
        parts.append(int(token))
    return tuple(parts)


def parse_pzi_map(path: Path) -> tuple[dict[str, dict[str, str]], dict[str, str], dict[str, str]]:
    rows = read_tsv(path)
    by_copy1: dict[str, dict[str, str]] = {}
    by_copy2: dict[str, str] = {}
    display_to_copy1: dict[str, str] = {}

    if not rows:
        raise ValueError(f"PZI chromosome map is empty: {path}")

    # Current provenance output: one row per corrected copy1/copy2 pair.
    # Unlike the retired chromosome_rename_map.tsv, this table is independent
    # of the former C01--C08 quartet naming.
    provenance_required = {"copy1_chr", "copy2_chr", "original_copy1_chr", "original_copy2_chr"}
    if not provenance_required.issubset(rows[0]):
        raise ValueError(
            "This revised workflow requires the current provenance table with "
            f"columns {sorted(provenance_required)}; legacy chromosome_rename_map.tsv "
            f"tables are intentionally rejected: {path}"
        )

    # Rows are pair01..pair16.  The copy1/copy2 chromosome names and
    # original_copy*_chr columns are retained as provenance audit fields;
    # sequence extraction below uses the resolved copy1/copy2 chr16 FASTAs.
    rows = sorted(rows, key=lambda row: natural_key(row.get("pair_id", "")))
    for index, row in enumerate(rows, 1):
        normalized = f"chr{index:02d}"
        original_copy1 = norm_chr(row["original_copy1_chr"])
        original_copy2 = norm_chr(row["original_copy2_chr"])
        if normalized in by_copy1:
            raise ValueError(f"Duplicate normalized copy chromosome in {path}: {normalized}")
        by_copy1[normalized] = {
            "original_chr32": original_copy1,
            "copy1": normalized,
            "copy2": normalized,
            "copy1_source": norm_chr(row["copy1_chr"]),
            "copy2_source": norm_chr(row["copy2_chr"]),
            **row,
        }
        by_copy2[normalized] = original_copy2
    if set(by_copy1) != {f"chr{i:02d}" for i in range(1, 17)}:
        raise ValueError("PZI provenance map must contain copy1 chr01--chr16 exactly once")
    if set(by_copy2) != set(by_copy1):
        raise ValueError("PZI provenance map must contain one copy2 partner for every copy1 chromosome")
    return by_copy1, display_to_copy1, by_copy2


def group_for_copy1(copy1: str) -> tuple[str, int]:
    copy1 = norm_chr(copy1)
    if copy1 in A_COPY1:
        return "A", A_COPY1.index(copy1) + 1
    if copy1 in B_COPY1:
        return "B", B_COPY1.index(copy1) + 1
    raise ValueError(f"copy1 chromosome is not in the revised A/B partition: {copy1}")


PZI_TARGET_RE = re.compile(
    r"^C0?([1-8])_?([AB])(?:_?([AB])?([12]))?$", re.IGNORECASE
)


def copy1_for_pzi_target(target: str) -> str | None:
    """Resolve C1A_A1/C01_A1/C1A targets to the revised copy1 chromosome."""
    match = PZI_TARGET_RE.fullmatch(target.strip())
    if match is None:
        return None
    rank, subgenome, phase_subgenome, _phase = match.groups()
    subgenome = subgenome.upper()
    if phase_subgenome and phase_subgenome.upper() != subgenome:
        raise ValueError(f"PZI target has inconsistent subgenome/phase: {target}")
    group = A_COPY1 if subgenome == "A" else B_COPY1
    return group[int(rank) - 1]


def normalize_copy_records(path: Path) -> dict[str, tuple[str, str]]:
    """Load a copy1/copy2 FASTA and require normalized chr01--chr16 IDs."""
    records = load_fasta(path)
    normalized: dict[str, tuple[str, str]] = {}
    for name, value in records.items():
        chrom = norm_chr(name)
        if chrom in normalized:
            raise ValueError(f"Duplicate normalized chromosome ID in {path}: {chrom}")
        normalized[chrom] = value
    expected = {f"chr{i:02d}" for i in range(1, 17)}
    if set(normalized) != expected:
        missing = sorted(expected - set(normalized))
        extra = sorted(set(normalized) - expected)
        raise ValueError(
            f"{path} must contain exactly chr01--chr16; missing={missing}, extra={extra}"
        )
    return normalized


def prepare_pzi(copy1_fasta: Path, copy2_fasta: Path, pzi_map: Path, out_root: Path) -> list[dict[str, object]]:
    by_copy1, _, _ = parse_pzi_map(pzi_map)
    copy1_records = normalize_copy_records(copy1_fasta)
    copy2_records = normalize_copy_records(copy2_fasta)
    rows: list[dict[str, object]] = []
    for phase, records, source_fasta, copy_name in (
        ("A1", copy1_records, copy1_fasta, "copy1"),
        ("A2", copy2_records, copy2_fasta, "copy2"),
        ("B1", copy1_records, copy1_fasta, "copy1"),
        ("B2", copy2_records, copy2_fasta, "copy2"),
    ):
        subgenome = phase[0]
        phase_number = phase[1]
        output_rows: list[tuple[str, str, str, str]] = []
        order = A_COPY1 if subgenome == "A" else B_COPY1
        for copy1 in order:
            mapping = by_copy1[copy1]
            source_id = copy1
            header, sequence = records[source_id]
            new_id = f"C{order.index(copy1) + 1}{subgenome}"
            output_rows.append((new_id, source_id, header, sequence))
            original_chr32 = (
                norm_chr(mapping["original_copy1_chr"])
                if copy_name == "copy1"
                else norm_chr(mapping["original_copy2_chr"])
            )
            rows.append(
                {
                    "species": "pzi",
                    "haplotype": f"hap{phase}",
                    "new_seqid": new_id,
                    "copy1_chr": copy1,
                    "copy2_chr": copy1,
                    "original_chr32": original_chr32,
                    "copy1_source_chr32": norm_chr(mapping["original_copy1_chr"]),
                    "copy2_source_chr32": norm_chr(mapping["original_copy2_chr"]),
                    "pair_id": mapping.get("pair_id", ""),
                    "subgenome": subgenome,
                    "phase": phase_number,
                    "source_fasta": str(source_fasta),
                }
            )
        write_fasta(out_root / "inputs" / "pzi" / "genome" / f"pzi_hap{phase}.fa", output_rows)
    return rows


def prepare_pzi_from_full32(full_fasta: Path, pzi_map: Path, out_root: Path) -> list[dict[str, object]]:
    """Build partition genomes directly from the complete 2n=32 assembly.

    This mode is required when annotations were generated on the complete
    assembly: selecting the mapped ``original_chr32`` records preserves the
    exact coordinate system used by the GFF3.  The copy1/copy2 mode above is
    retained for phase-specific genome analyses, but its sequences are not
    necessarily coordinate-identical to the complete assembly.
    """
    by_copy1, _, _ = parse_pzi_map(pzi_map)
    records = load_fasta(full_fasta)
    rows: list[dict[str, object]] = []
    for phase in ("A1", "A2", "B1", "B2"):
        subgenome = phase[0]
        phase_number = phase[1]
        order = A_COPY1 if subgenome == "A" else B_COPY1
        output_rows: list[tuple[str, str, str, str]] = []
        for rank, copy1 in enumerate(order, 1):
            mapping = by_copy1[copy1]
            source_id = norm_chr(
                mapping["original_copy1_chr"] if phase_number == "1" else mapping["original_copy2_chr"]
            )
            if source_id not in records:
                raise ValueError(
                    f"Complete pzi32 FASTA lacks mapped chromosome {source_id} "
                    f"for {phase}/{copy1}"
                )
            header, sequence = records[source_id]
            new_id = f"C{rank}{subgenome}"
            output_rows.append((new_id, source_id, header, sequence))
            rows.append(
                {
                    "species": "pzi",
                    "haplotype": f"hap{phase}",
                    "new_seqid": new_id,
                    "copy1_chr": copy1,
                    "copy2_chr": copy1,
                    "original_chr32": source_id,
                    "copy1_source_chr32": norm_chr(mapping["original_copy1_chr"]),
                    "copy2_source_chr32": norm_chr(mapping["original_copy2_chr"]),
                    "pair_id": mapping.get("pair_id", ""),
                    "subgenome": subgenome,
                    "phase": phase_number,
                    "source_fasta": str(full_fasta),
                }
            )
        write_fasta(out_root / "inputs" / "pzi" / "genome" / f"pzi_hap{phase}.fa", output_rows)
    return rows


def prepare_ppa(genome: Path, assignment: Path, out_root: Path) -> list[dict[str, object]]:
    records = load_fasta(genome)
    assignments = read_tsv(assignment)
    required = {"ppa_source_seqid", "pzi_copy1_chr"}
    if not assignments or not required.issubset(assignments[0]):
        raise ValueError(f"PPA assignment lacks {sorted(required)}: {assignment}")
    if len(assignments) != 16:
        raise ValueError(f"PPA assignment must contain 16 primary chromosomes: {assignment}")
    by_hap: dict[str, list[tuple[str, str, str, str]]] = {"A": [], "B": []}
    rows: list[dict[str, object]] = []
    for row in assignments:
        source = row["ppa_source_seqid"]
        if source not in records:
            raise ValueError(f"PPA FASTA lacks assigned chromosome {source}")
        subgenome, rank = group_for_copy1(row["pzi_copy1_chr"])
        new_id = f"C{rank}{subgenome}"
        header, sequence = records[source]
        by_hap[subgenome].append((new_id, source, header, sequence))
        rows.append(
            {
                "species": "ppa",
                "haplotype": f"hap{subgenome}",
                "new_seqid": new_id,
                "source_seqid": source,
                "copy1_chr": norm_chr(row["pzi_copy1_chr"]),
                "subgenome": subgenome,
                "source_assignment": str(assignment),
            }
        )
    for subgenome in ("A", "B"):
        if len(by_hap[subgenome]) != 8:
            raise ValueError(f"PPA {subgenome} group has {len(by_hap[subgenome])} chromosomes, expected 8")
        by_hap[subgenome].sort(key=lambda item: item[0])
        write_fasta(out_root / "inputs" / "ppa" / "genome" / f"ppa_hap{subgenome}.fa", by_hap[subgenome])
    return rows


def prepare_psp(genome: Path, assignment: Path, pzi_display_map: dict[str, str], out_root: Path) -> list[dict[str, object]]:
    records = load_fasta(genome)
    assignments = read_tsv(assignment)
    required = {"new_seqid", "pzi_target", "phase_suffix"}
    if not assignments or not required.issubset(assignments[0]):
        raise ValueError(f"PSP assignment lacks {sorted(required)}: {assignment}")
    by_hap: dict[str, list[tuple[str, str, str, str]]] = {"A1": [], "A2": [], "B1": [], "B2": []}
    rows: list[dict[str, object]] = []
    seen_new: set[tuple[str, str]] = set()
    for row in assignments:
        old_id = row["new_seqid"]
        if old_id not in records:
            raise ValueError(f"PSP renamed FASTA lacks assigned sequence {old_id}")
        target = row["pzi_target"]
        copy1 = pzi_display_map.get(target) or copy1_for_pzi_target(target)
        if not copy1:
            # Some maps use a chrXX target directly rather than a CxA/CxB label.
            copy1 = norm_chr(target)
        subgenome, rank = group_for_copy1(copy1)
        phase = row["phase_suffix"]
        if phase not in {"1", "2"}:
            raise ValueError(f"Invalid PSP phase_suffix for {old_id}: {phase}")
        hap = f"{subgenome}{phase}"
        new_id = f"C{rank}{subgenome}"
        if (hap, new_id) in seen_new:
            raise ValueError(f"Duplicate revised PSP chromosome assignment: {hap}/{new_id}")
        seen_new.add((hap, new_id))
        header, sequence = records[old_id]
        by_hap[hap].append((new_id, old_id, header, sequence))
        rows.append(
            {
                "species": "psp",
                "haplotype": f"hap{hap}",
                "new_seqid": new_id,
                "old_seqid": old_id,
                "pzi_target": target,
                "copy1_chr": copy1,
                "subgenome": subgenome,
                "phase": phase,
                "source_assignment": str(assignment),
            }
        )
    for hap in by_hap:
        if len(by_hap[hap]) != 8:
            raise ValueError(f"PSP {hap} group has {len(by_hap[hap])} chromosomes, expected 8")
        by_hap[hap].sort(key=lambda item: item[0])
        write_fasta(out_root / "inputs" / "psp" / "genome" / f"psp_hap{hap}.fa", by_hap[hap])
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pzi-copy1", type=Path, help="resolved 16-chromosome copy1 FASTA")
    parser.add_argument("--pzi-copy2", type=Path, help="resolved 16-chromosome copy2 FASTA")
    parser.add_argument(
        "--pzi-full32",
        type=Path,
        help="complete pzi32 chromosome FASTA; use this when annotations are on the complete assembly",
    )
    parser.add_argument("--pzi-map", type=Path, required=True, help="corrected pair provenance table")
    parser.add_argument("--out-root", type=Path, required=True)
    parser.add_argument("--ppa-genome", type=Path)
    parser.add_argument("--ppa-assignment", type=Path)
    parser.add_argument("--psp-genome", type=Path)
    parser.add_argument("--psp-assignment", type=Path)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.out_root.exists() and any(args.out_root.iterdir()) and not args.force:
        raise SystemExit(f"ERROR: output root is non-empty; use --force or a new directory: {args.out_root}")
    if args.pzi_full32 is None and (args.pzi_copy1 is None or args.pzi_copy2 is None):
        raise SystemExit("ERROR: provide --pzi-full32, or provide both --pzi-copy1 and --pzi-copy2")
    if args.pzi_full32 is not None and (args.pzi_copy1 is not None or args.pzi_copy2 is not None):
        raise SystemExit("ERROR: --pzi-full32 cannot be combined with --pzi-copy1/--pzi-copy2")
    input_files = [("PZI provenance map", args.pzi_map)]
    if args.pzi_full32 is not None:
        input_files.append(("PZI complete 2n=32 FASTA", args.pzi_full32))
    else:
        input_files.extend(
            [("PZI copy1 FASTA", args.pzi_copy1), ("PZI copy2 FASTA", args.pzi_copy2)]
        )
    for label, path in input_files:
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"ERROR: missing {label}: {path}")
    args.out_root.mkdir(parents=True, exist_ok=True)
    if args.pzi_full32 is not None:
        pzi_rows = prepare_pzi_from_full32(args.pzi_full32, args.pzi_map, args.out_root)
    else:
        pzi_rows = prepare_pzi(args.pzi_copy1, args.pzi_copy2, args.pzi_map, args.out_root)
    _, display_map, _ = parse_pzi_map(args.pzi_map)
    ppa_rows: list[dict[str, object]] = []
    psp_rows: list[dict[str, object]] = []
    if (args.ppa_genome is None) != (args.ppa_assignment is None):
        raise SystemExit("ERROR: provide both --ppa-genome and --ppa-assignment, or neither")
    if args.ppa_genome and args.ppa_assignment:
        ppa_rows = prepare_ppa(args.ppa_genome, args.ppa_assignment, args.out_root)
    if (args.psp_genome is None) != (args.psp_assignment is None):
        raise SystemExit("ERROR: provide both --psp-genome and --psp-assignment, or neither")
    if args.psp_genome and args.psp_assignment:
        psp_rows = prepare_psp(args.psp_genome, args.psp_assignment, display_map, args.out_root)

    map_dir = args.out_root / "maps"
    map_dir.mkdir(parents=True, exist_ok=True)
    pzi_fields = [
        "species", "haplotype", "new_seqid", "copy1_chr", "copy2_chr",
        "original_chr32", "copy1_source_chr32", "copy2_source_chr32",
        "pair_id", "subgenome", "phase", "source_fasta",
    ]
    write_tsv(map_dir / "pzi_repartition_map.tsv", pzi_rows, pzi_fields)
    if ppa_rows:
        write_tsv(map_dir / "ppa_repartition_map.tsv", ppa_rows, ["species", "haplotype", "new_seqid", "source_seqid", "copy1_chr", "subgenome", "source_assignment"])
    if psp_rows:
        write_tsv(map_dir / "psp_repartition_map.tsv", psp_rows, ["species", "haplotype", "new_seqid", "old_seqid", "pzi_target", "copy1_chr", "subgenome", "phase", "source_assignment"])
    write_tsv(
        map_dir / "revised_homology_edges.tsv",
        [
            {
                "a_copy1_chr": a,
                "a_chromosome": f"C{group_for_copy1(a)[1]}A",
                "b_copy1_chr": b,
                "b_chromosome": f"C{group_for_copy1(b)[1]}B",
            }
            for a, b in HOMOLOGY_EDGES
        ],
        ["a_copy1_chr", "a_chromosome", "b_copy1_chr", "b_chromosome"],
    )
    validation = []
    for species, rows in (("pzi", pzi_rows), ("ppa", ppa_rows), ("psp", psp_rows)):
        if not rows:
            continue
        key = "haplotype"
        counts = Counter(row[key] for row in rows)
        expected = 4 if species in {"pzi", "psp"} else 2
        validation.append({"species": species, "check": "haplotype_count", "observed": len(counts), "expected": expected, "status": "pass" if len(counts) == expected else "fail"})
        for hap, count in sorted(counts.items()):
            validation.append({"species": species, "check": f"{hap}_chromosome_count", "observed": count, "expected": 8, "status": "pass" if count == 8 else "fail"})
    write_tsv(args.out_root / "review" / "assignment_validation.tsv", validation, ["species", "check", "observed", "expected", "status"])
    manifest = []
    for path in sorted(args.out_root.rglob("*.fa")):
        manifest.append({"file": str(path.relative_to(args.out_root)), "bytes": path.stat().st_size, "md5": md5(path)})
    write_tsv(args.out_root / "MANIFEST.tsv", manifest, ["file", "bytes", "md5"])
    readme = args.out_root / "README.txt"
    pzi_source_note = (
        "PZI sequences are read directly from the complete pzi32 FASTA when "
        "--pzi-full32 is used; this preserves complete-assembly annotation coordinates.\n"
        if args.pzi_full32 is not None
        else "PZI sequences are read directly from the resolved copy1/copy2 chr16 FASTAs; "
        "the 2n=32 FASTA is not used for repartitioning.\n"
    )
    readme.write_text(
        "Revised A/B chromosome packages\n\n"
        "PZI copy1 partition: A=chr01,chr02,chr05,chr06,chr07,chr10,chr12,chr14; "
        "B=chr03,chr04,chr08,chr09,chr11,chr13,chr15,chr16.\n"
        + pzi_source_note
        + "Output sequence IDs are C1A-C8A or C1B-C8B; phase is carried by the file label.\n"
        + "The supplied homology table is retained as an edge list and is not forced into a bijection.\n"
        + "These genome-only inputs are intended for a fresh independent EDTA and gene-annotation run.\n",
        encoding="utf-8",
    )
    failures = [row for row in validation if row["status"] == "fail"]
    if failures:
        raise SystemExit(f"ERROR: validation failed: {failures}")
    print(f"[OK] Revised packages written to {args.out_root}")


if __name__ == "__main__":
    main()
