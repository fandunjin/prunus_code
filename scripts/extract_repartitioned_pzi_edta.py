#!/usr/bin/env python3
"""Extract the frozen 2n=32 PZI EDTA annotation for the revised haplotypes.

The revised PZI haplotypes contain unmodified chromosomes selected from the
frozen 2n=32 assembly.  This program verifies that relationship at the
sequence level before rewriting chromosome identifiers in the EDTA GFF3 and
extracting the matching masked sequence records.  All four outputs therefore
use the same TE discovery and classification library.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import shutil
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator


HAPLOTYPES = ("hapA1", "hapA2", "hapB1", "hapB2")


def die(message: str) -> None:
    raise RuntimeError(message)


def require_file(path: Path, label: str) -> Path:
    if not path.is_file() or path.stat().st_size == 0:
        die(f"Missing or empty {label}: {path}")
    return path


def canonical_seqid(value: str) -> str:
    """Normalize numeric chromosome aliases while preserving non-numeric IDs."""
    value = value.strip().split(None, 1)[0]
    match = re.fullmatch(r"(?:chr)?0*(\d+)", value, flags=re.IGNORECASE)
    if match:
        return f"chr{int(match.group(1)):02d}"
    return value


def fasta_records(path: Path, uppercase: bool = True) -> Iterator[tuple[str, str, str]]:
    header: str | None = None
    chunks: list[str] = []
    with path.open("r", encoding="ascii", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(chunks)
                    yield header.split(None, 1)[0], header, sequence.upper() if uppercase else sequence
                header = line[1:]
                chunks = []
            else:
                if header is None:
                    die(f"Sequence data before FASTA header in {path}")
                chunks.append(line)
    if header is not None:
        sequence = "".join(chunks)
        yield header.split(None, 1)[0], header, sequence.upper() if uppercase else sequence


def md5_text(sequence: str) -> str:
    return hashlib.md5(sequence.encode("ascii")).hexdigest()


def write_fasta_record(handle, header: str, sequence: str) -> None:
    handle.write(f">{header}\n")
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start : start + 80] + "\n")


def read_map(path: Path) -> dict[str, dict[str, str]]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"species", "haplotype", "new_seqid", "original_chr32"}
    if not rows or not required.issubset(rows[0]):
        die(f"PZI repartition map lacks required columns {sorted(required)}: {path}")

    by_haplotype: dict[str, dict[str, str]] = {hap: {} for hap in HAPLOTYPES}
    assigned_sources: dict[str, tuple[str, str]] = {}
    for row in rows:
        if row.get("species") != "pzi" or row.get("haplotype") not in by_haplotype:
            continue
        hap = row["haplotype"]
        target = row["new_seqid"].strip()
        source = canonical_seqid(row["original_chr32"])
        if not target or not source:
            die(f"Invalid PZI assignment row: {row}")
        if target in by_haplotype[hap]:
            die(f"Duplicate target chromosome {hap}/{target} in {path}")
        if source in assigned_sources:
            other_hap, other_target = assigned_sources[source]
            die(
                f"Source chromosome {source} is assigned twice: "
                f"{other_hap}/{other_target} and {hap}/{target}"
            )
        by_haplotype[hap][target] = source
        assigned_sources[source] = (hap, target)

    for hap, mapping in by_haplotype.items():
        expected_prefix = "C"
        if len(mapping) != 8 or any(not target.startswith(expected_prefix) for target in mapping):
            die(f"Expected exactly eight revised PZI chromosomes for {hap}, found {len(mapping)}")
    if len(assigned_sources) != 32:
        die(f"Expected 32 uniquely assigned source chromosomes, found {len(assigned_sources)}")
    return by_haplotype


def output_paths(out_root: Path, haplotype: str) -> dict[str, Path]:
    prefix = f"pzi_{haplotype}.fa"
    base = out_root / f"{haplotype}_edta"
    return {
        "dir": base,
        "genome": base / prefix,
        "masked": base / f"{prefix}.mod.MAKER.masked",
        "gff": base / f"{prefix}.mod.EDTA.TEanno.gff3",
        "sum": base / f"{prefix}.mod.EDTA.TEanno.sum",
        "table": base / f"{prefix}.mod.EDTA.TEanno.extracted.tsv",
        "library": base / f"{prefix}.mod.EDTA.TElib.fa",
        "completed": base / "edta.completed.tsv",
        "chromosome_map": base / "chromosome_map.tsv",
    }


def replace_output(path: Path, force: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() or path.is_symlink():
        if not force:
            die(f"Refusing to replace existing output (use --force): {path}")
        path.unlink()


def link_or_copy(source: Path, target: Path, force: bool) -> None:
    replace_output(target, force)
    try:
        target.symlink_to(source.resolve())
    except OSError:
        shutil.copy2(source, target)
    if not target.is_file() or target.stat().st_size == 0:
        die(f"Failed to materialize {target}")


def parse_attributes(text: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for part in text.rstrip(";\n").split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            values[key] = value
    return values


def classify_repeat(attributes: str) -> tuple[str, str]:
    """Map EDTA GFF3 classification strings to reportable class/subclass keys."""
    values = parse_attributes(attributes)
    value = (
        values.get("Classification")
        or values.get("classification")
        or values.get("Class")
        or values.get("class")
        or ""
    ).strip()
    if not value:
        name = values.get("Name", "")
        if "#" in name:
            value = name.rsplit("#", 1)[1]
    value = value.replace(" ", "")
    tokens = [token for token in re.split(r"[/|]", value) if token]
    parent = tokens[0] if tokens else "unknown"
    child = tokens[1] if len(tokens) > 1 else "unknown"
    parent_lower = parent.lower()

    if parent_lower == "ltr":
        return "LTR", child
    if parent_lower == "line":
        return "LINE", child
    if parent_lower == "sine":
        return "SINE", child
    if parent_lower in {"dna", "tir", "mite"}:
        return "TIR", child if child != "unknown" else parent
    if parent_lower in {"rc", "helitron"}:
        return "nonTIR", "helitron" if parent_lower == "helitron" else child
    if parent_lower == "nonltr":
        return "nonLTR", child
    if parent_lower in {"rdna", "rrna"}:
        return "rDNA", child
    if parent_lower in {"low_complexity", "simple_repeat", "satellite", "unknown", "unclassified"}:
        return parent, child
    return parent, child


@dataclass
class UnionState:
    start: int
    end: int


class RepeatSummary:
    """Summarize sorted EDTA GFF3 intervals without double-counting overlap."""

    def __init__(self) -> None:
        self.counts: Counter[tuple[str, str]] = Counter()
        self.bp: Counter[tuple[str, str]] = Counter()
        self.total_count = 0
        self.total_bp = 0
        self._states: dict[tuple[str, str, str], UnionState] = {}
        self._total_states: dict[str, UnionState] = {}

    @staticmethod
    def _add(state: UnionState | None, start: int, end: int) -> tuple[UnionState, int]:
        if state is None:
            return UnionState(start, end), 0
        if start < state.start:
            die("EDTA GFF3 features are not coordinate-sorted; cannot calculate an exact union summary")
        if start <= state.end + 1:
            state.end = max(state.end, end)
            return state, 0
        return UnionState(start, end), state.end - state.start + 1

    def add(self, chromosome: str, start: int, end: int, category: str, subclass: str) -> None:
        key = (chromosome, category, subclass)
        state, completed = self._add(self._states.get(key), start, end)
        self._states[key] = state
        self.bp[(category, subclass)] += completed
        total_state, total_completed = self._add(self._total_states.get(chromosome), start, end)
        self._total_states[chromosome] = total_state
        self.total_bp += total_completed
        self.counts[(category, subclass)] += 1
        self.total_count += 1

    def finish(self) -> None:
        for (_, category, subclass), state in self._states.items():
            self.bp[(category, subclass)] += state.end - state.start + 1
        for state in self._total_states.values():
            self.total_bp += state.end - state.start + 1


def write_summary(paths: dict[str, Path], stats: RepeatSummary, genome_length: int) -> None:
    rows = sorted(stats.counts, key=lambda item: (item[0].lower(), item[1].lower()))
    with paths["table"].open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("class\tsubclass\tcount\tbp_union\tpercent_genome\n")
        for category, subclass in rows:
            bp = stats.bp[(category, subclass)]
            handle.write(
                f"{category}\t{subclass}\t{stats.counts[(category, subclass)]}\t{bp}\t"
                f"{100 * bp / genome_length:.4f}\n"
            )

    grouped: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    for category, subclass in rows:
        grouped[category].append((subclass, stats.counts[(category, subclass)], stats.bp[(category, subclass)]))
    ordered_categories = ["LINE", "SINE", "LTR", "TIR", "nonLTR", "nonTIR", "rDNA"]
    ordered_categories.extend(sorted(set(grouped) - set(ordered_categories), key=str.lower))
    interspersed_bp = sum(
        sum(bp for _, _, bp in grouped[category])
        for category in ("LINE", "SINE", "LTR", "TIR", "nonLTR", "nonTIR")
        if category in grouped
    )
    interspersed_count = sum(
        sum(count for _, count, _ in grouped[category])
        for category in ("LINE", "SINE", "LTR", "TIR", "nonLTR", "nonTIR")
        if category in grouped
    )
    with paths["sum"].open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("Repeat Classes\n==============\n")
        handle.write(f"Total Sequences: 8\nTotal Length: {genome_length} bp\n")
        handle.write("Class                  Count        bpMasked    %masked\n")
        handle.write("=====                  =====        ========     =======\n")
        for category in ordered_categories:
            entries = grouped.get(category, [])
            if not entries:
                continue
            if category in {"LINE", "SINE", "LTR", "TIR", "nonLTR", "nonTIR", "rDNA"}:
                handle.write(f"{category:<23} --           --           --\n")
            for subclass, count, bp in sorted(entries, key=lambda item: item[0].lower()):
                label = subclass if category in {"LINE", "SINE", "LTR", "TIR", "nonLTR", "nonTIR", "rDNA"} else category
                handle.write(f"    {label:<19} {count:<12} {bp:<12} {100 * bp / genome_length:.2f}%\n")
                if category not in {"LINE", "SINE", "LTR", "TIR", "nonLTR", "nonTIR", "rDNA"}:
                    break
        handle.write("                      ---------------------------------\n")
        handle.write(
            f"    total interspersed {interspersed_count:<12} {interspersed_bp:<12} "
            f"{100 * interspersed_bp / genome_length:.2f}%\n\n"
        )
        handle.write("---------------------------------------------------------\n")
        handle.write(
            f"Total                  {stats.total_count:<12} {stats.total_bp:<12} "
            f"{100 * stats.total_bp / genome_length:.2f}%\n\n"
        )
        handle.write("Extraction note\n===============\n")
        handle.write(
            "This summary was recalculated from the chromosome-subset TEanno.gff3. "
            "bpMasked is the union of TE intervals and is not a de novo EDTA run.\n"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-fasta", type=Path, required=True, help="unmasked 2n=32 FASTA used by EDTA")
    parser.add_argument("--source-masked", type=Path, required=True, help="2n=32 EDTA .mod.MAKER.masked FASTA")
    parser.add_argument("--source-gff", type=Path, required=True, help="2n=32 EDTA TEanno.gff3")
    parser.add_argument("--source-library", type=Path, required=True, help="2n=32 EDTA TElib.fa")
    parser.add_argument("--pzi-map", type=Path, required=True, help="current pzi_repartition_map.tsv")
    parser.add_argument("--input-genome-dir", type=Path, required=True, help="repartitioned pzi genome directory")
    parser.add_argument("--out-root", type=Path, required=True, help="annotation/repartition/repeats/pzi")
    parser.add_argument("--force", action="store_true", help="replace known extracted outputs")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    for path, label in (
        (args.source_fasta, "2n=32 source FASTA"),
        (args.source_masked, "2n=32 masked FASTA"),
        (args.source_gff, "2n=32 TE annotation GFF3"),
        (args.source_library, "2n=32 TE library"),
        (args.pzi_map, "PZI repartition map"),
    ):
        require_file(path, label)

    mappings = read_map(args.pzi_map)
    paths_by_hap = {hap: output_paths(args.out_root, hap) for hap in HAPLOTYPES}
    source_to_target: dict[str, tuple[str, str]] = {}
    target_genomes: dict[str, Path] = {}
    for hap in HAPLOTYPES:
        genome = args.input_genome_dir / f"pzi_{hap}.fa"
        require_file(genome, f"repartitioned PZI {hap} genome")
        target_genomes[hap] = genome
        for target, source in mappings[hap].items():
            source_to_target[source] = (hap, target)

    # Index the frozen 2n=32 source and require all selected sequences to be present.
    source_info: dict[str, tuple[str, int]] = {}
    for name, _, sequence in fasta_records(args.source_fasta):
        key = canonical_seqid(name)
        if key in source_info:
            die(f"Duplicate canonical chromosome ID in source FASTA: {key}")
        source_info[key] = (md5_text(sequence), len(sequence))
    missing_source = sorted(set(source_to_target) - set(source_info))
    if missing_source:
        die(f"Mapped source chromosomes missing from 2n=32 FASTA: {','.join(missing_source)}")

    # Verify the newly partitioned genomes contain the exact selected 2n=32 sequences.
    verification_rows: list[dict[str, str]] = []
    genome_lengths: dict[str, int] = {}
    for hap in HAPLOTYPES:
        expected = mappings[hap]
        seen: set[str] = set()
        genome_length = 0
        for target, _, sequence in fasta_records(target_genomes[hap]):
            if target not in expected:
                die(f"Unexpected chromosome {target} in {target_genomes[hap]}")
            source = expected[target]
            source_md5, source_length = source_info[source]
            target_md5 = md5_text(sequence)
            if target_md5 != source_md5 or len(sequence) != source_length:
                die(
                    f"Sequence mismatch for {hap}/{target}: target does not exactly match "
                    f"2n=32 source {source}; EDTA extraction is unsafe"
                )
            seen.add(target)
            genome_length += len(sequence)
            verification_rows.append(
                {
                    "haplotype": hap,
                    "new_seqid": target,
                    "original_chr32": source,
                    "length": str(len(sequence)),
                    "source_md5": source_md5,
                    "target_md5": target_md5,
                    "status": "pass",
                }
            )
        missing_target = sorted(set(expected) - seen)
        if missing_target or len(seen) != 8:
            die(f"Genome {target_genomes[hap]} does not contain its expected eight chromosomes: {missing_target}")
        genome_lengths[hap] = genome_length

    # Materialize the unmasked working FASTA links required by downstream LAI preparation.
    for hap, paths in paths_by_hap.items():
        paths["dir"].mkdir(parents=True, exist_ok=True)
        link_or_copy(target_genomes[hap], paths["genome"], args.force)

    # Extract the masked source sequences under the new CxA/CxB chromosome IDs.
    masked_handles = {}
    masked_seen: set[str] = set()
    try:
        for hap, paths in paths_by_hap.items():
            replace_output(paths["masked"], args.force)
            masked_handles[hap] = paths["masked"].open("w", encoding="ascii", newline="\n")
        for name, _, sequence in fasta_records(args.source_masked, uppercase=False):
            source = canonical_seqid(name)
            if source not in source_to_target:
                continue
            hap, target = source_to_target[source]
            _, source_length = source_info[source]
            if len(sequence) != source_length:
                die(f"Masked source sequence length differs from unmasked source for {source}")
            write_fasta_record(masked_handles[hap], f"{target} source_chr32={source}", sequence)
            masked_seen.add(source)
    finally:
        for handle in masked_handles.values():
            handle.close()
    missing_masked = sorted(set(source_to_target) - masked_seen)
    if missing_masked:
        die(f"Mapped source chromosomes missing from masked 2n=32 FASTA: {','.join(missing_masked)}")

    # Filter and rename TE features. The same pass gathers an interval-union summary.
    gff_handles = {}
    summaries = {hap: RepeatSummary() for hap in HAPLOTYPES}
    try:
        for hap, paths in paths_by_hap.items():
            replace_output(paths["gff"], args.force)
            handle = paths["gff"].open("w", encoding="utf-8", newline="\n")
            handle.write("##gff-version 3\n")
            handle.write("# Extracted from the frozen PZI 2n=32 EDTA annotation; sequence IDs were renamed only.\n")
            gff_handles[hap] = handle
        with args.source_gff.open("r", encoding="utf-8", errors="replace") as source_gff:
            for raw in source_gff:
                if raw.startswith("##FASTA"):
                    break
                if raw.startswith("##sequence-region"):
                    fields = raw.rstrip("\n").split()
                    if len(fields) >= 4:
                        source = canonical_seqid(fields[1])
                        if source in source_to_target:
                            hap, target = source_to_target[source]
                            gff_handles[hap].write(f"##sequence-region {target} {fields[2]} {fields[3]}\n")
                    continue
                if raw.startswith("#") or not raw.strip():
                    continue
                fields = raw.rstrip("\n").split("\t")
                if len(fields) != 9:
                    continue
                source = canonical_seqid(fields[0])
                if source not in source_to_target:
                    continue
                try:
                    start, end = int(fields[3]), int(fields[4])
                except ValueError as exc:
                    raise RuntimeError(f"Invalid EDTA feature coordinates: {raw.rstrip()}") from exc
                hap, target = source_to_target[source]
                category, subclass = classify_repeat(fields[8])
                summaries[hap].add(source, start, end, category, subclass)
                fields[0] = target
                gff_handles[hap].write("\t".join(fields) + "\n")
    finally:
        for handle in gff_handles.values():
            handle.close()

    for summary in summaries.values():
        summary.finish()

    all_manifest = args.out_root / "pzi_chr32_edta_extraction_manifest.tsv"
    replace_output(all_manifest, args.force)
    with all_manifest.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(verification_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(verification_rows)

    for hap, paths in paths_by_hap.items():
        write_summary(paths, summaries[hap], genome_lengths[hap])
        link_or_copy(args.source_library, paths["library"], args.force)
        replace_output(paths["chromosome_map"], args.force)
        with paths["chromosome_map"].open("w", encoding="utf-8", newline="\n") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(verification_rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows([row for row in verification_rows if row["haplotype"] == hap])
        replace_output(paths["completed"], args.force)
        with paths["completed"].open("w", encoding="utf-8", newline="\n") as handle:
            handle.write("species\tlabel\tgenome\tte_gff3\tte_summary\tte_library\tmasked_genome\tstatus\n")
            handle.write(
                f"pzi\t{hap}\t{paths['genome']}\t{paths['gff']}\t{paths['sum']}\t"
                f"{paths['library']}\t{paths['masked']}\tpass_extracted_from_chr32\n"
            )

    print(f"[OK] Extracted and sequence-verified PZI EDTA annotations: {args.out_root}")


if __name__ == "__main__":
    try:
        main()
    except RuntimeError as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        raise SystemExit(2)
