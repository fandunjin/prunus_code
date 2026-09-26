#!/usr/bin/env python3
"""Split current PZI, PSP and PPA complete-assembly annotations by haplotype.

The program reads the postprocessed EviAnn GFF3, CDS and protein FASTA files
from the complete PZI 2n=32, PSP 2n=32 and PPA 2n=16 assemblies.  It uses the
corresponding ``haploid_repartition/maps`` table to select the eight source
chromosomes assigned to each output haplotype, renames only the sequence IDs
to C1A--C8A or C1B--C8B, and retains original gene/transcript IDs for
traceability.  No gene prediction or sequence reconstruction is performed.
"""

import argparse
import csv
import hashlib
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Set, Tuple


DNA_COMPLEMENT = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")
STANDARD_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def parse_args() -> argparse.Namespace:
    project_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=project_root)
    parser.add_argument("--force", action="store_true", help="replace the known annotation outputs")
    parser.add_argument(
        "--species", nargs="+", choices=("pzi", "psp", "ppa"),
        help="species to extract; defaults to SPECIES_SET or all three species",
    )
    return parser.parse_args()


def die(message: str) -> None:
    raise RuntimeError(message)


def require_file(path: Path, label: str) -> Path:
    if not path.is_file() or path.stat().st_size == 0:
        die("Missing or empty {}: {}".format(label, path))
    return path


def require_first_file(paths: Sequence[Path], label: str) -> Path:
    for path in paths:
        if path.is_file() and path.stat().st_size > 0:
            return path
    die("Missing or empty {}. Checked: {}".format(label, "; ".join(str(path) for path in paths)))


def parse_attributes(text: str) -> Dict[str, str]:
    values: Dict[str, str] = {}
    for item in text.rstrip(";\n").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            values[key] = value
    return values


def fasta_records(path: Path) -> Iterator[Tuple[str, str]]:
    header = None
    chunks: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks).upper()
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks).upper()


def fasta_identifier(header: str) -> str:
    return header.split(None, 1)[0]


def read_fasta_lengths(path: Path) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    for header, sequence in fasta_records(path):
        name = fasta_identifier(header)
        if name in lengths:
            die("Duplicate sequence ID in {}: {}".format(path, name))
        lengths[name] = len(sequence)
    if not lengths:
        die("No FASTA records in {}".format(path))
    return lengths


def write_fasta_record(handle, header: str, sequence: str) -> None:
    handle.write(">" + header + "\n")
    for offset in range(0, len(sequence), 80):
        handle.write(sequence[offset:offset + 80] + "\n")


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def gff_sequence_ids(path: Path) -> Set[str]:
    identifiers: Set[str] = set()
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) == 9:
                identifiers.add(fields[0])
    if not identifiers:
        die("No feature sequence IDs in {}".format(path))
    return identifiers


def canonical_psp_seqid(value: str) -> str:
    """Normalize C01_A1 and C1A_A1 to the same PSP chromosome key."""
    value = value.strip()
    duplicate_subgenome = re.fullmatch(r"C0*(\d+)([AB])_?([AB])([12])", value, flags=re.IGNORECASE)
    if duplicate_subgenome is not None:
        group, subgenome, phase_subgenome, phase = duplicate_subgenome.groups()
        if subgenome.upper() != phase_subgenome.upper():
            die("PSP sequence ID has incompatible subgenome/phase: {}".format(value))
        return "C{}{}_{}{}".format(int(group), subgenome.upper(), phase_subgenome.upper(), phase)
    compact = re.fullmatch(r"C0*(\d+)_?([AB])([12])", value, flags=re.IGNORECASE)
    if compact is not None:
        group, subgenome, phase = compact.groups()
        return "C{}{}_{}{}".format(int(group), subgenome.upper(), subgenome.upper(), phase)
    return value


def psp_source_map_for_haplotype(
    psp_rows: Sequence[Mapping[str, str]],
    haplotype: str,
    source_gff: Path,
    source_genome: Path,
    target_genome: Path,
    final_rename_candidates: Sequence[Path],
    logical_rename_candidates: Sequence[Path],
) -> Tuple[Dict[str, str], Dict[str, str], Mapping[str, int], str]:
    """Resolve the PSP source IDs and required coordinate orientation.

    PSP 2n=32 EviAnn may use either the raw HapHiC ``scaffold_*`` IDs or
    orientation-corrected C01_A1-style IDs.  The latter are the IDs retained
    in the repartition map.  Raw scaffold records need coordinate/strand
    conversion when the target chromosome was reverse complemented.
    """
    selected = [row for row in psp_rows if row["haplotype"] == haplotype]
    if len(selected) != 8:
        die("PSP {}: expected eight map rows, found {}".format(haplotype, len(selected)))
    targets = {canonical_psp_seqid(row["old_seqid"]): row["new_seqid"] for row in selected}
    source_ids = gff_sequence_ids(source_gff)

    # Build source-GFF-ID -> logical C#A_A#/C#B_B# using the final assembly
    # rename table and the constrained PSP haplotype assignment table.
    source_to_logical: Dict[str, str] = {}
    for source_id in source_ids:
        canonical = canonical_psp_seqid(source_id)
        if canonical in targets:
            source_to_logical[source_id] = canonical
    if len(source_to_logical) != 8:
        final_map = require_first_file(final_rename_candidates, "PSP final chr32 rename map")
        logical_map = require_first_file(logical_rename_candidates, "PSP haplotype assignment map")
        final_rows = read_tsv(final_map)
        logical_rows = read_tsv(logical_map)
        final_to_raw = {row.get("new_name", ""): row.get("old_name", "") for row in final_rows}
        raw_to_logical = {
            row.get("old_seqid", ""): canonical_psp_seqid(row.get("new_seqid", ""))
            for row in logical_rows
            if row.get("old_seqid") and row.get("new_seqid")
        }
        source_to_logical = {}
        for source_id in source_ids:
            raw_id = final_to_raw.get(source_id, source_id)
            logical = raw_to_logical.get(raw_id)
            if logical in targets:
                source_to_logical[source_id] = logical
    if len(source_to_logical) != 8 or len(set(source_to_logical.values())) != 8:
        available = ",".join(sorted(source_ids)[:12])
        die("PSP {}: cannot resolve eight source GFF IDs through chr32 rename maps; first IDs: {}".format(haplotype, available))

    source_map = {source: targets[logical] for source, logical in source_to_logical.items()}
    source_records = {
        fasta_identifier(header): sequence
        for header, sequence in fasta_records(source_genome)
        if fasta_identifier(header) in source_map
    }
    source_lengths = {name: len(sequence) for name, sequence in source_records.items()}
    if set(source_records) != set(source_map):
        missing = sorted(set(source_map) - set(source_records))
        die("PSP {}: source genome lacks selected record(s): {}".format(haplotype, ",".join(missing)))
    target_records = {fasta_identifier(header): sequence for header, sequence in fasta_records(target_genome)}
    orientations: Dict[str, str] = {}
    for source_id, target_id in source_map.items():
        source_sequence = source_records.get(source_id)
        target_sequence = target_records.get(target_id)
        if source_sequence is None or target_sequence is None:
            die("PSP {}: source/target genome record missing for {} -> {}".format(haplotype, source_id, target_id))
        if source_sequence == target_sequence:
            orientations[source_id] = "forward"
        elif reverse_complement(source_sequence) == target_sequence:
            orientations[source_id] = "reverse_complement"
        else:
            die("PSP {}: target sequence is not identical or reverse-complemented for {} -> {}".format(haplotype, source_id, target_id))
    mode = "psp32_named" if all(canonical_psp_seqid(source) in targets for source in source_map) else "psp32_chr32_rename_verified"
    return source_map, orientations, source_lengths, mode


def ensure_output_path(path: Path, force: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if not force:
            die("Refusing to overwrite existing output (use --force): {}".format(path))
        path.unlink()


def relative(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def validate_selection(label: str, sequence_map: Mapping[str, str], target_lengths: Mapping[str, int]) -> None:
    if len(sequence_map) != 8:
        die("{}: expected 8 source chromosomes, found {}".format(label, len(sequence_map)))
    if len(set(sequence_map.values())) != 8:
        die("{}: target chromosome IDs are not one-to-one".format(label))
    missing = sorted(set(sequence_map.values()) - set(target_lengths))
    extra = sorted(set(target_lengths) - set(sequence_map.values()))
    if missing or extra:
        die("{}: target genome IDs disagree with map; missing={}, extra={}".format(label, missing, extra))


def selected_gff(
    source_gff: Path,
    sequence_map: Mapping[str, str],
    target_lengths: Mapping[str, int],
    output_gff: Path,
    force: bool,
    append: bool = False,
    source_orientations: Optional[Mapping[str, str]] = None,
    source_lengths: Optional[Mapping[str, int]] = None,
) -> Dict[str, object]:
    """Filter a GFF3 and return transcript-to-target chromosome metadata."""
    if append:
        if not output_gff.is_file() or output_gff.stat().st_size == 0:
            die("Cannot append GFF3 before creating it: {}".format(output_gff))
    else:
        ensure_output_path(output_gff, force)
    transcripts: Dict[str, Dict[str, str]] = {}
    source_seen: Set[str] = set()
    feature_counts: Counter = Counter()
    coordinate_errors: List[str] = []
    region_seen: Set[str] = set()

    with source_gff.open("r", encoding="utf-8", errors="replace") as source, output_gff.open(
        "a" if append else "w", encoding="utf-8", newline="\n"
    ) as destination:
        if not append:
            destination.write("##gff-version 3\n")
            destination.write("# Extracted for revised haploid_repartition; original feature IDs are retained.\n")
        for raw in source:
            if raw.startswith("##sequence-region"):
                fields = raw.rstrip("\n").split()
                if len(fields) >= 4 and fields[1] in sequence_map:
                    source_seqid = fields[1]
                    target_seqid = sequence_map[source_seqid]
                    destination.write("##sequence-region {} {} {}\n".format(target_seqid, fields[2], fields[3]))
                    region_seen.add(source_seqid)
                continue
            if raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            source_seqid = fields[0]
            if source_seqid not in sequence_map:
                continue
            source_seen.add(source_seqid)
            target_seqid = sequence_map[source_seqid]
            try:
                end = int(fields[4])
            except ValueError:
                coordinate_errors.append("{}:invalid_end".format(source_seqid))
                continue
            if end > target_lengths[target_seqid]:
                coordinate_errors.append("{}:{}>{}".format(source_seqid, end, target_lengths[target_seqid]))
                continue
            orientation = (source_orientations or {}).get(source_seqid, "forward")
            if orientation not in {"forward", "reverse_complement"}:
                die("{}: unsupported orientation for {}: {}".format(source_gff, source_seqid, orientation))
            if orientation == "reverse_complement":
                source_length = (source_lengths or {}).get(source_seqid, target_lengths[target_seqid])
                if source_length != target_lengths[target_seqid]:
                    coordinate_errors.append(
                        "{}:length {} != target {}".format(source_seqid, source_length, target_lengths[target_seqid])
                    )
                    continue
                start = int(fields[3])
                fields[3] = str(source_length - end + 1)
                fields[4] = str(source_length - start + 1)
                if fields[6] == "+":
                    fields[6] = "-"
                elif fields[6] == "-":
                    fields[6] = "+"
            fields[0] = target_seqid
            feature = fields[2].lower()
            attrs = parse_attributes(fields[8])
            if feature in {"mrna", "transcript"}:
                transcript_id = attrs.get("ID")
                if transcript_id:
                    transcripts[transcript_id] = {
                        "source_seqid": source_seqid,
                        "target_seqid": target_seqid,
                        "gene_id": attrs.get("Parent", transcript_id).split(",")[0],
                        "strand": fields[6],
                        "source_orientation": orientation,
                    }
            feature_counts[feature] += 1
            destination.write("\t".join(fields) + "\n")

    missing_source = sorted(set(sequence_map) - source_seen)
    if missing_source:
        die("{}: selected source chromosomes absent from GFF3: {}".format(source_gff, ",".join(missing_source)))
    if coordinate_errors:
        die("{}: output coordinates exceed target genome: {}".format(source_gff, ";".join(coordinate_errors[:10])))
    if not transcripts:
        die("{}: no selected mRNA/transcript features".format(source_gff))
    return {
        "transcripts": transcripts,
        "feature_counts": feature_counts,
        "source_seen": source_seen,
        "region_seen": region_seen,
    }


def rewritten_header(header: str, target_seqid: str, target_strand: str) -> str:
    if re.search(r"(?:^|\s)seqid=[^\s;]+", header):
        header = re.sub(r"(?<!\S)seqid=[^\s;]+", "seqid=" + target_seqid, header, count=1)
    else:
        header = header + " seqid=" + target_seqid
    if target_strand in {"+", "-"}:
        if re.search(r"(?:^|\s)strand=[+-](?=\s|$)", header):
            header = re.sub(r"(?<!\S)strand=[+-](?=\s|$)", "strand=" + target_strand, header, count=1)
        else:
            header = header + " strand=" + target_strand
    return header


def selected_fasta(
    source_fasta: Path,
    transcript_info: Mapping[str, Mapping[str, str]],
    output_fasta: Path,
    force: bool,
    append: bool = False,
) -> Set[str]:
    if append:
        if not output_fasta.is_file() or output_fasta.stat().st_size == 0:
            die("Cannot append FASTA before creating it: {}".format(output_fasta))
    else:
        ensure_output_path(output_fasta, force)
    wanted = set(transcript_info)
    found: Set[str] = set()
    with output_fasta.open("a" if append else "w", encoding="utf-8", newline="\n") as destination:
        for header, sequence in fasta_records(source_fasta):
            transcript_id = fasta_identifier(header)
            if transcript_id not in wanted:
                continue
            if transcript_id in found:
                die("Duplicate transcript {} in {}".format(transcript_id, source_fasta))
            if not sequence:
                continue
            write_fasta_record(
                destination,
                rewritten_header(
                    header,
                    transcript_info[transcript_id]["target_seqid"],
                    transcript_info[transcript_id]["strand"],
                ),
                sequence,
            )
            found.add(transcript_id)
    return found


def filter_fasta_to_identifiers(path: Path, retained_identifiers: Set[str]) -> int:
    """Ensure an extracted FASTA contains exactly the model set in its GFF3."""
    temporary = path.with_name(path.name + ".filtering_tmp")
    found: Set[str] = set()
    with temporary.open("w", encoding="utf-8", newline="\n") as destination:
        for header, sequence in fasta_records(path):
            identifier = fasta_identifier(header)
            if identifier not in retained_identifiers:
                continue
            if identifier in found:
                die("Duplicate transcript {} in {}".format(identifier, path))
            write_fasta_record(destination, header, sequence)
            found.add(identifier)
    if found != retained_identifiers:
        missing = sorted(retained_identifiers - found)
        die("{}: retained transcript(s) absent from output FASTA: {}".format(path, ",".join(missing[:10])))
    temporary.replace(path)
    return len(found)


def prune_gff_to_transcripts(path: Path, retained_transcripts: Set[str]) -> Tuple[int, int]:
    """Keep only gene models represented by both output CDS and protein FASTA."""
    retained_genes: Set[str] = set()
    with path.open("r", encoding="utf-8", errors="replace") as source:
        for raw in source:
            if raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2].lower() not in {"mrna", "transcript"}:
                continue
            attrs = parse_attributes(fields[8])
            if attrs.get("ID") in retained_transcripts:
                retained_genes.update(item for item in attrs.get("Parent", "").split(",") if item)
    temporary = path.with_name(path.name + ".pruning_tmp")
    retained_mrnas = 0
    retained_gene_records = 0
    with path.open("r", encoding="utf-8", errors="replace") as source, temporary.open(
        "w", encoding="utf-8", newline="\n"
    ) as destination:
        for raw in source:
            if raw.startswith("#"):
                destination.write(raw)
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            feature = fields[2].lower()
            attrs = parse_attributes(fields[8])
            keep = False
            if feature == "gene":
                keep = attrs.get("ID") in retained_genes
                retained_gene_records += int(keep)
            elif feature in {"mrna", "transcript"}:
                keep = attrs.get("ID") in retained_transcripts
                retained_mrnas += int(keep)
            else:
                parents = set(item for item in attrs.get("Parent", "").split(",") if item)
                keep = bool(parents.intersection(retained_transcripts))
            if keep:
                destination.write("\t".join(fields) + "\n")
    temporary.replace(path)
    return retained_gene_records, retained_mrnas


def fasta_md5_index(path: Path) -> Dict[str, str]:
    """Index source sequence IDs by exact full-length sequence MD5."""
    index: Dict[str, str] = {}
    for header, sequence in fasta_records(path):
        digest = hashlib.md5(sequence.encode("ascii")).hexdigest()
        source_id = fasta_identifier(header)
        if digest in index:
            die("Non-unique source sequence MD5 in {}: {} and {}".format(path, index[digest], source_id))
        index[digest] = source_id
    return index


def psp_sequence_map_from_md5(
    target_genome: Path, source_md5: Mapping[str, str], haplotype: str
) -> Dict[str, str]:
    """Return old-PSP-seqid -> revised-seqid, verified against sequence content.

    The revised PSP packages deliberately recombine chromosomes from the old
    A/B haplotype files.  Their historical ``old_seqid`` column is a PZI
    display label, not an annotation sequence identifier, so only exact
    sequence matching can safely select the appropriate old annotation.
    """
    expected_phase = haplotype.replace("hap", "")[-1]
    sequence_map: Dict[str, str] = {}
    for header, sequence in fasta_records(target_genome):
        target_id = fasta_identifier(header)
        digest = hashlib.md5(sequence.encode("ascii")).hexdigest()
        source_id = source_md5.get(digest)
        if source_id is None:
            die("{}: {} has no exact match in original PSP chr32 FASTA".format(target_genome, target_id))
        match = re.match(r"^C\d+_([AB])([12])$", source_id)
        if match is None or match.group(2) != expected_phase:
            die("{}: {} matched incompatible old PSP source {}".format(target_genome, target_id, source_id))
        if source_id in sequence_map:
            die("{}: old PSP source {} occurs more than once".format(target_genome, source_id))
        sequence_map[source_id] = target_id
    if len(sequence_map) != 8:
        die("{}: expected 8 exact PSP chromosome matches, found {}".format(target_genome, len(sequence_map)))
    return sequence_map


def reverse_complement(sequence: str) -> str:
    return sequence.translate(DNA_COMPLEMENT)[::-1]


def translate_standard(coding: str) -> str:
    amino = []
    for offset in range(0, len(coding), 3):
        amino.append(STANDARD_CODE.get(coding[offset:offset + 3], "X"))
    return "".join(amino)


def selected_cds_segments(source_gff: Path, wanted_transcripts: Set[str]) -> Dict[str, List[Tuple[int, int, str, str]]]:
    segments: Dict[str, List[Tuple[int, int, str, str]]] = defaultdict(list)
    with source_gff.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2].lower() != "cds":
                continue
            attrs = parse_attributes(fields[8])
            parents = [item for item in attrs.get("Parent", "").split(",") if item in wanted_transcripts]
            if not parents:
                continue
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                continue
            for parent in parents:
                segments[parent].append((start, end, fields[6], fields[7]))
    return segments


def extract_copy2_fastas(
    source_gff: Path,
    source_genome: Path,
    sequence_map: Mapping[str, str],
    transcript_info: Mapping[str, Mapping[str, str]],
    output_cds: Path,
    output_pep: Path,
    force: bool,
) -> Tuple[Set[str], Counter]:
    """Rebuild copy2 CDS/proteins with one source chromosome resident at once."""
    ensure_output_path(output_cds, force)
    ensure_output_path(output_pep, force)
    segments = selected_cds_segments(source_gff, set(transcript_info))
    by_chromosome: Dict[str, List[str]] = defaultdict(list)
    for transcript_id, info in transcript_info.items():
        by_chromosome[info["source_seqid"]].append(transcript_id)
    for chromosome, ids in by_chromosome.items():
        ids.sort(key=lambda item: (min(segment[0] for segment in segments.get(item, [(0, 0, "+", "0")])), item))

    found_chromosomes: Set[str] = set()
    written: Set[str] = set()
    skipped: Counter = Counter()
    with output_cds.open("w", encoding="utf-8", newline="\n") as cds_handle, output_pep.open(
        "w", encoding="utf-8", newline="\n"
    ) as pep_handle:
        for header, genomic_sequence in fasta_records(source_genome):
            chromosome = fasta_identifier(header)
            if chromosome not in sequence_map:
                continue
            found_chromosomes.add(chromosome)
            for transcript_id in by_chromosome.get(chromosome, []):
                transcript_segments = segments.get(transcript_id, [])
                if not transcript_segments:
                    skipped["no_cds"] += 1
                    continue
                strand = transcript_info[transcript_id]["strand"]
                ordered = sorted(transcript_segments, key=lambda row: (row[0], row[1]), reverse=(strand == "-"))
                parts: List[str] = []
                invalid = False
                coding_length_before_segment = 0
                for start, end, segment_strand, phase in ordered:
                    if segment_strand != strand or start < 1 or end < start or end > len(genomic_sequence):
                        invalid = True
                        break
                    fragment = genomic_sequence[start - 1:end]
                    if strand == "-":
                        fragment = reverse_complement(fragment)
                    try:
                        phase_value = int(phase) if phase != "." else 0
                    except ValueError:
                        invalid = True
                        break
                    if phase_value not in {0, 1, 2}:
                        invalid = True
                        break
                    # GFF3 phase describes how this segment joins the prior
                    # segment.  Its leading bases complete a codon and must
                    # remain in the concatenated CDS; they are not trimmed.
                    expected_phase = (3 - (coding_length_before_segment % 3)) % 3
                    if phase != "." and phase_value != expected_phase:
                        skipped["phase_discordant_segments"] += 1
                    parts.append(fragment)
                    coding_length_before_segment += len(fragment)
                if invalid:
                    skipped["invalid_feature"] += 1
                    continue
                coding = "".join(parts).upper()
                if len(coding) < 3 or len(coding) % 3:
                    skipped["non_triplet_cds"] += 1
                    continue
                protein = translate_standard(coding)
                if protein.endswith("*"):
                    protein = protein[:-1]
                if not protein or "*" in protein:
                    skipped["internal_stop_or_empty_protein"] += 1
                    continue
                info = transcript_info[transcript_id]
                header_text = "{} gene={} seqid={} strand={} cds_len={} source=reconstructed_from_copy2_gff3".format(
                    transcript_id, info["gene_id"], info["target_seqid"], strand, len(coding)
                )
                write_fasta_record(cds_handle, header_text, coding)
                write_fasta_record(pep_handle, header_text, protein)
                written.add(transcript_id)
    missing_chromosomes = sorted(set(sequence_map) - found_chromosomes)
    if missing_chromosomes:
        die("{}: selected chromosomes absent from copy2 genome: {}".format(source_genome, ",".join(missing_chromosomes)))
    return written, skipped


def map_for_haplotype(rows: Iterable[Mapping[str, str]], haplotype: str, source_column: str) -> Dict[str, str]:
    selected = [row for row in rows if row["haplotype"] == haplotype]
    sequence_map: Dict[str, str] = {}
    for row in selected:
        source = row[source_column]
        target = row["new_seqid"]
        if source in sequence_map and sequence_map[source] != target:
            die("{} maps to more than one target chromosome".format(source))
        sequence_map[source] = target
    return sequence_map


def output_paths(root: Path, species: str, haplotype: str) -> Dict[str, Path]:
    base = root / "haploid_repartition" / "inputs" / species
    prefix = "{}_{}".format(species, haplotype)
    return {
        "genome": base / "genome" / (prefix + ".fa"),
        "gff": base / "gff" / (prefix + ".gff3"),
        "cds": base / "cds" / (prefix + ".cds.fa"),
        "pep": base / "pep" / (prefix + ".pep.fa"),
    }


def write_transcript_rows(
    output_rows: List[Dict[str, str]], species: str, haplotype: str, annotation_mode: str, transcript_info: Mapping[str, Mapping[str, str]]
) -> None:
    for transcript_id in sorted(transcript_info):
        info = transcript_info[transcript_id]
        output_rows.append({
            "species": species,
            "haplotype": haplotype,
            "annotation_mode": annotation_mode,
            "source_seqid": info["source_seqid"],
            "target_seqid": info["target_seqid"],
            "gene_id": info["gene_id"],
            "transcript_id": transcript_id,
        })


def run_direct_export(
    root: Path,
    species: str,
    haplotype: str,
    sequence_map: Mapping[str, str],
    source_gff: Path,
    source_cds: Path,
    source_pep: Path,
    force: bool,
    manifest: List[Dict[str, str]],
    transcript_rows: List[Dict[str, str]],
) -> None:
    paths = output_paths(root, species, haplotype)
    target_lengths = read_fasta_lengths(paths["genome"])
    validate_selection("{}_{}".format(species, haplotype), sequence_map, target_lengths)
    result = selected_gff(source_gff, sequence_map, target_lengths, paths["gff"], force)
    transcript_info = result["transcripts"]
    cds_ids = selected_fasta(source_cds, transcript_info, paths["cds"], force)
    pep_ids = selected_fasta(source_pep, transcript_info, paths["pep"], force)
    retained_ids = cds_ids.intersection(pep_ids)
    if not retained_ids:
        die("{}_{}: no transcript has both a CDS and protein sequence".format(species, haplotype))
    output_cds_count = filter_fasta_to_identifiers(paths["cds"], retained_ids)
    output_pep_count = filter_fasta_to_identifiers(paths["pep"], retained_ids)
    retained_gene_count, retained_transcript_count = prune_gff_to_transcripts(paths["gff"], retained_ids)
    if retained_transcript_count != output_cds_count or output_cds_count != output_pep_count:
        die(
            "{}_{}: GFF3/CDS/protein transcript counts disagree: {}/{}/{}".format(
                species, haplotype, retained_transcript_count, output_cds_count, output_pep_count
            )
        )
    write_transcript_rows(
        transcript_rows, species, haplotype, "source_fasta", {item: transcript_info[item] for item in retained_ids}
    )
    manifest.append({
        "species": species,
        "haplotype": haplotype,
        "annotation_mode": "source_fasta",
        "source_gff": relative(source_gff, root),
        "source_cds": relative(source_cds, root),
        "source_pep": relative(source_pep, root),
        "source_genome": "",
        "output_gff": relative(paths["gff"], root),
        "output_cds": relative(paths["cds"], root),
        "output_pep": relative(paths["pep"], root),
        "chromosomes": str(len(sequence_map)),
        "gff_genes": str(retained_gene_count),
        "gff_transcripts": str(retained_transcript_count),
        "source_gff_transcripts": str(len(transcript_info)),
        "cds_records": str(output_cds_count),
        "pep_records": str(output_pep_count),
        "missing_cds": str(len(set(transcript_info) - cds_ids)),
        "missing_pep": str(len(set(transcript_info) - pep_ids)),
        "gff_models_removed": str(len(set(transcript_info) - retained_ids)),
        "invalid_cds": "0",
        "status": "pass" if retained_ids == set(transcript_info) else "pass_after_source_model_filter",
    })


def run_direct_export_multi(
    root: Path,
    species: str,
    haplotype: str,
    sequence_map: Mapping[str, str],
    source_groups: Sequence[Tuple[Path, Path, Path, Mapping[str, str]]],
    force: bool,
    manifest: List[Dict[str, str]],
    transcript_rows: List[Dict[str, str]],
) -> None:
    """Export one new package from multiple old annotation packages."""
    paths = output_paths(root, species, haplotype)
    target_lengths = read_fasta_lengths(paths["genome"])
    validate_selection("{}_{}".format(species, haplotype), sequence_map, target_lengths)
    combined_transcripts: Dict[str, Dict[str, str]] = {}
    all_cds: Set[str] = set()
    all_pep: Set[str] = set()
    source_gffs: List[str] = []
    source_cdss: List[str] = []
    source_peps: List[str] = []
    for index, (source_gff, source_cds, source_pep, group_map) in enumerate(source_groups):
        result = selected_gff(source_gff, group_map, target_lengths, paths["gff"], force, append=(index > 0))
        transcript_info = result["transcripts"]
        duplicates = set(combined_transcripts).intersection(transcript_info)
        if duplicates:
            die("{} {} has duplicate old transcript IDs: {}".format(species, haplotype, ",".join(sorted(duplicates)[:5])))
        cds_ids = selected_fasta(source_cds, transcript_info, paths["cds"], force, append=(index > 0))
        pep_ids = selected_fasta(source_pep, transcript_info, paths["pep"], force, append=(index > 0))
        combined_transcripts.update(transcript_info)
        all_cds.update(cds_ids)
        all_pep.update(pep_ids)
        source_gffs.append(relative(source_gff, root))
        source_cdss.append(relative(source_cds, root))
        source_peps.append(relative(source_pep, root))
    retained_ids = all_cds.intersection(all_pep)
    retained_gene_count, retained_transcript_count = prune_gff_to_transcripts(paths["gff"], retained_ids)
    write_transcript_rows(
        transcript_rows, species, haplotype, "source_fasta_md5_verified",
        {item: combined_transcripts[item] for item in retained_ids},
    )
    manifest.append({
        "species": species,
        "haplotype": haplotype,
        "annotation_mode": "source_fasta_md5_verified",
        "source_gff": ";".join(source_gffs),
        "source_cds": ";".join(source_cdss),
        "source_pep": ";".join(source_peps),
        "source_genome": "data/psp/haploid/genome/psp_A1A2B1B2.renamed.chr32.fa",
        "output_gff": relative(paths["gff"], root),
        "output_cds": relative(paths["cds"], root),
        "output_pep": relative(paths["pep"], root),
        "chromosomes": str(len(sequence_map)),
        "gff_genes": str(retained_gene_count),
        "gff_transcripts": str(retained_transcript_count),
        "source_gff_transcripts": str(len(combined_transcripts)),
        "cds_records": str(len(all_cds)),
        "pep_records": str(len(all_pep)),
        "missing_cds": str(len(set(combined_transcripts) - all_cds)),
        "missing_pep": str(len(set(combined_transcripts) - all_pep)),
        "gff_models_removed": str(len(set(combined_transcripts) - retained_ids)),
        "invalid_cds": "0",
        "status": "pass" if retained_ids == set(combined_transcripts) else "pass_after_source_model_filter",
    })


def run_copy2_export(
    root: Path,
    haplotype: str,
    sequence_map: Mapping[str, str],
    source_gff: Path,
    source_genome: Path,
    force: bool,
    manifest: List[Dict[str, str]],
    transcript_rows: List[Dict[str, str]],
) -> None:
    paths = output_paths(root, "pzi", haplotype)
    target_lengths = read_fasta_lengths(paths["genome"])
    validate_selection("pzi_{}".format(haplotype), sequence_map, target_lengths)
    result = selected_gff(source_gff, sequence_map, target_lengths, paths["gff"], force)
    transcript_info = result["transcripts"]
    record_ids, skipped = extract_copy2_fastas(
        source_gff, source_genome, sequence_map, transcript_info, paths["cds"], paths["pep"], force
    )
    retained_gene_count, retained_transcript_count = prune_gff_to_transcripts(paths["gff"], record_ids)
    write_transcript_rows(
        transcript_rows, "pzi", haplotype, "reconstructed_copy2", {item: transcript_info[item] for item in record_ids}
    )
    manifest.append({
        "species": "pzi",
        "haplotype": haplotype,
        "annotation_mode": "reconstructed_copy2",
        "source_gff": relative(source_gff, root),
        "source_cds": "",
        "source_pep": "",
        "source_genome": relative(source_genome, root),
        "output_gff": relative(paths["gff"], root),
        "output_cds": relative(paths["cds"], root),
        "output_pep": relative(paths["pep"], root),
        "chromosomes": str(len(sequence_map)),
        "gff_genes": str(retained_gene_count),
        "gff_transcripts": str(retained_transcript_count),
        "source_gff_transcripts": str(len(transcript_info)),
        "cds_records": str(len(record_ids)),
        "pep_records": str(len(record_ids)),
        "missing_cds": str(len(set(transcript_info) - record_ids)),
        "missing_pep": str(len(set(transcript_info) - record_ids)),
        "gff_models_removed": str(len(set(transcript_info) - record_ids)),
        "invalid_cds": ";".join("{}={}".format(key, skipped[key]) for key in sorted(skipped)),
        "status": "pass" if len(record_ids) else "failed_no_valid_copy2_records",
    })


def write_table(path: Path, rows: Sequence[Mapping[str, str]], fieldnames: Sequence[str], force: bool) -> None:
    ensure_output_path(path, force)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    maps_dir = root / "haploid_repartition" / "maps"
    requested_species = args.species
    if requested_species is None:
        requested_species = os.environ.get("SPECIES_SET", "pzi psp ppa").split()
    requested_species = list(dict.fromkeys(requested_species))
    if not requested_species:
        die("No species selected; choose one or more of pzi, psp or ppa")
    species_set = set(requested_species)
    pzi_map = read_tsv(require_file(maps_dir / "pzi_repartition_map.tsv", "PZI repartition map")) if "pzi" in species_set else []
    psp_map = read_tsv(require_file(maps_dir / "psp_repartition_map.tsv", "PSP repartition map")) if "psp" in species_set else []
    ppa_map = read_tsv(require_file(maps_dir / "ppa_repartition_map.tsv", "PPA repartition map")) if "ppa" in species_set else []
    manifest: List[Dict[str, str]] = []
    transcript_rows: List[Dict[str, str]] = []

    def annotation_triplet(candidates: Sequence[Path], label: str) -> Tuple[Path, Path, Path]:
        """Resolve one current postprocessed annotation triplet by GFF path."""
        gff = require_first_file(candidates, label + " GFF3")
        stem = gff.name[:-len(".genes.gff3")]
        base = gff.parent
        cds = require_file(base / (stem + ".cds.fa"), label + " CDS")
        pep = require_file(base / (stem + ".protein.fa"), label + " protein")
        return gff, cds, pep

    # These are the current complete-assembly postprocess products.  Retired
    # copy1/copy2 and per-haplotype paths are deliberately not candidates.
    pzi_gff = pzi_cds = pzi_pep = None
    psp_gff = psp_cds = psp_pep = None
    ppa_gff = ppa_cds = ppa_pep = None
    if "pzi" in species_set:
        pzi_gff, pzi_cds, pzi_pep = annotation_triplet(
            [root / "annotation" / "genes" / "pzi_eviann" / "postprocess" / "pzi32.eviann.genes.gff3"],
            "PZI 2n=32 EviAnn",
        )
    if "psp" in species_set:
        psp_gff, psp_cds, psp_pep = annotation_triplet(
            [root / "data" / "psp" / "annotation" / "genes" / "psp32_eviann" / "postprocess" / "psp32.eviann.genes.gff3"],
            "PSP 2n=32 EviAnn",
        )
    if "ppa" in species_set:
        ppa_gff, ppa_cds, ppa_pep = annotation_triplet(
            [
                root / "data" / "ppa" / "annotation" / "eviann" / "postprocess" / "ppa.eviann.genes.gff3",
                root / "data" / "ppa" / "postprocess" / "ppa.eviann.genes.gff3",
            ],
            "PPA 2n=16 EviAnn",
        )

    # PZI full assembly IDs are chr01--chr32; the repartition map provides
    # the exact chr32 source for each target C1A--C8B chromosome.
    if "pzi" in species_set:
        for haplotype in ("hapA1", "hapA2", "hapB1", "hapB2"):
            run_direct_export(
                root, "pzi", haplotype, map_for_haplotype(pzi_map, haplotype, "original_chr32"),
                pzi_gff, pzi_cds, pzi_pep, args.force, manifest, transcript_rows,
            )

    # PSP 2n=32 annotation can be keyed by C01_A1-style IDs or raw scaffold
    # IDs. Resolve and verify against the actual revised haplotype FASTA so
    # reverse-complemented chromosomes receive the corresponding GFF3
    # coordinate/strand liftover.
    if "psp" in species_set:
        psp_source_genome = require_file(
            root / "data" / "psp" / "haphic" / "merge" / "04.build" / "final_assembly_qc_outputs" / "psp_2n32.chromosomes.fa",
            "PSP source 2n=32 chromosome FASTA",
        )
        psp_final_rename = [
            root / "data" / "psp" / "haphic" / "merge" / "04.build" / "final_assembly_qc_outputs" / "psp_2n32.rename_map.tsv",
            root / "data" / "psp" / "haphic" / "merge" / "04.build" / "psp_2n32.rename_map.tsv",
        ]
        psp_logical_rename = [root / "data" / "psp" / "haploid" / "psp_A1A2B1B2.rename_map.tsv"]
        for haplotype in ("hapA1", "hapA2", "hapB1", "hapB2"):
            target_genome = output_paths(root, "psp", haplotype)["genome"]
            psp_seq_map, psp_orientations, psp_source_lengths, psp_mode = psp_source_map_for_haplotype(
                psp_map, haplotype, psp_gff, psp_source_genome, target_genome, psp_final_rename, psp_logical_rename
            )
            result = selected_gff(
                psp_gff,
                psp_seq_map,
                read_fasta_lengths(target_genome),
                output_paths(root, "psp", haplotype)["gff"],
                args.force,
                source_orientations=psp_orientations,
                source_lengths=psp_source_lengths,
            )
        # selected_gff has already written the orientation-corrected GFF3;
        # export the corresponding CDS/protein records and apply the same
        # common-model filtering as the other species.
            paths = output_paths(root, "psp", haplotype)
            transcript_info = result["transcripts"]
        # CDS/protein FASTA records are already in transcript 5'-to-3'
        # orientation. A genomic reverse-complement liftover changes the GFF3
        # coordinates and strand, but not the spliced CDS or peptide sequence.
            cds_ids = selected_fasta(psp_cds, transcript_info, paths["cds"], args.force)
            pep_ids = selected_fasta(psp_pep, transcript_info, paths["pep"], args.force)
            retained_ids = cds_ids.intersection(pep_ids)
            if not retained_ids:
                die("PSP {}: no transcript has both a CDS and protein sequence".format(haplotype))
            output_cds_count = filter_fasta_to_identifiers(paths["cds"], retained_ids)
            output_pep_count = filter_fasta_to_identifiers(paths["pep"], retained_ids)
            retained_gene_count, retained_transcript_count = prune_gff_to_transcripts(paths["gff"], retained_ids)
            if retained_transcript_count != output_cds_count or output_cds_count != output_pep_count:
                die("PSP {}: GFF3/CDS/protein transcript counts disagree".format(haplotype))
            write_transcript_rows(transcript_rows, "psp", haplotype, psp_mode, {item: transcript_info[item] for item in retained_ids})
            manifest.append({
                "species": "psp", "haplotype": haplotype, "annotation_mode": psp_mode,
                "source_gff": relative(psp_gff, root), "source_cds": relative(psp_cds, root),
                "source_pep": relative(psp_pep, root), "source_genome": relative(psp_source_genome, root),
                "output_gff": relative(paths["gff"], root), "output_cds": relative(paths["cds"], root),
                "output_pep": relative(paths["pep"], root), "chromosomes": str(len(psp_seq_map)),
                "gff_genes": str(retained_gene_count), "gff_transcripts": str(retained_transcript_count),
                "source_gff_transcripts": str(len(transcript_info)), "cds_records": str(output_cds_count),
                "pep_records": str(output_pep_count), "missing_cds": str(len(set(transcript_info) - cds_ids)),
                "missing_pep": str(len(set(transcript_info) - pep_ids)), "gff_models_removed": str(len(set(transcript_info) - retained_ids)),
                "invalid_cds": "0", "status": "pass" if retained_ids == set(transcript_info) else "pass_after_source_model_filter",
            })
    # PPA current annotation contains the OZ060... primary chromosomes plus
    # unplaced contigs; the map selects exactly the eight chromosomes per A/B.
    if "ppa" in species_set:
        for haplotype in ("hapA", "hapB"):
            run_direct_export(
                root, "ppa", haplotype, map_for_haplotype(ppa_map, haplotype, "source_seqid"),
                ppa_gff, ppa_cds, ppa_pep, args.force, manifest, transcript_rows,
            )

    maps_dir = root / "haploid_repartition" / "maps"
    suffix = "" if species_set == {"pzi", "psp", "ppa"} else "_" + "_".join(requested_species)
    manifest_path = root / "haploid_repartition" / "inputs" / ("annotation_extraction_manifest" + suffix + ".tsv")
    transcript_path = maps_dir / ("annotation_extraction_transcript_map" + suffix + ".tsv")
    write_table(
        manifest_path, manifest,
        [
            "species", "haplotype", "annotation_mode", "source_gff", "source_cds", "source_pep", "source_genome",
            "output_gff", "output_cds", "output_pep", "chromosomes", "gff_genes", "gff_transcripts",
            "source_gff_transcripts", "cds_records", "pep_records", "missing_cds", "missing_pep",
            "gff_models_removed", "invalid_cds", "status",
        ], args.force,
    )
    write_table(
        transcript_path, transcript_rows,
        ["species", "haplotype", "annotation_mode", "source_seqid", "target_seqid", "gene_id", "transcript_id"], args.force,
    )
    complete = root / "haploid_repartition" / "inputs" / ("ANNOTATIONS_COMPLETE" if not suffix else "ANNOTATIONS_" + "_".join(requested_species).upper() + "_COMPLETE")
    ensure_output_path(complete, args.force)
    complete.write_text(
        "Annotations were extracted from the current complete-assembly EviAnn postprocess outputs.\n"
        "Selected complete-assembly GFF3/CDS/protein files were partitioned\n"
        "according to the corresponding haploid_repartition map tables.\n",
        encoding="ascii",
    )
    print("[OK] Species selected: {}".format(",".join(requested_species)))
    print("[OK] Extracted revised annotation packages: {}".format(root / "haploid_repartition" / "inputs"))
    print("[OK] Audit manifest: {}".format(manifest_path))
    for row in manifest:
        print(
            "[OK] {species}_{haplotype}: GFF transcripts={gff_transcripts}; CDS={cds_records}; PEP={pep_records}; status={status}".format(
                **row
            )
        )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as error:
        print("[ERROR] {}".format(error), file=sys.stderr)
        raise SystemExit(2)
