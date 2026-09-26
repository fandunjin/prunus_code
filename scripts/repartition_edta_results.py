#!/usr/bin/env python3
"""Repartition existing EDTA annotations according to revised haplotypes.

The script verifies chromosome identity by exact MD5, accepting either the
forward sequence or its reverse complement.  When a chromosome was reversed
after EDTA was run, masked FASTA, GFF3 coordinates/strands, and intact-element
FASTA records are transformed consistently.  It deliberately does not run
EDTA.
"""

import argparse
import csv
import hashlib
import os
import re
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterator, List, Mapping, Optional, Sequence, Set, Tuple


ROOT = Path(__file__).resolve().parents[1]
HAPS = {
    "pzi": ("hapA1", "hapA2", "hapB1", "hapB2"),
    "ppa": ("hapA", "hapB"),
    "psp": ("hapA1", "hapA2", "hapB1", "hapB2"),
}

# Used when a repartitioned haplotype chromosome is the reverse complement
# of the chromosome that was supplied to EDTA.  Keep case so soft-masking in
# the EDTA masked FASTA is preserved.
DNA_COMPLEMENT = str.maketrans(
    "ACGTRYKMSWBDHVNacgtrykmswbdhvn",
    "TGCAYRMKSWVHDBNtgcayrmkswvhdbn",
)


def reverse_complement(sequence):
    return sequence.translate(DNA_COMPLEMENT)[::-1]


def sequence_digest(sequence):
    """MD5 used for sequence identity (case-insensitive)."""
    return hashlib.md5(sequence.upper().encode("ascii")).hexdigest()


def oriented_match(
    target_record,
    source_records,
    source_by_digest=None,
    source_by_reverse_digest=None,
):
    """Return (source_record, orientation) for an exact match.

    The returned orientation describes how the target is related to the EDTA
    source chromosome.  A digest collision or an ambiguous forward/reverse
    match is rejected rather than silently assigning an annotation.
    """
    target_digest = target_record.digest
    if source_by_digest is None:
        forward = [record for record in source_records if record.digest == target_digest]
    else:
        record = source_by_digest.get(target_digest)
        forward = [record] if record is not None else []
    if source_by_reverse_digest is None:
        reverse_digest = sequence_digest(reverse_complement(target_record.sequence))
        reverse = [record for record in source_records if record.digest == reverse_digest]
    else:
        # The index key is digest(reverse_complement(source)); for a target
        # that is reverse-complemented, this key equals digest(target).
        reverse = list(source_by_reverse_digest.get(target_digest, ()))
    # The reverse digest above is the source digest when target is the reverse
    # complement of source.  The target itself can also be a forward match.
    candidates = [(record, "forward") for record in forward]
    candidates.extend((record, "reverse_complement") for record in reverse
                      if record not in forward)
    if len(candidates) > 1:
        names = ", ".join(f"{record.seqid} ({orientation})" for record, orientation in candidates)
        die(f"Ambiguous exact sequence match for {target_record.seqid}: {names}")
    return candidates[0] if candidates else (None, None)


def reverse_digest_index(records):
    """Index reverse-complement MD5 values without changing source records."""
    index = defaultdict(list)
    for record in records:
        index[sequence_digest(reverse_complement(record.sequence))].append(record)
    return index


class Record:
    def __init__(self, seqid, header, sequence):
        self.seqid = seqid
        self.header = header
        self.sequence = sequence

    @property
    def digest(self):
        return sequence_digest(self.sequence)


def die(message):
    raise RuntimeError(message)


def require(path, label):
    if not path.is_file() or path.stat().st_size == 0:
        die(f"Missing or empty {label}: {path}")
    return path


def fasta_records(path, uppercase=False):
    header = None
    chunks = []
    with path.open("r", encoding="ascii", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(chunks)
                    yield Record(header.split(None, 1)[0], header,
                                 sequence.upper() if uppercase else sequence)
                header = line[1:]
                chunks = []
            else:
                if header is None:
                    die(f"Sequence data before FASTA header in {path}")
                chunks.append(line)
    if header is not None:
        sequence = "".join(chunks)
        yield Record(header.split(None, 1)[0], header,
                     sequence.upper() if uppercase else sequence)


def aliases(value):
    """Return conservative aliases for common Prunus FASTA/GFF identifiers."""
    token = value.strip().split(None, 1)[0].lstrip(">")
    out = {token, token.split(".", 1)[0]}
    base = token.split(".", 1)[0]
    # EDTA on the PPA assembly drops the CAXKWO accession prefix and version.
    m = re.match(r"^[A-Za-z]+(\d{6,})$", base)
    if m:
        out.add(m.group(1))
    # Numeric chromosome aliases: chr01, chr1, 01, 1.
    m = re.fullmatch(r"(?:chr)?0*(\d+)", base, flags=re.I)
    if m:
        n = str(int(m.group(1)))
        out.update({n, f"{int(m.group(1)):02d}", f"chr{int(m.group(1)):02d}", f"chr{int(m.group(1))}"})
    # Old PSP names (C01_A1) and revised names (C1A_A1).
    m = re.fullmatch(r"C0*(\d+)_([AB])(\d+)", base)
    if m:
        n, sg, phase = m.groups()
        out.update({f"C{int(n)}_{sg}{phase}", f"C{int(n)}{sg}_{sg}{phase}",
                    f"C{int(n):02d}_{sg}{phase}", f"C{int(n):02d}{sg}_{sg}{phase}"})
    m = re.fullmatch(r"C0*(\d+)([AB])_([AB])(\d+)", base)
    if m:
        n, sg, sg2, phase = m.groups()
        out.update({f"C{int(n)}_{sg2}{phase}", f"C{int(n):02d}_{sg2}{phase}"})
    return {x for x in out if x}


def index_fasta(path):
    records = list(fasta_records(path, uppercase=True))
    by_digest = {}
    by_alias = {}
    for record in records:
        if record.digest in by_digest and by_digest[record.digest].seqid != record.seqid:
            die(f"Duplicate sequence content in {path}: {record.seqid} and {by_digest[record.digest].seqid}")
        by_digest[record.digest] = record
        for key in aliases(record.seqid) | aliases(record.header):
            previous = by_alias.get(key)
            if previous is not None and previous.seqid != record.seqid:
                # Ambiguous aliases are removed; MD5 matching remains authoritative.
                by_alias.pop(key, None)
            else:
                by_alias[key] = record
    return records, by_digest, by_alias


def write_fasta(handle, header, sequence):
    handle.write(f">{header}\n")
    for i in range(0, len(sequence), 80):
        handle.write(sequence[i:i + 80] + "\n")


def rewrite_intact_header(header, source_id, target_record, orientation):
    """Rewrite source chromosome/interval metadata in an intact-TE header."""
    source_pattern = re.escape(source_id)

    def replace_interval(match):
        left, right = int(match.group(1)), int(match.group(2))
        if orientation == "reverse_complement":
            left = len(target_record.sequence) - left + 1
            right = len(target_record.sequence) - right + 1
        return f"{target_record.seqid}:{left}..{right}"

    header = re.sub(
        rf"{source_pattern}:(\d+)\.\.(\d+)",
        replace_interval,
        header,
    )
    return header.replace(source_id, target_record.seqid)


def parse_attr(text):
    out = {}
    for field in text.rstrip(";\n").split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            out[key] = value
    return out


def classify(attrs):
    values = parse_attr(attrs)
    value = (values.get("classification") or values.get("Classification") or
             values.get("Class") or values.get("class") or "")
    if not value and "#" in values.get("Name", ""):
        value = values["Name"].rsplit("#", 1)[1]
    value = value.replace(" ", "")
    parts = [p for p in re.split(r"[/|]", value) if p]
    parent = parts[0] if parts else "unknown"
    child = parts[1] if len(parts) > 1 else "unknown"
    p = parent.lower()
    if p == "ltr": return "LTR", child
    if p == "line": return "LINE", child
    if p == "sine": return "SINE", child
    if p in {"dna", "tir", "mite"}: return "TIR", child
    if p in {"rc", "helitron"}: return "nonTIR", "helitron" if p == "helitron" else child
    if p == "nonltr": return "nonLTR", child
    if p in {"rdna", "rrna"}: return "rDNA", child
    return parent, child


class Summary:
    def __init__(self):
        self.count = Counter()
        self.bp = Counter()
        self.intervals = defaultdict(list)
        self.total = defaultdict(list)

    def add(self, seqid, start, end, category, subclass):
        self.count[(category, subclass)] += 1
        self.intervals[(seqid, category, subclass)].append((start, end))
        self.total[seqid].append((start, end))

    def register_sequence(self, seqid):
        """Keep chromosomes with no annotated TE in the summary denominator."""
        self.total.setdefault(seqid, [])

    @staticmethod
    def union(intervals):
        if not intervals:
            return 0
        intervals = sorted(intervals)
        total = 0
        start, end = intervals[0]
        for left, right in intervals[1:]:
            if left <= end + 1:
                end = max(end, right)
            else:
                total += end - start + 1
                start, end = left, right
        return total + end - start + 1

    def finish(self):
        for key, intervals in self.intervals.items():
            self.bp[(key[1], key[2])] += self.union(intervals)
        self.total_bp = sum(self.union(x) for x in self.total.values())
        self.total_count = sum(self.count.values())

    def write(self, path, genome_length):
        self.finish()
        rows = sorted(self.count, key=lambda x: (x[0].lower(), x[1].lower()))
        with path.open("w", encoding="utf-8", newline="\n") as out:
            out.write("Repeat Classes\n==============\n")
            out.write(f"Total Sequences: {len(self.total)}\nTotal Length: {genome_length} bp\n")
            out.write("Class                  Count        bpMasked    %masked\n")
            out.write("=====                  =====        ========     =======\n")
            for category, subclass in rows:
                bp = self.bp[(category, subclass)]
                out.write(f"    {category}/{subclass:<15} {self.count[(category, subclass)]:<12} {bp:<12} {100 * bp / genome_length:.2f}%\n")
            out.write("                      ---------------------------------\n")
            out.write(f"    total interspersed {self.total_count:<12} {self.total_bp:<12} {100 * self.total_bp / genome_length:.2f}%\n\n")
            out.write("Extraction note\n===============\n")
            out.write("Recalculated from the chromosome-subset EDTA TEanno.gff3 after revised haplotype assignment; no de novo EDTA run was performed.\n")


def read_map(path, species):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        die(f"Empty repartition map: {path}")
    result = {hap: {} for hap in HAPS[species]}
    assigned_sources = {}
    for row in rows:
        if row.get("species") != species:
            continue
        hap = row.get("haplotype", "")
        if hap not in result:
            continue
        target = row.get("new_seqid", "").strip()
        source = (row.get("original_chr32") if species == "pzi" else
                  (row.get("source_seqid") or row.get("old_seqid"))) or ""
        source = source.strip()
        if not target or not source:
            die(f"Invalid map row in {path}: {row}")
        if target in result[hap]:
            die(f"Duplicate target {hap}/{target} in {path}")
        previous = assigned_sources.get(source)
        if previous is not None and previous != (hap, target):
            die(f"Source chromosome {source} is assigned twice: {previous[0]}/{previous[1]} and {hap}/{target}")
        result[hap][target] = source
        assigned_sources[source] = (hap, target)
    for hap, mapping in result.items():
        if len(mapping) != 8:
            die(f"Expected eight chromosomes for {species}/{hap}, found {len(mapping)}")
    return result


def output_dir(root, species, hap):
    """Return the stable haplotype EDTA input directory.

    EDTA is run once on each complete assembly.  The haplotype results are
    chromosome subsets of those frozen annotations, so they belong beside the
    other haplotype inputs (genome/gff/cds/pep), not under the independent-run
    ``repeat`` tree.
    """
    return root / species / "edta" / hap


def materialize(source, target, force):
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() or target.is_symlink():
        if not force:
            die(f"Output exists (use --force): {target}")
        target.unlink()
    try:
        os.link(source, target)
    except OSError:
        shutil.copy2(source, target)


def run_one(species, hap_sources, maps, input_root, out_root, force):
    manifest = []
    # The complete EDTA files are shared by all haplotypes of this species.
    # Indexing them once avoids rereading multi-gigabyte FASTA files four times.
    source_fa, source_masked, source_gff, source_intact_gff, source_intact_fa, source_lib = next(iter(hap_sources.values()))
    for path, label in ((source_fa, "source EDTA FASTA"), (source_masked, "source masked FASTA"),
                        (source_gff, "source TEanno GFF3"), (source_intact_gff, "source intact GFF3"),
                        (source_intact_fa, "source intact FASTA"), (source_lib, "source EDTA library")):
        require(path, f"{species} {label}")
    source_records, source_by_digest, source_by_alias = index_fasta(source_fa)
    source_by_reverse_digest = reverse_digest_index(source_records)
    masked_records = {r.seqid: r for r in fasta_records(source_masked, uppercase=False)}

    for hap in HAPS[species]:
        target_fa = input_root / species / "genome" / f"{species}_{hap}.fa"
        require(target_fa, f"target genome {target_fa}")
        target_records, target_by_digest, _ = index_fasta(target_fa)
        target_by_id = {r.seqid: r for r in target_records}
        expected = maps[hap]
        selected = {}
        selected_orientation = {}
        for target_id, map_source in expected.items():
            if target_id not in target_by_id:
                die(f"{species}/{hap}: target chromosome {target_id} not found in {target_fa}")
            target_record = target_by_id[target_id]
            source_record, orientation = oriented_match(
                target_record, source_records, source_by_digest, source_by_reverse_digest
            )
            if source_record is None:
                die(
                    f"{species}/{hap}/{target_id}: no exact forward or reverse-complement "
                    f"MD5 match in source EDTA FASTA"
                )
            # The map source is checked against the matched source sequence through aliases.
            mapped_record = next((source_by_alias.get(key) for key in aliases(map_source)
                                  if source_by_alias.get(key) is not None), None)
            if mapped_record is not None and mapped_record.seqid != source_record.seqid:
                print(
                    f"[WARN] {species}/{hap}/{target_id}: map source {map_source} "
                    f"does not match exact sequence source; using "
                    f"{source_record.seqid} ({orientation})", file=sys.stderr
                )
            if source_record.seqid in selected:
                previous_target = selected[source_record.seqid][0].seqid
                die(
                    f"{species}/{hap}: source EDTA chromosome {source_record.seqid} "
                    f"matches both {previous_target} and {target_id}"
                )
            selected[source_record.seqid] = (target_record, source_record)
            selected_orientation[source_record.seqid] = orientation
            print(
                f"[MATCH] {species}/{hap}/{target_id} <- {source_record.seqid} "
                f"({orientation})", file=sys.stderr
            )
        source_to_target = {}
        selected_by_alias = {}
        for source_id, (target_record, source_record) in selected.items():
            source_to_target[source_id] = target_record.seqid
            for key in aliases(source_id):
                selected_by_alias[key] = (target_record.seqid, source_record)

        out = output_dir(out_root, species, hap)
        out.mkdir(parents=True, exist_ok=True)
        prefix = f"{species}_{hap}.fa"
        genome_out = out / prefix
        materialize(target_fa, genome_out, force)
        masked_out = out / f"{prefix}.mod.MAKER.masked"
        gff_out = out / f"{prefix}.mod.EDTA.TEanno.gff3"
        sum_out = out / f"{prefix}.mod.EDTA.TEanno.sum"
        intact_gff_out = out / f"{prefix}.mod.EDTA.intact.gff3"
        intact_fa_out = out / f"{prefix}.mod.EDTA.intact.fa"
        lib_out = out / f"{prefix}.mod.EDTA.TElib.fa"
        map_out = out / "chromosome_map.tsv"
        if map_out.exists() and not force: die(f"Output exists (use --force): {map_out}")
        with map_out.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write("target_seqid\tmap_source_id\tedta_source_seqid\torientation\tlength\ttarget_md5\tsource_md5\n")
            for target_id, map_source in expected.items():
                target_record = next(x for x, _ in selected.values() if x.seqid == target_id)
                source_record = next(y for x, y in selected.values() if x.seqid == target_id)
                orientation = selected_orientation[source_record.seqid]
                handle.write(
                    f"{target_id}\t{map_source}\t{source_record.seqid}\t{orientation}\t"
                    f"{len(target_record.sequence)}\t{target_record.digest}\t{source_record.digest}\n"
                )

        with masked_out.open("w", encoding="ascii", newline="\n") as out_handle:
            seen = set()
            for source_id, (target_record, source_record) in selected.items():
                record = masked_records.get(source_id)
                if record is None:
                    # Resolve accession/version aliases used by PPA masked FASTA.
                    for key in aliases(source_id):
                        record = next((x for x in masked_records.values() if key in aliases(x.seqid)), None)
                        if record is not None:
                            break
                if record is None or len(record.sequence) != len(source_record.sequence):
                    die(f"Masked FASTA record missing or length mismatch for {species}/{hap}/{source_id}")
                sequence = (
                    reverse_complement(record.sequence)
                    if selected_orientation[source_id] == "reverse_complement"
                    else record.sequence
                )
                if len(sequence) != len(target_record.sequence):
                    die(f"Masked FASTA transformed length mismatch for {species}/{hap}/{source_id}")
                write_fasta(out_handle, f"{target_record.seqid} source_seqid={source_id} orientation={selected_orientation[source_id]}", sequence)
                seen.add(source_id)
            if seen != set(selected):
                die(f"Masked FASTA extraction incomplete for {species}/{hap}")

        summary = Summary()
        for target_record, _ in selected.values():
            summary.register_sequence(target_record.seqid)
        selected_feature_ids = set()
        feature_orientations = {}
        with gff_out.open("w", encoding="utf-8", newline="\n") as out_handle, source_gff.open("r", encoding="utf-8", errors="replace") as source_handle:
            out_handle.write("##gff-version 3\n")
            out_handle.write(f"# Extracted from existing {species} EDTA annotation; revised haplotype assignment only.\n")
            for raw in source_handle:
                if raw.startswith("##FASTA"):
                    break
                if raw.startswith("##sequence-region"):
                    fields = raw.rstrip().split()
                    if len(fields) >= 4:
                        hit = selected_by_alias.get(fields[1])
                        if hit:
                            out_handle.write(f"##sequence-region {hit[0]} {fields[2]} {fields[3]}\n")
                    continue
                if raw.startswith("#") or not raw.strip():
                    continue
                fields = raw.rstrip("\n").split("\t")
                if len(fields) != 9:
                    continue
                hit = selected_by_alias.get(fields[0])
                if not hit:
                    continue
                target_id = hit[0]
                try:
                    start, end = int(fields[3]), int(fields[4])
                except ValueError:
                    die(f"Invalid coordinates in {source_gff}: {raw.rstrip()}")
                if start < 1 or end < start or end > len(hit[1].sequence):
                    die(f"Out-of-range coordinates in {source_gff}: {raw.rstrip()}")
                orientation = selected_orientation[hit[1].seqid]
                if orientation == "reverse_complement":
                    start, end = len(hit[1].sequence) - end + 1, len(hit[1].sequence) - start + 1
                    if fields[6] == "+":
                        fields[6] = "-"
                    elif fields[6] == "-":
                        fields[6] = "+"
                attrs = parse_attr(fields[8])
                feature_ids = [x for x in (attrs.get("ID", "").split(",") + attrs.get("Parent", "").split(",")) if x]
                selected_feature_ids.update(feature_ids)
                feature_orientations.update({feature_id: orientation for feature_id in feature_ids})
                fields[3], fields[4] = str(start), str(end)
                category, subclass = classify(fields[8])
                summary.add(target_id, start, end, category, subclass)
                fields[0] = target_id
                out_handle.write("\t".join(fields) + "\n")
        summary.write(sum_out, sum(len(x.sequence) for x, _ in selected.values()))

        # Intact GFF3 is filtered by chromosome and rewritten in the same way.
        with intact_gff_out.open("w", encoding="utf-8", newline="\n") as out_handle, source_intact_gff.open("r", encoding="utf-8", errors="replace") as source_handle:
            out_handle.write("##gff-version 3\n")
            out_handle.write(f"# Extracted from existing {species} EDTA intact annotation; revised haplotype assignment and orientation liftover.\n")
            for raw in source_handle:
                if raw.startswith("##FASTA"):
                    break
                if raw.startswith("##sequence-region"):
                    fields = raw.rstrip().split()
                    if len(fields) >= 4 and fields[1] in selected_by_alias:
                        out_handle.write(f"##sequence-region {selected_by_alias[fields[1]][0]} {fields[2]} {fields[3]}\n")
                    continue
                if raw.startswith("#") or not raw.strip():
                    continue
                fields = raw.rstrip("\n").split("\t")
                if len(fields) != 9:
                    continue
                hit = selected_by_alias.get(fields[0])
                if not hit:
                    continue
                try:
                    start, end = int(fields[3]), int(fields[4])
                except ValueError:
                    die(f"Invalid coordinates in {source_intact_gff}: {raw.rstrip()}")
                if start < 1 or end < start or end > len(hit[1].sequence):
                    die(f"Out-of-range coordinates in {source_intact_gff}: {raw.rstrip()}")
                orientation = selected_orientation[hit[1].seqid]
                if orientation == "reverse_complement":
                    start, end = len(hit[1].sequence) - end + 1, len(hit[1].sequence) - start + 1
                    if fields[6] == "+":
                        fields[6] = "-"
                    elif fields[6] == "-":
                        fields[6] = "+"
                attrs = parse_attr(fields[8])
                feature_ids = [x for x in (attrs.get("ID", "").split(",") + attrs.get("Parent", "").split(",")) if x]
                selected_feature_ids.update(feature_ids)
                feature_orientations.update({feature_id: orientation for feature_id in feature_ids})
                fields[3], fields[4] = str(start), str(end)
                fields[0] = hit[0]
                out_handle.write("\t".join(fields) + "\n")

        with intact_fa_out.open("w", encoding="ascii", newline="\n") as out_handle:
            for record in fasta_records(source_intact_fa, uppercase=False):
                match = re.search(r"(?:^|\s)ID=([^\s]+)", record.header)
                if match and match.group(1) in selected_feature_ids:
                    header = record.header
                    for source_id, (target_record, _) in selected.items():
                        header = rewrite_intact_header(
                            header, source_id, target_record, selected_orientation[source_id]
                        )
                    # EDTA intact.fa is already strand-oriented (its sequence
                    # follows the source GFF strand).  A whole-chromosome
                    # reverse-complement flips both the genomic coordinates
                    # and GFF strand, so the element sequence itself remains
                    # unchanged; reverse-complementing it here would flip it
                    # twice.
                    write_fasta(out_handle, header, record.sequence)
        materialize(source_lib, lib_out, force)

        completion = out / "edta.completed.tsv"
        with completion.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write("species\thaplotype\tgenome\tte_gff3\tte_summary\tmasked_genome\tintact_gff3\tintact_fasta\tte_library\tchromosome_map\tstatus\n")
            handle.write(f"{species}\t{hap}\t{genome_out}\t{gff_out}\t{sum_out}\t{masked_out}\t{intact_gff_out}\t{intact_fa_out}\t{lib_out}\t{map_out}\trepartitioned_oriented_md5_verified\n")
        manifest.append({"species": species, "haplotype": hap, "target_genome": str(target_fa), "output_dir": str(out), "chromosomes": str(len(selected)), "te_features": str(summary.total_count), "te_union_bp": str(summary.total_bp), "status": "pass"})
        print(f"[OK] {species}/{hap}: {len(selected)} chromosomes, {summary.total_count} TE features -> {out}")
    return manifest


def run_psp_pooled(hap_sources, maps, input_root, out_root, force):
    """Repartition PSP chromosomes from a pooled index of the four old haps."""
    datasets = []
    for old_hap, paths in hap_sources.items():
        source_fa, source_masked, source_gff, source_intact_gff, source_intact_fa, source_lib = paths
        for path, label in ((source_fa, "source FASTA"), (source_masked, "masked FASTA"),
                            (source_gff, "TEanno GFF3"), (source_intact_gff, "intact GFF3"),
                            (source_intact_fa, "intact FASTA"), (source_lib, "TE library")):
            require(path, f"PSP {old_hap} {label}")
        records, by_digest, by_alias = index_fasta(source_fa)
        datasets.append({"old_hap": old_hap, "paths": paths, "records": records,
                         "by_digest": by_digest, "by_alias": by_alias,
                         "by_reverse_digest": reverse_digest_index(records)})

    manifest = []
    for hap in HAPS["psp"]:
        target_fa = input_root / "psp" / "genome" / f"psp_{hap}.fa"
        require(target_fa, f"target genome {target_fa}")
        target_records, _, _ = index_fasta(target_fa)
        target_by_id = {r.seqid: r for r in target_records}
        selected = {}
        selected_orientation = {}
        for target_id, map_source in maps[hap].items():
            target_record = target_by_id.get(target_id)
            if target_record is None:
                die(f"PSP/{hap}: target chromosome {target_id} not found")
            pool_hits = []
            for dataset in datasets:
                source_record, orientation = oriented_match(
                    target_record,
                    dataset["records"],
                    dataset["by_digest"],
                    dataset["by_reverse_digest"],
                )
                if source_record is not None:
                    pool_hits.append((source_record, orientation, dataset["paths"]))
            if not pool_hits:
                die(
                    f"PSP/{hap}/{target_id}: no exact forward or reverse-complement "
                    f"MD5 match in the four old PSP FASTAs"
                )
            if len(pool_hits) > 1:
                names = ", ".join(f"{record.seqid} ({orientation})" for record, orientation, _ in pool_hits)
                die(f"Ambiguous PSP sequence match for {target_id}: {names}")
            source_record, orientation, paths = pool_hits[0]
            mapped = next((d["by_alias"].get(key) for d in datasets for key in aliases(map_source)
                           if d["by_alias"].get(key) is not None), None)
            if mapped is not None and mapped.seqid != source_record.seqid:
                print(
                    f"[WARN] psp/{hap}/{target_id}: map source {map_source} differs from "
                    f"exact source {source_record.seqid}; using exact sequence match",
                    file=sys.stderr,
                )
            if source_record.seqid in selected:
                previous_target = selected[source_record.seqid][0].seqid
                die(
                    f"PSP/{hap}: source EDTA chromosome {source_record.seqid} "
                    f"matches both {previous_target} and {target_id}"
                )
            selected[source_record.seqid] = (target_record, source_record, paths)
            selected_orientation[source_record.seqid] = orientation

        out = output_dir(out_root, "psp", hap); out.mkdir(parents=True, exist_ok=True)
        prefix = f"psp_{hap}.fa"
        genome_out = out / prefix; materialize(target_fa, genome_out, force)
        masked_out = out / f"{prefix}.mod.MAKER.masked"
        gff_out = out / f"{prefix}.mod.EDTA.TEanno.gff3"
        sum_out = out / f"{prefix}.mod.EDTA.TEanno.sum"
        intact_gff_out = out / f"{prefix}.mod.EDTA.intact.gff3"
        intact_fa_out = out / f"{prefix}.mod.EDTA.intact.fa"
        lib_out = out / f"{prefix}.mod.EDTA.TElib.fa"
        map_out = out / "chromosome_map.tsv"
        if map_out.exists() and not force: die(f"Output exists (use --force): {map_out}")
        with map_out.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write("target_seqid\tmap_source_id\tedta_source_seqid\torientation\tlength\ttarget_md5\tsource_md5\n")
            for target_id, map_source in maps[hap].items():
                target_record, source_record, _ = next(v for v in selected.values() if v[0].seqid == target_id)
                orientation = selected_orientation[source_record.seqid]
                handle.write(
                    f"{target_id}\t{map_source}\t{source_record.seqid}\t{orientation}\t"
                    f"{len(target_record.sequence)}\t{target_record.digest}\t{source_record.digest}\n"
                )
        selected_by_alias = {}
        for source_id, (target_record, source_record, paths) in selected.items():
            for key in aliases(source_id): selected_by_alias[key] = (target_record.seqid, source_record, paths)

        with masked_out.open("w", encoding="ascii", newline="\n") as out_handle:
            for dataset in datasets:
                selected_ids = {s for s, (_, _, p) in selected.items() if p == dataset["paths"]}
                if not selected_ids: continue
                for record in fasta_records(dataset["paths"][1], uppercase=False):
                    hit = next((selected_by_alias.get(k) for k in aliases(record.seqid) if selected_by_alias.get(k)), None)
                    if hit and hit[2] == dataset["paths"] and hit[1].seqid in selected_ids:
                        if len(record.sequence) != len(hit[1].sequence): die(f"PSP masked length mismatch for {record.seqid}")
                        sequence = (
                            reverse_complement(record.sequence)
                            if selected_orientation[hit[1].seqid] == "reverse_complement"
                            else record.sequence
                        )
                        write_fasta(
                            out_handle,
                            f"{hit[0]} source_seqid={record.seqid} orientation={selected_orientation[hit[1].seqid]}",
                            sequence,
                        )

        summary = Summary(); selected_feature_ids = set(); feature_orientations = {}
        for target_record, _, _ in selected.values():
            summary.register_sequence(target_record.seqid)
        with gff_out.open("w", encoding="utf-8", newline="\n") as out_handle:
            out_handle.write("##gff-version 3\n# Extracted from existing PSP EDTA annotations; revised haplotype assignment only.\n")
            for dataset in datasets:
                selected_ids = {s for s, (_, _, p) in selected.items() if p == dataset["paths"]}
                if not selected_ids: continue
                with dataset["paths"][2].open("r", encoding="utf-8", errors="replace") as source_handle:
                    for raw in source_handle:
                        if raw.startswith("##FASTA"): break
                        if raw.startswith("#") or not raw.strip(): continue
                        fields = raw.rstrip("\n").split("\t")
                        if len(fields) != 9: continue
                        hit = next((selected_by_alias.get(k) for k in aliases(fields[0]) if selected_by_alias.get(k)), None)
                        if not hit or hit[2] != dataset["paths"]: continue
                        try: start, end = int(fields[3]), int(fields[4])
                        except ValueError: die(f"Invalid PSP coordinates: {raw.rstrip()}")
                        if start < 1 or end < start or end > len(hit[1].sequence): die(f"Out-of-range PSP coordinates: {raw.rstrip()}")
                        orientation = selected_orientation[hit[1].seqid]
                        if orientation == "reverse_complement":
                            start, end = len(hit[1].sequence) - end + 1, len(hit[1].sequence) - start + 1
                            if fields[6] == "+": fields[6] = "-"
                            elif fields[6] == "-": fields[6] = "+"
                        attrs = parse_attr(fields[8])
                        feature_ids = [x for x in (attrs.get("ID", "").split(",") + attrs.get("Parent", "").split(",")) if x]
                        selected_feature_ids.update(feature_ids)
                        feature_orientations.update({feature_id: orientation for feature_id in feature_ids})
                        category, subclass = classify(fields[8]); summary.add(hit[0], start, end, category, subclass)
                        fields[0], fields[3], fields[4] = hit[0], str(start), str(end)
                        out_handle.write("\t".join(fields) + "\n")
        summary.write(sum_out, sum(len(x.sequence) for x, _, _ in selected.values()))

        with intact_gff_out.open("w", encoding="utf-8", newline="\n") as out_handle:
            out_handle.write("##gff-version 3\n# Extracted from existing PSP EDTA intact annotations; revised haplotype assignment and orientation liftover.\n")
            for dataset in datasets:
                if not any(p == dataset["paths"] for _, _, p in selected.values()): continue
                with dataset["paths"][3].open("r", encoding="utf-8", errors="replace") as source_handle:
                    for raw in source_handle:
                        if raw.startswith("##FASTA"): break
                        if raw.startswith("#") or not raw.strip(): continue
                        fields = raw.rstrip("\n").split("\t")
                        if len(fields) != 9: continue
                        hit = next((selected_by_alias.get(k) for k in aliases(fields[0]) if selected_by_alias.get(k)), None)
                        if not hit or hit[2] != dataset["paths"]: continue
                        try: start, end = int(fields[3]), int(fields[4])
                        except ValueError: die(f"Invalid PSP intact coordinates: {raw.rstrip()}")
                        if start < 1 or end < start or end > len(hit[1].sequence):
                            die(f"Out-of-range PSP intact coordinates: {raw.rstrip()}")
                        orientation = selected_orientation[hit[1].seqid]
                        if orientation == "reverse_complement":
                            start, end = len(hit[1].sequence) - end + 1, len(hit[1].sequence) - start + 1
                            if fields[6] == "+": fields[6] = "-"
                            elif fields[6] == "-": fields[6] = "+"
                        attrs = parse_attr(fields[8])
                        feature_ids = [x for x in (attrs.get("ID", "").split(",") + attrs.get("Parent", "").split(",")) if x]
                        selected_feature_ids.update(feature_ids)
                        feature_orientations.update({feature_id: orientation for feature_id in feature_ids})
                        fields[0], fields[3], fields[4] = hit[0], str(start), str(end)
                        out_handle.write("\t".join(fields) + "\n")
        with intact_fa_out.open("w", encoding="ascii", newline="\n") as out_handle:
            for dataset in datasets:
                if not any(p == dataset["paths"] for _, _, p in selected.values()): continue
                for record in fasta_records(dataset["paths"][4], uppercase=False):
                    match = re.search(r"(?:^|\s)ID=([^\s]+)", record.header)
                    if not match or match.group(1) not in selected_feature_ids: continue
                    header = record.header
                    for source_id, (target_record, _, paths) in selected.items():
                        if paths == dataset["paths"]:
                            header = rewrite_intact_header(
                                header, source_id, target_record, selected_orientation[source_id]
                            )
                    write_fasta(out_handle, header, record.sequence)
        # Merge the four old PSP libraries, removing exact duplicate sequences.
        seen_sequences = set()
        with lib_out.open("w", encoding="ascii", newline="\n") as out_handle:
            for dataset in datasets:
                for record in fasta_records(dataset["paths"][5], uppercase=False):
                    digest = hashlib.md5(record.sequence.upper().encode("ascii")).hexdigest()
                    if digest in seen_sequences: continue
                    seen_sequences.add(digest); write_fasta(out_handle, record.header, record.sequence)
        completion = out / "edta.completed.tsv"
        with completion.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write("species\thaplotype\tgenome\tte_gff3\tte_summary\tmasked_genome\tintact_gff3\tintact_fasta\tte_library\tchromosome_map\tstatus\n")
            handle.write(f"psp\t{hap}\t{genome_out}\t{gff_out}\t{sum_out}\t{masked_out}\t{intact_gff_out}\t{intact_fa_out}\t{lib_out}\t{map_out}\trepartitioned_oriented_md5_verified\n")
        manifest.append({"species": "psp", "haplotype": hap, "target_genome": str(target_fa), "output_dir": str(out), "chromosomes": str(len(selected)), "te_features": str(summary.total_count), "te_union_bp": str(summary.total_bp), "status": "pass"})
        print(f"[OK] psp/{hap}: {len(selected)} chromosomes, {summary.total_count} TE features -> {out}")
    return manifest


def build_sources(source_dir, species):
    def paths(stem):
        return (source_dir / f"{stem}.fa", source_dir / f"{stem}.fa.mod.MAKER.masked",
                source_dir / f"{stem}.fa.mod.EDTA.TEanno.gff3", source_dir / f"{stem}.fa.mod.EDTA.intact.gff3",
                source_dir / f"{stem}.fa.mod.EDTA.intact.fa", source_dir / f"{stem}.fa.mod.EDTA.TElib.fa")
    if species == "pzi":
        return {hap: paths("pzi_2n32.chromosomes") for hap in HAPS[species]}
    if species == "ppa":
        return {hap: paths("ppa") for hap in HAPS[species]}
    return {hap: paths("psp_2n32.chromosomes") for hap in HAPS[species]}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", nargs="+", choices=["pzi", "ppa", "psp", "all"], default=["all"])
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--inputs", type=Path, default=ROOT / "haploid_repartition" / "inputs")
    parser.add_argument(
        "--pzi-source-dir", type=Path,
        default=ROOT / "annotation" / "chr32_edta",
        help="complete PZI 2n=32 EDTA directory",
    )
    parser.add_argument(
        "--psp-source-dir", type=Path,
        default=ROOT / "data" / "psp" / "annotation" / "repeats" / "chr32_edta",
        help="complete PSP 2n=32 EDTA directory",
    )
    parser.add_argument(
        "--ppa-source-dir", type=Path,
        default=ROOT.parent / "prunus_old" / "analysis" / "prunus_pan_te" / "discovery" / "ppa",
        help="complete PPA 2n=16 EDTA directory",
    )
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if "all" in args.species:
        species_list = list(HAPS)
    else:
        species_list = list(dict.fromkeys(args.species))
    all_manifest = []
    for species in species_list:
        map_path = args.root / "haploid_repartition" / "maps" / f"{species}_repartition_map.tsv"
        require(map_path, f"{species} repartition map")
        maps = read_map(map_path, species)
        source_dir = {
            "pzi": args.pzi_source_dir,
            "psp": args.psp_source_dir,
            "ppa": args.ppa_source_dir,
        }[species]
        sources = build_sources(source_dir, species)
        all_manifest.extend(run_one(species, sources, maps, args.inputs, args.inputs, args.force))
    manifest_path = args.inputs / "edta_repartition_manifest.tsv"
    if manifest_path.exists() and not args.force:
        die(f"Output exists (use --force): {manifest_path}")
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(all_manifest[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(all_manifest)
    print(f"[DONE] Wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(2)
