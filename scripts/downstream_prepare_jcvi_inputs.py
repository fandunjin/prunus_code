#!/usr/bin/env python3
"""Prepare JCVI-compatible BED/PEP/CDS files from canonical inputs.

The canonical files in this project come from several annotation sources. This
script only normalizes IDs so that each BED feature and FASTA record share the
same accession, which is required by `jcvi.compara.catalog ortholog`.
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import unquote


TRANSCRIPT_TYPES = {"mrna", "transcript"}
CDS_TYPES = {"cds"}
ALIAS_KEYS = {
    "ID",
    "Name",
    "Alias",
    "Parent",
    "gene",
    "geneID",
    "locus_tag",
    "transcript_id",
    "protein_id",
    "orig_transcript_id",
    "orig_protein_id",
}


@dataclass
class CoordRecord:
    accn: str
    seqid: str
    start0: int
    end: int
    strand: str
    feature: str
    aliases: set[str] = field(default_factory=set)


@dataclass
class SeqRecord:
    source_id: str
    header: str
    seq: str
    aliases: set[str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", required=True, help="Project analysis directory")
    parser.add_argument("--outdir", required=True, help="Directory for JCVI input files")
    parser.add_argument("--summary", required=True, help="Output preparation summary TSV")
    parser.add_argument("--labels", nargs="+", required=True, help="Canonical species labels")
    return parser.parse_args()


def natural_key(value: str) -> list[object]:
    parts = re.split(r"(\d+)", value)
    return [int(x) if x.isdigit() else x for x in parts]


def clean_alias(value: str | None) -> str | None:
    if value is None:
        return None
    value = unquote(value.strip())
    if not value or value == ".":
        return None
    return value


def add_alias(aliases: set[str], value: str | None) -> None:
    value = clean_alias(value)
    if not value:
        return
    aliases.add(value)
    for token in re.split(r"[,|]", value):
        token = clean_alias(token)
        if token:
            aliases.add(token)
            for prefix in ("rna-", "gene-", "cds-"):
                if token.startswith(prefix) and len(token) > len(prefix):
                    aliases.add(token[len(prefix) :])
    for prefix in ("rna-", "gene-", "cds-"):
        if value.startswith(prefix) and len(value) > len(prefix):
            aliases.add(value[len(prefix) :])


def parse_gff_attributes(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in text.strip().split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
        else:
            continue
        attrs[key.strip()] = unquote(value.strip())
    return attrs


def aliases_from_attrs(attrs: dict[str, str]) -> set[str]:
    aliases: set[str] = set()
    for key in ALIAS_KEYS:
        value = attrs.get(key)
        if not value:
            continue
        for token in value.split(","):
            add_alias(aliases, token)
    dbxref = attrs.get("Dbxref")
    if dbxref:
        for token in dbxref.split(","):
            add_alias(aliases, token)
            if ":" in token:
                add_alias(aliases, token.split(":", 1)[1])
    return aliases


def parse_gff(gff: Path) -> tuple[list[CoordRecord], dict[str, int]]:
    transcript_records: dict[str, CoordRecord] = {}
    gene_records: dict[str, CoordRecord] = {}
    cds_groups: dict[str, CoordRecord] = {}
    cds_aliases: dict[str, set[str]] = defaultdict(set)
    stats = defaultdict(int)

    with open_text_auto(gff) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                stats["malformed_gff_lines"] += 1
                continue
            seqid, _source, feature, start, end, _score, strand, _phase, attr_text = fields
            feature_lc = feature.lower()
            attrs = parse_gff_attributes(attr_text)
            feature_id = clean_alias(attrs.get("ID"))
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                stats["bad_coordinate_lines"] += 1
                continue
            coord_start = min(start_i, end_i)
            coord_end = max(start_i, end_i)
            if coord_start != start_i or coord_end != end_i:
                stats["reversed_coordinate_lines"] += 1
            start0 = max(0, coord_start - 1)
            end_i = coord_end
            aliases = aliases_from_attrs(attrs)

            if feature_lc in TRANSCRIPT_TYPES:
                if not feature_id:
                    stats["transcript_without_id"] += 1
                    continue
                parent = clean_alias(attrs.get("Parent"))
                rec = CoordRecord(feature_id, seqid, start0, end_i, strand, feature, aliases)
                add_alias(rec.aliases, feature_id)
                if parent:
                    add_alias(rec.aliases, parent)
                    add_alias(rec.aliases, parent + "|" + feature_id)
                transcript_records[feature_id] = rec
                stats["transcript_features"] += 1
                continue

            if feature_lc == "gene":
                if feature_id:
                    rec = CoordRecord(feature_id, seqid, start0, end_i, strand, feature, aliases)
                    add_alias(rec.aliases, feature_id)
                    gene_records[feature_id] = rec
                    stats["gene_features"] += 1
                continue

            if feature_lc in CDS_TYPES:
                parents = attrs.get("Parent", "")
                if not parents and feature_id:
                    parents = feature_id
                for parent in [clean_alias(x) for x in parents.split(",") if clean_alias(x)]:
                    if parent not in cds_groups:
                        rec = CoordRecord(parent, seqid, start0, end_i, strand, "CDS_group", set())
                        add_alias(rec.aliases, parent)
                        add_alias(rec.aliases, parent + "|" + parent)
                        cds_groups[parent] = rec
                    rec = cds_groups[parent]
                    rec.start0 = min(rec.start0, start0)
                    rec.end = max(rec.end, end_i)
                    rec.aliases.update(aliases)
                    add_alias(rec.aliases, parent)
                    cds_aliases[parent].update(aliases)
                    stats["cds_features"] += 1

    for parent, aliases in cds_aliases.items():
        if parent in transcript_records:
            transcript_records[parent].aliases.update(aliases)
        elif parent in gene_records:
            gene_records[parent].aliases.update(aliases)

    records: dict[str, CoordRecord] = {}
    records.update(transcript_records)
    for parent, rec in cds_groups.items():
        if parent not in records:
            records[parent] = rec
    for gene_id, rec in gene_records.items():
        if gene_id not in records:
            records[gene_id] = rec

    ordered = sorted(records.values(), key=lambda r: (natural_key(r.seqid), r.start0, r.end, r.accn))
    stats["coordinate_records"] = len(ordered)
    return ordered, dict(stats)


def fasta_aliases(header: str) -> set[str]:
    aliases: set[str] = set()
    first = header.split()[0]
    add_alias(aliases, first)
    for token in first.split("|"):
        add_alias(aliases, token)
    if first.endswith(".p"):
        add_alias(aliases, first[:-2])
    for key, value in re.findall(r"([A-Za-z0-9_.:-]+)=([^\s\];]+)", header):
        add_alias(aliases, value)
        if key in {"ID", "transcript", "polypeptide", "locus", "gene", "Parent", "protein_id"}:
            add_alias(aliases, value)
    for value in re.findall(r"\[([A-Za-z0-9_.:-]+)=([^\]]+)\]", header):
        add_alias(aliases, value[1])
    return aliases


def open_text_auto(path: Path):
    with path.open("rb") as probe:
        is_gzip = probe.read(2) == b"\x1f\x8b"
    if is_gzip:
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def parse_fasta(path: Path, protein: bool = False) -> list[SeqRecord]:
    records: list[SeqRecord] = []
    header: str | None = None
    chunks: list[str] = []
    with open_text_auto(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append(make_seq_record(header, chunks, protein))
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        records.append(make_seq_record(header, chunks, protein))
    return records


def make_seq_record(header: str, chunks: list[str], protein: bool) -> SeqRecord:
    seq = "".join(chunks).replace(" ", "").replace("\t", "").upper()
    if protein:
        seq = seq.replace("*", "")
    source_id = header.split()[0]
    return SeqRecord(source_id=source_id, header=header, seq=seq, aliases=fasta_aliases(header))


def build_alias_index(records: list[CoordRecord]) -> dict[str, set[str]]:
    alias_to_ids: dict[str, set[str]] = defaultdict(set)
    for rec in records:
        rec.aliases.add(rec.accn)
        for alias in rec.aliases:
            alias_to_ids[alias].add(rec.accn)
    return alias_to_ids


def feature_rank(feature: str) -> int:
    feature_lc = feature.lower()
    if feature_lc in TRANSCRIPT_TYPES:
        return 0
    if feature_lc == "cds_group":
        return 1
    if feature_lc == "gene":
        return 2
    return 3


def match_sequences(
    seq_records: list[SeqRecord],
    alias_index: dict[str, set[str]],
    record_by_id: dict[str, CoordRecord],
) -> tuple[dict[str, SeqRecord], dict[str, str], int, int]:
    matched: dict[str, SeqRecord] = {}
    matched_alias: dict[str, str] = {}
    ambiguous = 0
    unmatched = 0
    for seq in seq_records:
        candidates: list[tuple[tuple[int, int, int, int, str], str, str]] = []
        for alias in seq.aliases:
            if alias not in alias_index:
                continue
            ids = alias_index[alias]
            if len(ids) > 20:
                ambiguous += 1
                continue
            for candidate in ids:
                rec = record_by_id[candidate]
                exact = 0 if alias == candidate else 1
                # Prefer aliases that identify a single feature (for example an
                # NCBI protein_id) over gene-level aliases shared by isoforms.
                rank = (exact, len(ids), feature_rank(rec.feature), -len(alias), candidate)
                candidates.append((rank, candidate, alias))
        if not candidates:
            unmatched += 1
            continue
        candidates.sort(key=lambda x: x[0])
        hit = candidates[0][1]
        hit_alias = candidates[0][2]
        if hit not in matched or len(seq.seq) > len(matched[hit].seq):
            matched[hit] = seq
            matched_alias[hit] = hit_alias
    return matched, matched_alias, ambiguous, unmatched


def wrap(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def write_species(label: str, analysis_dir: Path, outdir: Path) -> dict[str, object]:
    canonical = analysis_dir / "canonical" / label
    gff = canonical / "primary.gff3"
    pep = canonical / "primary.pep.fa"
    cds = canonical / "primary.cds.fa"
    for path in (gff, pep):
        if not path.is_file():
            raise FileNotFoundError(f"{label}: missing required file {path}")

    species_out = outdir / label
    species_out.mkdir(parents=True, exist_ok=True)
    records, gff_stats = parse_gff(gff)
    alias_index = build_alias_index(records)
    record_by_id = {rec.accn: rec for rec in records}
    pep_records = parse_fasta(pep, protein=True)
    pep_by_id, pep_alias, pep_ambiguous, pep_unmatched = match_sequences(
        pep_records, alias_index, record_by_id
    )
    cds_by_id: dict[str, SeqRecord] = {}
    cds_alias: dict[str, str] = {}
    cds_ambiguous = 0
    cds_unmatched = 0
    if cds.is_file():
        cds_records = parse_fasta(cds, protein=False)
        cds_by_id, cds_alias, cds_ambiguous, cds_unmatched = match_sequences(
            cds_records, alias_index, record_by_id
        )

    usable = [rec for rec in records if rec.accn in pep_by_id]
    bed_out = species_out / f"{label}.bed"
    pep_out = species_out / f"{label}.pep"
    cds_out = species_out / f"{label}.cds"
    map_out = species_out / f"{label}.id_map.tsv"
    seqids_out = species_out / f"{label}.seqids"

    with bed_out.open("w", encoding="utf-8") as bed_fh, pep_out.open(
        "w", encoding="utf-8"
    ) as pep_fh, cds_out.open("w", encoding="utf-8") as cds_fh, map_out.open(
        "w", encoding="utf-8"
    ) as map_fh:
        print(
            "\t".join(
                [
                    "output_id",
                    "gff_feature",
                    "seqid",
                    "start0",
                    "end",
                    "strand",
                    "pep_source_id",
                    "pep_match_alias",
                    "cds_source_id",
                    "cds_match_alias",
                ]
            ),
            file=map_fh,
        )
        for rec in usable:
            print(
                "\t".join([rec.seqid, str(rec.start0), str(rec.end), rec.accn, "0", rec.strand]),
                file=bed_fh,
            )
            pep_rec = pep_by_id[rec.accn]
            print(f">{rec.accn}\n{wrap(pep_rec.seq)}", file=pep_fh)
            cds_rec = cds_by_id.get(rec.accn)
            if cds_rec is not None:
                print(f">{rec.accn}\n{wrap(cds_rec.seq)}", file=cds_fh)
            print(
                "\t".join(
                    [
                        rec.accn,
                        rec.feature,
                        rec.seqid,
                        str(rec.start0),
                        str(rec.end),
                        rec.strand,
                        pep_rec.source_id,
                        pep_alias.get(rec.accn, ""),
                        cds_rec.source_id if cds_rec else "",
                        cds_alias.get(rec.accn, ""),
                    ]
                ),
                file=map_fh,
            )

    seqids = []
    seen_seqids = set()
    for rec in usable:
        if rec.seqid not in seen_seqids:
            seqids.append(rec.seqid)
            seen_seqids.add(rec.seqid)
    seqids_out.write_text(",".join(seqids) + "\n", encoding="utf-8")

    return {
        "label": label,
        "gff_records": len(records),
        "pep_records": len(pep_records),
        "usable_bed_pep_records": len(usable),
        "cds_matched_records": sum(1 for rec in usable if rec.accn in cds_by_id),
        "pep_unmatched_records": pep_unmatched,
        "pep_ambiguous_alias_records": pep_ambiguous,
        "cds_unmatched_records": cds_unmatched,
        "cds_ambiguous_alias_records": cds_ambiguous,
        "transcript_features": gff_stats.get("transcript_features", 0),
        "gene_features": gff_stats.get("gene_features", 0),
        "cds_features": gff_stats.get("cds_features", 0),
        "bed": bed_out,
        "pep": pep_out,
        "cds": cds_out,
    }


def main() -> int:
    args = parse_args()
    analysis_dir = Path(args.analysis_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    summaries = []
    for label in args.labels:
        print(f"[INFO] Preparing JCVI inputs for {label}", file=sys.stderr)
        summaries.append(write_species(label, analysis_dir, outdir))

    summary_path = Path(args.summary)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "label",
        "gff_records",
        "pep_records",
        "usable_bed_pep_records",
        "cds_matched_records",
        "pep_unmatched_records",
        "pep_ambiguous_alias_records",
        "cds_unmatched_records",
        "cds_ambiguous_alias_records",
        "transcript_features",
        "gene_features",
        "cds_features",
        "bed",
        "pep",
        "cds",
    ]
    with summary_path.open("w", encoding="utf-8") as out:
        print("\t".join(fields), file=out)
        for row in summaries:
            print("\t".join(str(row.get(field, "")) for field in fields), file=out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
