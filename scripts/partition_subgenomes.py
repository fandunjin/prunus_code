#!/usr/bin/env python3
"""Create reviewed PPA_A/PPA_B genome and annotation packages.

This script deliberately does not infer subgenomes.  It only applies a
manually approved correspondence table produced from the PPA--PZI macro-
synteny review.  Inputs and validation are kept strict so a partial or
ambiguous table cannot silently produce downstream phylogenetic inputs.
"""

import argparse
import csv
import gzip
import hashlib
import re
from collections import Counter
from pathlib import Path


REQUIRED_COLUMNS = {
    "ppa_source_seqid",
    "pzi_copy1_chr",
    "subgenome",
    "renamed_seqid",
    "ppa_union_covered_bp",
    "best_to_second_ratio",
    "review_status",
}
RENAME_RE = re.compile(r"^C(0[1-8])_([AB])$")


def open_text(path: Path):
    return (
        gzip.open(path, "rt", encoding="utf-8", errors="replace")
        if path.suffix == ".gz"
        else path.open("r", encoding="utf-8", errors="replace")
    )


def fasta_records(path: Path):
    name = None
    header = ""
    chunks = []
    with open_text(path) as handle:
        for raw in handle:
            if raw.startswith(">"):
                if name is not None:
                    yield name, header, "".join(chunks)
                header = raw.rstrip("\r\n")[1:]
                name = header.split(None, 1)[0]
                chunks = []
            else:
                chunks.append(raw.strip())
    if name is not None:
        yield name, header, "".join(chunks)


def parse_attrs(text):
    attrs = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
    return attrs


def md5sum(path):
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def tsv(path, rows, names):
    with path.open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(
            out, fieldnames=names, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def read_primary_map(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        needed = {"ppa_source_seqid", "ppa_neutral_seqid"}
        if not reader.fieldnames or not needed <= set(reader.fieldnames):
            raise SystemExit(
                f"ERROR: primary chromosome map lacks {sorted(needed)}: {path}"
            )
        rows = list(reader)
    sources = {row["ppa_source_seqid"] for row in rows}
    if len(rows) != 16 or len(sources) != 16:
        raise SystemExit(
            f"ERROR: primary chromosome map must contain 16 unique rows: {path}"
        )
    return sources


def read_assignment(path, primary_sources):
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not REQUIRED_COLUMNS <= set(reader.fieldnames):
            raise SystemExit(
                f"ERROR: assignment table lacks required columns {sorted(REQUIRED_COLUMNS)}: {path}"
            )
        rows = list(reader)
    if len(rows) != 16:
        raise SystemExit(
            f"ERROR: assignment table must contain exactly 16 rows, found {len(rows)}"
        )
    assignments = {}
    for line_no, row in enumerate(rows, 2):
        source = row["ppa_source_seqid"].strip()
        status = row["review_status"].strip().lower()
        sg = row["subgenome"].strip().upper()
        renamed = row["renamed_seqid"].strip()
        match = RENAME_RE.fullmatch(renamed)
        if not source or source in assignments:
            raise SystemExit(
                f"ERROR: duplicate/empty ppa_source_seqid at assignment row {line_no}"
            )
        if source not in primary_sources:
            raise SystemExit(
                f"ERROR: assignment source is not a PPA primary chromosome: {source}"
            )
        if status != "approved":
            raise SystemExit(
                f"ERROR: review_status must be approved at assignment row {line_no}: {source}={row['review_status']!r}"
            )
        if sg not in {"A", "B"} or match is None or match.group(2) != sg:
            raise SystemExit(
                f"ERROR: inconsistent subgenome/renamed_seqid at assignment row {line_no}: {source}"
            )
        assignments[source] = {key: value.strip() for key, value in row.items()}
    if set(assignments) != primary_sources:
        missing = sorted(primary_sources - set(assignments))
        extra = sorted(set(assignments) - primary_sources)
        raise SystemExit(
            f"ERROR: assignment does not cover exactly the primary chromosomes; missing={missing}, extra={extra}"
        )
    sub_counts = Counter(row["subgenome"].upper() for row in assignments.values())
    if sub_counts != Counter({"A": 8, "B": 8}):
        raise SystemExit(
            f"ERROR: expected 8 A and 8 B chromosomes, observed {dict(sub_counts)}"
        )
    for sg in ("A", "B"):
        labels = [
            row["renamed_seqid"]
            for row in assignments.values()
            if row["subgenome"].upper() == sg
        ]
        expected = {f"C{i:02d}_{sg}" for i in range(1, 9)}
        if set(labels) != expected or len(labels) != len(set(labels)):
            raise SystemExit(
                f"ERROR: {sg} must contain C01--C08 exactly once; observed {sorted(labels)}"
            )
    return assignments


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome", type=Path, required=True)
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--cds", type=Path, required=True)
    parser.add_argument("--pep", type=Path, required=True)
    parser.add_argument("--primary-map", type=Path, required=True)
    parser.add_argument("--assignment", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    for label, path in (
        ("genome", args.genome),
        ("GFF3", args.gff),
        ("CDS", args.cds),
        ("protein", args.pep),
        ("primary map", args.primary_map),
        ("assignment", args.assignment),
    ):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"ERROR: missing {label}: {path}")
    if args.out_dir.exists() and any(args.out_dir.iterdir()) and not args.force:
        raise SystemExit(
            f"ERROR: output directory exists and is nonempty: {args.out_dir}; review it or rerun with --force"
        )
    primary_sources = read_primary_map(args.primary_map)
    assignments = read_assignment(args.assignment, primary_sources)
    aliases = {source: row["renamed_seqid"] for source, row in assignments.items()}
    args.out_dir.mkdir(parents=True, exist_ok=True)
    dirs = {kind: args.out_dir / kind for kind in ("genome", "gff", "cds", "pep")}
    for directory in dirs.values():
        directory.mkdir(parents=True, exist_ok=True)

    genome_out = {sg: dirs["genome"] / f"ppa_hap{sg}.fa" for sg in ("A", "B")}
    genome_handles = {
        sg: path.open("w", encoding="utf-8", newline="\n")
        for sg, path in genome_out.items()
    }
    genome_counts = Counter()
    genome_bp = Counter()
    try:
        for source, header, sequence in fasta_records(args.genome):
            if source not in aliases:
                continue
            renamed = aliases[source]
            sg = renamed[-1]
            genome_handles[sg].write(
                f">{renamed} source_seqid={source} source_header={header}\n"
            )
            for start in range(0, len(sequence), 80):
                genome_handles[sg].write(sequence[start : start + 80] + "\n")
            genome_counts[sg] += 1
            genome_bp[sg] += len(sequence)
    finally:
        for handle in genome_handles.values():
            handle.close()
    if dict(genome_counts) != {"A": 8, "B": 8}:
        raise SystemExit(
            f"ERROR: genome extraction failed, records={dict(genome_counts)}"
        )

    gff_out = {sg: dirs["gff"] / f"ppa_hap{sg}.gff3" for sg in ("A", "B")}
    gff_handles = {
        sg: path.open("w", encoding="utf-8", newline="\n")
        for sg, path in gff_out.items()
    }
    gff_features = Counter()
    mrnas = {"A": set(), "B": set()}
    mrna_aliases = {}
    excluded_features = 0
    try:
        for raw in open_text(args.gff):
            line = raw.rstrip("\r\n")
            if line.startswith("##sequence-region"):
                fields = line.split()
                if len(fields) >= 2 and fields[1] in aliases:
                    renamed = aliases[fields[1]]
                    fields[1] = renamed
                    gff_handles[renamed[-1]].write(" ".join(fields) + "\n")
                continue
            if line.startswith("#") or not line.strip():
                for handle in gff_handles.values():
                    handle.write(line + "\n")
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise SystemExit(f"ERROR: malformed GFF3 row: {line[:160]}")
            source = fields[0]
            if source not in aliases:
                excluded_features += 1
                continue
            renamed = aliases[source]
            sg = renamed[-1]
            fields[0] = renamed
            gff_handles[sg].write("\t".join(fields) + "\n")
            gff_features[sg] += 1
            if fields[2] in {"mRNA", "transcript"}:
                ident = parse_attrs(fields[8]).get("ID", "")
                if not ident:
                    raise SystemExit(f"ERROR: {fields[2]} lacks ID: {line[:160]}")
                if ident in mrna_aliases:
                    raise SystemExit(
                        f"ERROR: duplicate mRNA ID in selected chromosomes: {ident}"
                    )
                mrnas[sg].add(ident)
                mrna_aliases[ident] = (sg, source, renamed)
    finally:
        for handle in gff_handles.values():
            handle.close()

    fasta_outputs = {
        "cds": (args.cds, dirs["cds"], "cds"),
        "pep": (args.pep, dirs["pep"], "pep"),
    }
    fasta_counts = {"cds": Counter(), "pep": Counter()}
    fasta_seqid_rewrites = {"cds": Counter(), "pep": Counter()}
    for kind, (source_path, out_path, suffix) in fasta_outputs.items():
        handles = {
            sg: (out_path / f"ppa_hap{sg}.{suffix}.fa").open(
                "w", encoding="utf-8", newline="\n"
            )
            for sg in ("A", "B")
        }
        try:
            for record_id, header, sequence in fasta_records(source_path):
                assignment = mrna_aliases.get(record_id)
                if assignment is None:
                    continue
                sg, source_seqid, renamed_seqid = assignment
                old_token = "seqid=" + source_seqid
                new_token = "seqid=" + renamed_seqid
                if old_token not in header:
                    raise SystemExit(
                        f"ERROR: selected {kind} FASTA header lacks expected {old_token}: {record_id}"
                    )
                new_header = header.replace(old_token, new_token)
                if old_token in new_header:
                    raise SystemExit(
                        f"ERROR: failed to rewrite {kind} seqid in header: {record_id}"
                    )
                handles[sg].write(f">{new_header}\n")
                for start in range(0, len(sequence), 80):
                    handles[sg].write(sequence[start : start + 80] + "\n")
                fasta_counts[kind][sg] += 1
                fasta_seqid_rewrites[kind][sg] += 1
        finally:
            for handle in handles.values():
                handle.close()

    map_rows = []
    for source in sorted(
        assignments, key=lambda key: assignments[key]["renamed_seqid"]
    ):
        row = dict(assignments[source])
        row["ppa_source_seqid"] = source
        map_rows.append(row)
    output_map = args.out_dir / "ppa_A_B.rename_map.tsv"
    columns = [
        "ppa_source_seqid",
        "pzi_copy1_chr",
        "subgenome",
        "renamed_seqid",
        "ppa_union_covered_bp",
        "best_to_second_ratio",
        "review_status",
    ]
    extra = sorted(set().union(*(row.keys() for row in map_rows)) - set(columns))
    tsv(output_map, map_rows, columns + extra)

    validation_rows = []
    for sg in ("A", "B"):
        validation_rows.append(
            {
                "check": f"hap{sg}_chromosomes",
                "observed": genome_counts[sg],
                "expected": 8,
                "status": "pass" if genome_counts[sg] == 8 else "fail",
            }
        )
        validation_rows.append(
            {
                "check": f"hap{sg}_GFF_features",
                "observed": gff_features[sg],
                "expected": ">0",
                "status": "pass" if gff_features[sg] else "fail",
            }
        )
        validation_rows.append(
            {
                "check": f"hap{sg}_mRNA",
                "observed": len(mrnas[sg]),
                "expected": ">0",
                "status": "pass" if mrnas[sg] else "fail",
            }
        )
        for kind in ("cds", "pep"):
            validation_rows.append(
                {
                    "check": f"hap{sg}_{kind}_records",
                    "observed": fasta_counts[kind][sg],
                    "expected": f"<= {len(mrnas[sg])}",
                    "status": (
                        "pass" if fasta_counts[kind][sg] <= len(mrnas[sg]) else "fail"
                    ),
                }
            )
            validation_rows.append(
                {
                    "check": f"hap{sg}_{kind}_seqid_rewrites",
                    "observed": fasta_seqid_rewrites[kind][sg],
                    "expected": fasta_counts[kind][sg],
                    "status": (
                        "pass"
                        if fasta_seqid_rewrites[kind][sg] == fasta_counts[kind][sg]
                        else "fail"
                    ),
                }
            )
    validation_rows.append(
        {
            "check": "excluded_unplaced_or_nonchromosomal_GFF_features",
            "observed": excluded_features,
            "expected": "reported",
            "status": "info",
        }
    )
    tsv(
        args.out_dir / "validation.tsv",
        validation_rows,
        ["check", "observed", "expected", "status"],
    )

    outputs = [
        *genome_out.values(),
        *gff_out.values(),
        *(dirs["cds"] / f"ppa_hap{sg}.cds.fa" for sg in ("A", "B")),
        *(dirs["pep"] / f"ppa_hap{sg}.pep.fa" for sg in ("A", "B")),
        output_map,
        args.out_dir / "validation.tsv",
    ]
    manifest = []
    sources = {"genome": args.genome, "gff": args.gff, "cds": args.cds, "pep": args.pep}
    for path in outputs:
        kind = next((key for key in sources if path.parent.name == key), "metadata")
        manifest.append(
            {
                "file": path.relative_to(args.out_dir).as_posix(),
                "bytes": path.stat().st_size,
                "md5": md5sum(path),
                "source": str(sources.get(kind, args.assignment)),
            }
        )
    tsv(args.out_dir / "manifest.tsv", manifest, ["file", "bytes", "md5", "source"])
    (args.out_dir / "README.txt").write_text(
        "PPA reviewed subgenome package\n\n"
        "PPA_A/PPA_B are structural correspondence labels derived from a manually approved\n"
        "PPA-primary-to-PZI-copy1 assignment table. They are not claims that PPA is a\n"
        "direct parental species of PZI. Unplaced PPA sequences and their annotations\n"
        "are deliberately excluded from these chromosome-level packages.\n",
        encoding="utf-8",
    )
    failures = [row for row in validation_rows if row["status"] == "fail"]
    if failures:
        raise SystemExit(f"ERROR: validation failed: {failures}")
    print(f"[OK] Wrote reviewed PPA subgenome package: {args.out_dir}")


if __name__ == "__main__":
    main()
