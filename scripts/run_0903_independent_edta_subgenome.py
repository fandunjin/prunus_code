#!/usr/bin/env python3
"""Summarize independent EDTA annotations without using pan-TE results.

The script has two deliberate boundaries:

* abundance is read only from the independent EDTA GFF3 files listed in the
  manifest; no file under analysis/prunus_pan_te is searched or accepted;
* PZI/PSP/PPA use the revised A/B chromosome packages under
  0903/haploid_repartition (per-haplotype genomes plus md5-verified
  repartitioned EDTA annotations); PPE/PCAM/GTR/RCH use 0903/TE.

Run ``--mode all`` to create the manifest, validate inputs and write summary
tables. Missing libraries/intact files are reported explicitly. They are not
silently substituted with another library.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


TARGET_ORDER = [
    "pziA1", "pziB1",
    "pspA1", "pspB1",
    "ppaA", "ppaB", "ppe", "pcam", "gtr", "rch",
]


def first_existing(*paths: Path) -> Path:
    """Return the first existing path, otherwise keep the first candidate."""
    for path in paths:
        if path.exists():
            return path
    return paths[0]


def fasta_records(path: Path):
    """Yield (identifier, sequence) without requiring an external FASTA tool."""
    name = None
    chunks: list[str] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def write_alias_fasta(source: Path, output: Path, aliases: dict[str, str]) -> None:
    """Write selected source records under the exact IDs used by an EDTA GFF3."""
    output.parent.mkdir(parents=True, exist_ok=True)
    found: set[str] = set()
    with output.open("w", encoding="utf-8") as handle:
        for source_id, sequence in fasta_records(source):
            alias = aliases.get(source_id)
            if not alias:
                continue
            found.add(alias)
            handle.write(f">{alias}\n")
            for i in range(0, len(sequence), 80):
                handle.write(sequence[i:i + 80] + "\n")
    missing = sorted(set(aliases.values()) - found)
    if missing:
        raise RuntimeError(f"Could not materialize FASTA aliases from {source}: {missing[:8]}")


def default_rows(project_root: Path) -> list[dict[str, str]]:
    """Build an explicit, auditable manifest for the independent EDTA data."""
    root = project_root.resolve()
    d = root / "data"
    r = root / "0903"
    rows: list[dict[str, str]] = []

    def add(
        target: str,
        species: str,
        subgenome: str,
        haplotype: str,
        genome: Path,
        annotation_genome: Path,
        gff: Path,
        library: Path,
        intact_gff: Path,
        intact_fasta: Path,
        summary: Path,
        seq_regex: str = "",
        source: str = "independent_EDTA",
        assignment_basis: str = "existing_haplotype_label",
        annotation_scope: str = "independent_haplotype",
        library_owner: str = "",
    ) -> None:
        rows.append({
            "target": target,
            "species": species,
            "subgenome": subgenome,
            "haplotype": haplotype,
            "genome_fasta": str(genome),
            "annotation_genome_fasta": str(annotation_genome),
            "edta_gff3": str(gff),
            "edta_library_fasta": str(library),
            "edta_intact_gff3": str(intact_gff),
            "edta_intact_fasta": str(intact_fasta),
            "edta_summary": str(summary),
            "seq_regex": seq_regex,
            "source": source,
            "assignment_basis": assignment_basis,
            "annotation_scope": annotation_scope,
            "library_owner": library_owner or target,
        })

    # PZI: complete independent per-haplotype EDTA runs are stored in 0903.
    # Panel restriction (2026-09): A1/B1 only for the cross-species TE re-analysis.
    # Both per-haplotype TElib.fa files are symlinks to the same combined 2n=32
    # library, so the tree pipeline must treat them as one library owner
    # (pzi_combined) to search both per-haplotype GFF3 files for copies; with
    # per-haplotype owners the library is deduplicated on load and the B1 GFF3
    # is never queried (observed 2026-09-16: zero pziB1 copies in all trees).
    for hap in ("A1", "B1"):
        edta = r / "haploid_repartition" / "inputs" / "pzi" / "edta" / f"hap{hap}"
        fa = edta / f"pzi_hap{hap}.fa"
        add(
            f"pzi{hap}", "PZI", hap[0], hap, fa, fa,
            edta / f"pzi_hap{hap}.fa.mod.EDTA.TEanno.gff3",
            edta / f"pzi_hap{hap}.fa.mod.EDTA.TElib.fa",
            edta / f"pzi_hap{hap}.fa.mod.EDTA.intact.gff3",
            edta / f"pzi_hap{hap}.fa.mod.EDTA.intact.fa",
            edta / f"pzi_hap{hap}.fa.mod.EDTA.TEanno.sum",
            library_owner="pzi_combined",
        )

    # PSP / PPA: revised A/B chromosome packages (0903/haploid_repartition,
    # "Revised A/B chromosome packages").  The EDTA directories under
    # inputs/{psp,ppa}/edta/hap* contain repartitions of the phased combined
    # runs with md5-verified chromosome maps (chromosome_map.tsv) and are in
    # the same coordinate system as the per-haplotype genome FASTAs.
    for target, species, sub, hap, subdir, stem, owner in (
        ("pspA1", "PSP", "A", "A1", "psp/edta/hapA1", "psp_hapA1", "psp_combined"),
        ("pspB1", "PSP", "B", "B1", "psp/edta/hapB1", "psp_hapB1", "psp_combined"),
        ("ppaA", "PPA", "A", "A", "ppa/edta/hapA", "ppa_hapA", "ppa_combined"),
        ("ppaB", "PPA", "B", "B", "ppa/edta/hapB", "ppa_hapB", "ppa_combined"),
    ):
        edta = r / "haploid_repartition" / "inputs" / subdir
        fa = edta / f"{stem}.fa"
        add(
            target, species, sub, hap, fa, fa,
            edta / f"{stem}.fa.mod.EDTA.TEanno.gff3",
            edta / f"{stem}.fa.mod.EDTA.TElib.fa",
            edta / f"{stem}.fa.mod.EDTA.intact.gff3",
            edta / f"{stem}.fa.mod.EDTA.intact.fa",
            edta / f"{stem}.fa.mod.EDTA.TEanno.sum",
            source="repartitioned_from_combined_md5_verified",
            assignment_basis="0903 revised A/B packages (haploid_repartition chromosome_map.tsv)",
            annotation_scope="independent_haplotype",
            library_owner=owner,
        )

    # PPE, PCAM, GTR and RCH: independent single-genome EDTA outputs under 0903/TE.
    # gtr's EDTA GFF uses accession-stripped sequence IDs (031080, 010000010), so
    # its copy extraction uses gtr.gffmatch.fa when that renamed FASTA is present.
    for target, species, folder, stem in (
        ("ppe", "PPE", "ppe", "ppe"),
        ("pcam", "PCAM", "pcam", "pcam"),
        ("gtr", "GTR", "gtr", "gtr"),
        ("rch", "RCH", "rch", "rch"),
    ):
        edta = r / "TE" / folder
        fa = edta / f"{stem}.fa"
        gffmatch = edta / f"{stem}.gffmatch.fa"
        annotation_fa = gffmatch if gffmatch.is_file() else fa
        add(
            target, species, "single", "single", fa, annotation_fa,
            edta / f"{stem}.fa.mod.EDTA.TEanno.gff3",
            edta / f"{stem}.fa.mod.EDTA.TElib.fa",
            edta / f"{stem}.fa.mod.EDTA.intact.gff3",
            edta / f"{stem}.fa.mod.EDTA.intact.fa",
            edta / f"{stem}.fa.mod.EDTA.TEanno.sum",
            assignment_basis="single genome; no A/B partition",
        )
    return rows


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_attrs(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for field in text.strip().split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            out[key] = value
    return out


def fasta_length(path: Path) -> int:
    length = 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                length += len(line.strip())
    return length


def union_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    total = 0
    start, end = sorted(intervals)[0]
    for nxt_start, nxt_end in sorted(intervals)[1:]:
        if nxt_start <= end + 1:
            end = max(end, nxt_end)
        else:
            total += end - start + 1
            start, end = nxt_start, nxt_end
    return total + end - start + 1


def classify_family(attrs: dict[str, str]) -> tuple[str, str, str]:
    family = attrs.get("Name") or attrs.get("Target") or attrs.get("ID") or "unknown"
    classification = attrs.get("classification", "unknown")
    primary = classification.split("/", 1)[0] if classification else "unknown"
    return family, classification, primary


def validate_manifest(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    reports: list[dict[str, object]] = []
    for row in rows:
        required = {
            "genome_fasta": Path(row["genome_fasta"]),
            "annotation_genome_fasta": Path(row["annotation_genome_fasta"]),
            "edta_gff3": Path(row["edta_gff3"]),
        }
        optional = {
            "edta_library_fasta": Path(row["edta_library_fasta"]),
            "edta_intact_gff3": Path(row["edta_intact_gff3"]),
            "edta_intact_fasta": Path(row["edta_intact_fasta"]),
            "edta_summary": Path(row["edta_summary"]),
        }
        status = "pass" if all(path.is_file() for path in required.values()) else "missing_required"
        reports.append({
            "target": row["target"],
            "species": row["species"],
            "subgenome": row["subgenome"],
            "status": status,
            "genome_ok": int(required["genome_fasta"].is_file()),
            "annotation_genome_ok": int(required["annotation_genome_fasta"].is_file()),
            "gff_ok": int(required["edta_gff3"].is_file()),
            "library_ok": int(optional["edta_library_fasta"].is_file()),
            "intact_gff_ok": int(optional["edta_intact_gff3"].is_file()),
            "intact_fasta_ok": int(optional["edta_intact_fasta"].is_file()),
            "summary_ok": int(optional["edta_summary"].is_file()),
            "genome_fasta": row["genome_fasta"],
            "edta_gff3": row["edta_gff3"],
            "edta_library_fasta": row["edta_library_fasta"],
            "note": (
                "required independent EDTA inputs present"
                if status == "pass" else
                "stop: genome/GFF3 missing; do not substitute pan-TE files"
            ),
        })
    return reports


def parse_annotation(row: dict[str, str]) -> tuple[dict[str, object], list[dict[str, object]], list[dict[str, object]]]:
    """Parse one EDTA GFF3 and return target summary, family rows and class rows."""
    target = row["target"]
    genome = Path(row["genome_fasta"])
    annotation_genome = Path(row["annotation_genome_fasta"])
    gff = Path(row["edta_gff3"])
    genome_bp = fasta_length(genome)
    seq_regex = re.compile(row["seq_regex"]) if row.get("seq_regex") else None

    family_intervals: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    family_classes: dict[str, set[str]] = defaultdict(set)
    class_intervals: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    all_intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    feature_count = 0
    filtered_count = 0
    for line in gff.open(encoding="utf-8", errors="replace"):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            continue
        seqid, _source, _feature, start_s, end_s, _score, _strand, _phase, attr_text = fields
        attrs = parse_attrs(attr_text)
        if "classification" not in attrs:
            continue
        feature_count += 1
        if seq_regex and not seq_regex.search(seqid):
            continue
        try:
            start, end = int(start_s), int(end_s)
        except ValueError:
            continue
        if end < start:
            start, end = end, start
        family, classification, primary = classify_family(attrs)
        family_intervals[family].append((seqid, start, end))
        family_classes[family].add(classification)
        class_intervals[primary].append((seqid, start, end))
        all_intervals[seqid].append((start, end))
        filtered_count += 1

    fam_rows: list[dict[str, object]] = []
    for family in sorted(family_intervals):
        by_seq: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for seqid, start, end in family_intervals[family]:
            by_seq[seqid].append((start, end))
        union_bp = sum(union_length(vals) for vals in by_seq.values())
        classes = sorted(family_classes[family])
        fam_rows.append({
            "target": target,
            "species": row["species"],
            "subgenome": row["subgenome"],
            "haplotype": row["haplotype"],
            "independent_family_id": family,
            "classification": ";".join(classes),
            "primary_class": ";".join(sorted({c.split("/", 1)[0] for c in classes})),
            "feature_count": len(family_intervals[family]),
            "union_bp": union_bp,
            "genome_fraction": union_bp / genome_bp if genome_bp else 0,
            "classification_conflict": int(len(classes) > 1),
            "library_fasta": row["edta_library_fasta"],
            "edta_gff3": str(gff),
        })

    class_rows: list[dict[str, object]] = []
    for primary in sorted(class_intervals):
        by_seq = defaultdict(list)
        for seqid, start, end in class_intervals[primary]:
            by_seq[seqid].append((start, end))
        union_bp = sum(union_length(vals) for vals in by_seq.values())
        class_rows.append({
            "target": target,
            "species": row["species"],
            "subgenome": row["subgenome"],
            "haplotype": row["haplotype"],
            "primary_class": primary,
            "feature_count": sum(len(v) for v in by_seq.values()),
            "raw_bp": sum(end - start + 1 for vals in by_seq.values() for start, end in vals),
            "union_bp": union_bp,
            "genome_fraction": union_bp / genome_bp if genome_bp else 0,
        })

    te_union_bp = sum(union_length(vals) for vals in all_intervals.values())
    annotation_genome_bp = genome_bp if annotation_genome.resolve() == genome.resolve() else fasta_length(annotation_genome)
    summary = {
        "target": target,
        "species": row["species"],
        "subgenome": row["subgenome"],
        "haplotype": row["haplotype"],
        "genome_bp": genome_bp,
        "annotation_genome_bp": annotation_genome_bp,
        "te_union_bp": te_union_bp,
        "te_fraction": te_union_bp / genome_bp if genome_bp else 0,
        "feature_count": filtered_count,
        "gff_feature_count_before_seq_filter": feature_count,
        "family_count": len(family_intervals),
        "annotation_seq_filter": row.get("seq_regex", ""),
        "source": row["source"],
    }
    return summary, fam_rows, class_rows


def aggregate_group_summaries(summary_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    groups: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in summary_rows:
        groups[(str(row["species"]), str(row["subgenome"]))].append(row)
    out = []
    for (species, subgenome), rows in sorted(groups.items()):
        def vals(field: str) -> list[float]:
            return [float(r[field]) for r in rows]
        genomes, tes, fracs = vals("genome_bp"), vals("te_union_bp"), vals("te_fraction")
        out.append({
            "species": species,
            "subgenome": subgenome,
            "targets": ",".join(str(r["target"]) for r in rows),
            "n_targets": len(rows),
            "mean_genome_bp": sum(genomes) / len(genomes),
            "mean_te_union_bp": sum(tes) / len(tes),
            "mean_te_fraction": sum(fracs) / len(fracs),
            "min_te_fraction": min(fracs),
            "max_te_fraction": max(fracs),
        })
    return out


def run(args: argparse.Namespace) -> int:
    project_root = Path(args.project_root).resolve()
    out_dir = Path(args.output_dir).resolve() if args.output_dir else project_root / "0903" / "TE" / "independent_edta"
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = out_dir / "input_manifest.tsv"
    if args.manifest:
        manifest_path = Path(args.manifest).resolve()
        rows = read_tsv(manifest_path)
    else:
        rows = default_rows(project_root)
        fields = list(rows[0])
        write_tsv(manifest_path, rows, fields)

    reports = validate_manifest(rows)
    write_tsv(out_dir / "validation.tsv", reports, list(reports[0]))
    failed = [r for r in reports if r["status"] != "pass"]
    if failed and args.strict:
        print(f"ERROR: {len(failed)} target(s) missing required independent EDTA inputs", file=sys.stderr)
        print(f"See {out_dir / 'validation.tsv'}", file=sys.stderr)
        return 2

    if args.mode == "manifest":
        print(f"Manifest: {manifest_path}")
        print(f"Validation: {out_dir / 'validation.tsv'}")
        return 0

    summary_rows: list[dict[str, object]] = []
    family_rows: list[dict[str, object]] = []
    class_rows: list[dict[str, object]] = []
    for row in rows:
        if not Path(row["genome_fasta"]).is_file() or not Path(row["edta_gff3"]).is_file():
            continue
        summary, families, classes = parse_annotation(row)
        summary_rows.append(summary)
        family_rows.extend(families)
        class_rows.extend(classes)

    if not summary_rows:
        print("ERROR: no valid independent EDTA GFF3/genome pairs", file=sys.stderr)
        return 2
    write_tsv(out_dir / "subgenome_summary.tsv", summary_rows, list(summary_rows[0]))
    write_tsv(out_dir / "subgenome_group_summary.tsv", aggregate_group_summaries(summary_rows), [
        "species", "subgenome", "targets", "n_targets", "mean_genome_bp",
        "mean_te_union_bp", "mean_te_fraction", "min_te_fraction", "max_te_fraction",
    ])
    write_tsv(out_dir / "independent_family_catalogue.tsv", family_rows, list(family_rows[0]))
    write_tsv(out_dir / "independent_class_summary.tsv", class_rows, list(class_rows[0]))
    print(f"Wrote independent EDTA summaries to {out_dir}")
    print(f"Targets summarized: {len(summary_rows)}/{len(rows)}")
    if failed:
        print("WARNING: optional or required files are missing for some targets; see validation.tsv")
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=".", help="Prunus project root")
    parser.add_argument("--output-dir", default="", help="Output directory (default: 0903/TE/independent_edta)")
    parser.add_argument("--manifest", default="", help="Use an existing manifest TSV")
    parser.add_argument("--mode", choices=("manifest", "summarize", "all"), default="all")
    parser.add_argument("--strict", action=argparse.BooleanOptionalAction, default=True,
                        help="Fail if a genome or EDTA GFF3 is missing (default: true)")
    args = parser.parse_args()
    raise SystemExit(run(args))


if __name__ == "__main__":
    main()
