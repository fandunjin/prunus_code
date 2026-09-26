#!/usr/bin/env python3
"""Standardize corrected haplotype and reference inputs for five Ks panels.

No prior Figure 2 A/B input, anchor or Ks result is read.  The phased PZI,
PSP and PPA inputs originate solely from ``haploid_repartition/inputs``.
PCAM, PPE and GTR are copied from the project-wide canonical reference set.
Each output label receives BED, CDS, protein and sequence-order files with an
identical set of transcript identifiers.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument(
        "--repartition-root", type=Path,
        help="Directory containing haploid_repartition/inputs (defaults to project-root/haploid_repartition)",
    )
    parser.add_argument(
        "--canonical-root", type=Path,
        help="Canonical PCAM/PPE staging directory (defaults to project-root/downstream_analysis/ks/stage/canonical)",
    )
    parser.add_argument(
        "--gtr-root", type=Path,
        help="Directory containing the GTR annotation package (defaults to project-root/data/gtr)",
    )
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def require_file(path: Path, label: str) -> Path:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError("Missing or empty {}: {}".format(label, path))
    return path


def parse_attributes(text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for item in text.rstrip(";\n").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
    return attrs


def open_text(path: Path):
    """Open plain text or gzip payloads whose extension was not preserved."""
    with path.open("rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def fasta_records(path: Path) -> Iterator[Tuple[str, str]]:
    header = None
    chunks: List[str] = []
    with open_text(path) as handle:
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


def fasta_id(header: str) -> str:
    protein_id = re.search(r"\[protein_id=([^\]]+)\]", header)
    if protein_id:
        return protein_id.group(1)
    return header.split(None, 1)[0]


def fasta_sequence_order(path: Path) -> List[str]:
    return [fasta_id(header) for header, _ in fasta_records(path)]


def fasta_id_set(path: Path) -> Set[str]:
    values: Set[str] = set()
    for header, _ in fasta_records(path):
        identifier = fasta_id(header)
        if identifier in values:
            raise RuntimeError("Duplicate FASTA ID in {}: {}".format(path, identifier))
        values.add(identifier)
    return values


def read_mrna_features(path: Path) -> Tuple[List[Dict[str, object]], Dict[str, Set[str]]]:
    records: List[Dict[str, object]] = []
    seen: Set[str] = set()
    cds_protein_by_parent: Dict[str, Set[str]] = defaultdict(set)
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            attrs = parse_attributes(fields[8])
            feature = fields[2].lower()
            if feature == "cds":
                protein_id = attrs.get("protein_id", "") or attrs.get("Name", "")
                parents = [parent for parent in attrs.get("Parent", "").split(",") if parent]
                if protein_id:
                    for parent in parents:
                        cds_protein_by_parent[parent].add(protein_id)
                continue
            if feature not in {"mrna", "transcript"}:
                continue
            identifier = attrs.get("ID", "")
            if not identifier:
                raise RuntimeError("Missing transcript ID at {}:{}".format(path, line_number))
            if identifier in seen:
                raise RuntimeError("Duplicate transcript ID in {}: {}".format(path, identifier))
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as error:
                raise RuntimeError("Invalid transcript coordinate at {}:{}".format(path, line_number)) from error
            seen.add(identifier)
            records.append({
                "seqid": fields[0], "start0": min(start, end) - 1, "end": max(start, end),
                "gene": identifier, "strand": fields[6] if fields[6] in {"+", "-"} else ".",
            })
    if not records:
        raise RuntimeError("No mRNA/transcript features in {}".format(path))
    return records, cds_protein_by_parent


def harmonize_model_ids(
    records: List[Dict[str, object]],
    cds_ids: Set[str],
    pep_ids: Set[str],
    cds_protein_by_parent: Dict[str, Set[str]],
) -> Tuple[List[Dict[str, object]], int, int]:
    """Map GFF mRNA coordinates to IDs shared by CDS and protein FASTA files.

    Project annotations use mRNA IDs directly. NCBI annotations instead use an
    mRNA ID in the GFF and a distinct protein_id in both FASTA files. The
    latter mapping is only accepted when the CDS Parent relation yields one
    unambiguous FASTA identifier.
    """
    shared = cds_ids.intersection(pep_ids)
    retained: List[Dict[str, object]] = []
    direct = 0
    cds_linked = 0
    seen_output: Set[str] = set()
    for record in records:
        original = str(record["gene"])
        if original in shared:
            output_id = original
            direct += 1
        else:
            candidates = sorted(cds_protein_by_parent.get(original, set()).intersection(shared))
            if len(candidates) == 0:
                continue
            if len(candidates) != 1:
                raise RuntimeError("Ambiguous CDS-to-protein mapping for {}: {}".format(
                    original, ",".join(candidates)))
            output_id = candidates[0]
            cds_linked += 1
        if output_id in seen_output:
            raise RuntimeError("Multiple GFF mRNAs map to FASTA identifier: {}".format(output_id))
        selected = dict(record)
        selected["gene"] = output_id
        retained.append(selected)
        seen_output.add(output_id)
    return retained, direct, cds_linked


def write_fasta(path: Path, source: Path, allowed: Set[str]) -> int:
    found: Set[str] = set()
    with path.open("w", encoding="ascii", newline="\n") as destination:
        for header, sequence in fasta_records(source):
            identifier = fasta_id(header)
            if identifier not in allowed:
                continue
            if identifier in found:
                raise RuntimeError("Duplicate selected FASTA ID in {}: {}".format(source, identifier))
            destination.write(">{}\n".format(identifier))
            for offset in range(0, len(sequence), 80):
                destination.write(sequence[offset:offset + 80] + "\n")
            found.add(identifier)
    missing = sorted(allowed - found)
    if missing:
        raise RuntimeError("{} selected identifiers missing from {}: {}".format(
            len(missing), source, ",".join(missing[:5])))
    return len(found)


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def source_definitions(
    root: Path,
    repartition_root: Optional[Path] = None,
    canonical_root: Optional[Path] = None,
    gtr_root: Optional[Path] = None,
) -> List[Dict[str, object]]:
    repartition = (repartition_root or (root / "haploid_repartition")).resolve() / "inputs"
    output: List[Dict[str, object]] = []
    for species, haplotypes in (("pzi", ("hapA1", "hapA2", "hapB1", "hapB2")),
                                ("psp", ("hapA1", "hapA2", "hapB1", "hapB2")),
                                ("ppa", ("hapA", "hapB"))):
        for haplotype in haplotypes:
            label = "{}_{}".format(species, haplotype)
            prefix = label
            base = repartition / species
            output.append({
                "label": label,
                "kind": "corrected_haplotype",
                "genome": base / "genome" / (prefix + ".fa"),
                "gff": base / "gff" / (prefix + ".gff3"),
                "cds": base / "cds" / (prefix + ".cds.fa"),
                "pep": base / "pep" / (prefix + ".pep.fa"),
            })
    canonical = (canonical_root or (root / "downstream_analysis" / "ks" / "stage" / "canonical")).resolve()
    for label in ("pcam", "ppe"):
        base = canonical / label
        spec: Dict[str, object] = {
            "label": label,
            "kind": "canonical_reference",
            "genome": base / "genome.fa",
            "gff": base / "primary.gff3",
            "cds": base / "primary.cds.fa",
            "pep": base / "primary.pep.fa",
        }
        if label == "ppe":
            # The RefSeq package contains short unplaced sequences in addition
            # to its eight nuclear pseudomolecules. They are excluded from all
            # chromosome-level synteny and self-synteny analyses.
            spec["include_seqids_regex"] = r"^NC_0340(?:09|10|11|12|13|14|15|16)\.1$"
        output.append(spec)
    # GTR FASTA/GFF files in the older canonical staging directory originate
    # from different assembly name spaces. Use the matching raw annotation
    # package and retain its nine pseudomolecules only; JCVI requires BED/CDS/
    # protein inputs, so genomic sequence is not used for this reference.
    gtr = (gtr_root or (root / "data" / "gtr")).resolve()
    output.append({
        "label": "gtr",
        "kind": "canonical_reference_chromosome_only",
        "genome": None,
        "gff": gtr / "Gtr.NCBI.gene.gff3",
        "cds": gtr / "Gtr.NCBI.gene.modId.fna",
        "pep": gtr / "Gtr.NCBI.prot.faa",
        "include_seqids_regex": r"^Chr[0-9]{2}$",
    })
    return output


def build_plans(outdir: Path) -> None:
    external: List[Dict[str, object]] = []
    taxa = {
        "PZI": {"A": ("pzi_hapA1", "pzi_hapA2"), "B": ("pzi_hapB1", "pzi_hapB2")},
        "PSP": {"A": ("psp_hapA1", "psp_hapA2"), "B": ("psp_hapB1", "psp_hapB2")},
        "PPA": {"A": ("ppa_hapA",), "B": ("ppa_hapB",)},
    }
    for reference in ("pcam", "ppe", "gtr"):
        for species, groups in taxa.items():
            for subgenome, labels in groups.items():
                for query in labels:
                    external.append({
                        "comparison": "{}_vs_{}".format(query, reference),
                        "query_label": query,
                        "ref_label": reference,
                        "comparison_type": "cross_species_reference",
                        "panel": "orthologue_{}".format(reference),
                        "curve": "{}_{}".format(species, subgenome),
                    })
    ab = [
        ("pzi_hapA1", "pzi_hapB1", "PZI"), ("pzi_hapA2", "pzi_hapB2", "PZI"),
        ("psp_hapA1", "psp_hapB1", "PSP"), ("psp_hapA2", "psp_hapB2", "PSP"),
        ("ppa_hapA", "ppa_hapB", "PPA"),
    ]
    ab_rows = [{
        "comparison": "{}_vs_{}".format(query, ref), "query_label": query, "ref_label": ref,
        "comparison_type": "within_species_homeolog", "panel": "ab_paralogue", "curve": species,
    } for query, ref, species in ab]
    self_rows: List[Dict[str, object]] = []
    for species, groups in taxa.items():
        for subgenome, labels in groups.items():
            for label in labels:
                self_rows.append({"label": label, "panel": "self_paralogue", "curve": "{}_{}".format(species, subgenome)})
    for label, curve in (("pcam", "PCAM"), ("ppe", "PPE"), ("gtr", "GTR")):
        self_rows.append({"label": label, "panel": "self_paralogue", "curve": curve})
    write_tsv(outdir / "plans" / "external_pair_plan.tsv",
              ["comparison", "query_label", "ref_label", "comparison_type", "panel", "curve"], external)
    write_tsv(outdir / "plans" / "ab_pair_plan.tsv",
              ["comparison", "query_label", "ref_label", "comparison_type", "panel", "curve"], ab_rows)
    write_tsv(outdir / "plans" / "jcvi_pair_plan.tsv",
              ["comparison", "query_label", "ref_label", "comparison_type", "panel", "curve"], external + ab_rows)
    write_tsv(outdir / "plans" / "self_pair_plan.tsv", ["label", "panel", "curve"], self_rows)


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    outdir = args.outdir.resolve()
    if outdir.exists() and any(outdir.iterdir()):
        if not args.force:
            raise RuntimeError("Output directory is not empty: {}".format(outdir))
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    manifest: List[Dict[str, object]] = []
    validation: List[Dict[str, object]] = []
    repartition_root = args.repartition_root.resolve() if args.repartition_root else None
    canonical_root = args.canonical_root.resolve() if args.canonical_root else None
    gtr_root = args.gtr_root.resolve() if args.gtr_root else None
    for spec in source_definitions(root, repartition_root, canonical_root, gtr_root):
        label = str(spec["label"])
        genome_value = spec["genome"]
        genome = require_file(Path(genome_value), label + " genome") if genome_value else None
        gff = require_file(Path(spec["gff"]), label + " GFF3")
        cds = require_file(Path(spec["cds"]), label + " CDS")
        pep = require_file(Path(spec["pep"]), label + " protein")
        mrna, cds_protein_by_parent = read_mrna_features(gff)
        seqid_regex = str(spec.get("include_seqids_regex", ""))
        if seqid_regex:
            rule = re.compile(seqid_regex)
            mrna = [row for row in mrna if rule.fullmatch(str(row["seqid"]))]
            if not mrna:
                raise RuntimeError("{} has no mRNA models after sequence filter {}".format(label, seqid_regex))
        if genome is not None:
            chromosome_order = fasta_sequence_order(genome)
            if seqid_regex:
                chromosome_order = [seqid for seqid in chromosome_order if rule.fullmatch(seqid)]
            if len(chromosome_order) < 2:
                raise RuntimeError("{} has fewer than two genomic sequences".format(label))
            order_index = {seqid: index for index, seqid in enumerate(chromosome_order)}
            bad_seqids = sorted({str(row["seqid"]) for row in mrna if str(row["seqid"]) not in order_index})
            if bad_seqids:
                raise RuntimeError("{} GFF sequence IDs absent from genome: {}".format(label, ",".join(bad_seqids[:5])))
            order_source = "genome_fasta"
        else:
            chromosome_order = sorted({str(row["seqid"]) for row in mrna})
            order_index = {seqid: index for index, seqid in enumerate(chromosome_order)}
            order_source = "GFF_selected_sequence_ids"
        gff_ids = {str(row["gene"]) for row in mrna}
        cds_ids = fasta_id_set(cds)
        pep_ids = fasta_id_set(pep)
        bed_rows, direct_models, cds_linked_models = harmonize_model_ids(
            mrna, cds_ids, pep_ids, cds_protein_by_parent)
        retained_ids = {str(row["gene"]) for row in bed_rows}
        if len(retained_ids) < 100:
            raise RuntimeError("{} retains too few shared GFF/CDS/PEP models".format(label))
        bed_rows.sort(key=lambda row: (order_index[str(row["seqid"])], int(row["start0"]), int(row["end"]), str(row["gene"])))
        target = outdir / "inputs" / label
        target.mkdir(parents=True)
        with (target / (label + ".bed")).open("w", encoding="utf-8", newline="\n") as handle:
            for row in bed_rows:
                handle.write("{}\t{}\t{}\t{}\t0\t{}\n".format(
                    row["seqid"], row["start0"], row["end"], row["gene"], row["strand"]
                ))
        cds_n = write_fasta(target / (label + ".cds"), cds, retained_ids)
        pep_n = write_fasta(target / (label + ".pep"), pep, retained_ids)
        (target / (label + ".seqids")).write_text(",".join(chromosome_order) + "\n", encoding="ascii")
        manifest.append({
            "label": label, "kind": spec["kind"], "genome": genome or "", "gff": gff, "cds": cds, "pep": pep,
            "output_dir": target,
        })
        validation.append({
            "label": label, "chromosomes": len(chromosome_order), "gff_transcripts": len(gff_ids),
            "cds_records": len(cds_ids), "pep_records": len(pep_ids), "retained_models": len(retained_ids),
            "direct_id_models": direct_models, "cds_linked_models": cds_linked_models,
            "excluded_gff_only": len(gff_ids) - len(bed_rows), "excluded_cds_only": len(cds_ids - retained_ids),
            "excluded_pep_only": len(pep_ids - retained_ids), "output_cds": cds_n, "output_pep": pep_n,
            "sequence_order_source": order_source,
            "status": "pass" if cds_n == pep_n == len(bed_rows) else "fail",
        })
    build_plans(outdir)
    write_tsv(outdir / "input_manifest.tsv", ["label", "kind", "genome", "gff", "cds", "pep", "output_dir"], manifest)
    write_tsv(outdir / "input_validation.tsv", [
        "label", "chromosomes", "gff_transcripts", "cds_records", "pep_records", "retained_models",
        "direct_id_models", "cds_linked_models",
        "excluded_gff_only", "excluded_cds_only", "excluded_pep_only", "output_cds", "output_pep",
        "sequence_order_source", "status",
    ], validation)
    (outdir / "method_contract.txt").write_text(
        "Corrected Figure 2 Ks input contract. PZI, PSP and PPA derive only from haploid_repartition/inputs.\n"
        "Each label retains GFF coordinates linked unambiguously to matching CDS and protein FASTA IDs.\n"
        "NCBI PCAM/PPE identifiers use the audited GFF CDS Parent-to-protein_id relation.\n"
        "GTR retains only Chr01-Chr09 from its matching raw annotation package.\n"
        "PZI/PSP phased components remain independent inputs; pooling is deferred until post-Ks plotting.\n"
        "External comparisons and A-B homeolog comparisons use JCVI; self-synteny uses i-ADHoRe.\n",
        encoding="ascii",
    )
    (outdir / "COMPLETE").write_text("Corrected Figure 2 Ks inputs and plans prepared.\n", encoding="ascii")
    print("[OK] Prepared {} input tracks: {}".format(len(validation), outdir))
    for row in validation:
        print("[OK] {label}: retained={retained_models}; chromosomes={chromosomes}; status={status}".format(**row))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as error:
        print("[ERROR] {}".format(error))
        raise SystemExit(2)
