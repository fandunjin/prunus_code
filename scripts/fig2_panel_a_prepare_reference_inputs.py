#!/usr/bin/env python3
"""Prepare auditable Figure 2a inputs for the reference-method Ks workflow."""

from __future__ import print_function

import argparse
import csv
from collections import OrderedDict
from pathlib import Path


UPPER_SPECS = ("PZI", "PSP", "PPA")
LOWER_TRACKS = ("PZI_A", "PZI_B", "PSP_A", "PSP_B", "PPA_A", "PPA_B", "PCAM", "PPE")


def read_fasta(path):
    records = OrderedDict()
    identifier = None
    chunks = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    records[identifier] = "".join(chunks).upper()
                identifier = line[1:].split(None, 1)[0]
                if not identifier:
                    raise ValueError("Empty FASTA identifier at {}:{}".format(path, line_number))
                if identifier in records:
                    raise ValueError("Duplicate FASTA identifier {} in {}".format(identifier, path))
                chunks = []
            else:
                if identifier is None:
                    raise ValueError("Sequence before FASTA header in {}".format(path))
                chunks.append(line)
    if identifier is not None:
        records[identifier] = "".join(chunks).upper()
    if not records:
        raise ValueError("No FASTA records in {}".format(path))
    return records


def merge_records(paths):
    merged = OrderedDict()
    for path in paths:
        for identifier, sequence in read_fasta(path).items():
            if identifier in merged:
                raise ValueError("Duplicate identifier {} across {}".format(identifier, path))
            merged[identifier] = sequence
    return merged


def write_fasta(path, records):
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for identifier, sequence in records.items():
            handle.write(">{0}\n".format(identifier))
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def read_tsv(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def select_valid(cds_records, pep_records):
    shared = set(cds_records) & set(pep_records)
    invalid_frame = {identifier for identifier in shared if len(cds_records[identifier]) % 3 != 0}
    empty = {
        identifier for identifier in shared
        if not cds_records[identifier] or not pep_records[identifier]
    }
    valid = shared - invalid_frame - empty
    return valid, invalid_frame, empty, set(cds_records) - set(pep_records), set(pep_records) - set(cds_records)


def upper_pair_rows(audit_rows, species, valid):
    deduplicated = OrderedDict()
    for row in audit_rows:
        if row["species"] != species:
            continue
        gene_a, gene_b = row["gene_a_original"], row["gene_b_original"]
        if gene_a not in valid or gene_b not in valid:
            continue
        key = (gene_a, gene_b)
        if key not in deduplicated:
            deduplicated[key] = {
                "species": species,
                "pair_id": "{}_AB_{:06d}".format(species, len(deduplicated) + 1),
                "gene_a": gene_a,
                "gene_b": gene_b,
                "representative_multiplicon": row["multiplicon"],
                "duplicate_anchor_records": 1,
            }
        else:
            deduplicated[key]["duplicate_anchor_records"] += 1
    return list(deduplicated.values())


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--iadhore-audit", type=Path, required=True)
    parser.add_argument("--quota-anchors", type=Path, required=True)
    parser.add_argument("--tracks", type=Path, required=True)
    parser.add_argument("--pzi-a-cds", type=Path, required=True)
    parser.add_argument("--pzi-b-cds", type=Path, required=True)
    parser.add_argument("--pzi-a-pep", type=Path, required=True)
    parser.add_argument("--pzi-b-pep", type=Path, required=True)
    parser.add_argument("--psp-a-cds", type=Path, required=True)
    parser.add_argument("--psp-b-cds", type=Path, required=True)
    parser.add_argument("--psp-a-pep", type=Path, required=True)
    parser.add_argument("--psp-b-pep", type=Path, required=True)
    parser.add_argument("--ppa-cds", type=Path, required=True)
    parser.add_argument("--ppa-pep", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.outdir.exists() and any(args.outdir.iterdir()):
        raise SystemExit("ERROR: output directory is not empty: {}".format(args.outdir))
    args.outdir.mkdir(parents=True, exist_ok=True)
    upper_dir = args.outdir / "upper"
    lower_dir = args.outdir / "lower"
    upper_dir.mkdir()
    lower_dir.mkdir()

    for path in (args.iadhore_audit, args.quota_anchors):
        if not path.is_file() or not path.stat().st_size:
            raise SystemExit("ERROR: missing input: {}".format(path))
    audit_rows = read_tsv(args.iadhore_audit)
    quota_rows = read_tsv(args.quota_anchors)

    upper_sources = {
        "PZI": ([args.pzi_a_cds, args.pzi_b_cds], [args.pzi_a_pep, args.pzi_b_pep]),
        "PSP": ([args.psp_a_cds, args.psp_b_cds], [args.psp_a_pep, args.psp_b_pep]),
        "PPA": ([args.ppa_cds], [args.ppa_pep]),
    }
    manifest_rows = []
    for species in UPPER_SPECS:
        cds_paths, pep_paths = upper_sources[species]
        cds_records = merge_records(cds_paths)
        pep_records = merge_records(pep_paths)
        valid, invalid_frame, empty, cds_only, pep_only = select_valid(cds_records, pep_records)
        valid_cds = OrderedDict((key, value) for key, value in cds_records.items() if key in valid)
        valid_pep = OrderedDict((key, value) for key, value in pep_records.items() if key in valid)
        write_fasta(upper_dir / (species + "_AB.cds.fa"), valid_cds)
        write_fasta(upper_dir / (species + "_AB.pep.fa"), valid_pep)
        pairs = upper_pair_rows(audit_rows, species, valid)
        if not pairs:
            raise SystemExit("ERROR: no valid i-ADHoRe A-B pairs retained for {}".format(species))
        write_tsv(
            upper_dir / (species + "_AB.iadhore_pairs_pre_mcl.tsv"),
            ["species", "pair_id", "gene_a", "gene_b", "representative_multiplicon", "duplicate_anchor_records"],
            pairs,
        )
        manifest_rows.append({
            "dataset": species,
            "role": "upper_AB_collinear",
            "cds_records_input": len(cds_records),
            "pep_records_input": len(pep_records),
            "valid_records": len(valid),
            "cds_without_peptide": len(cds_only),
            "peptide_without_cds": len(pep_only),
            "nontriplet_cds_excluded": len(invalid_frame),
            "empty_records_excluded": len(empty),
            "candidate_pairs_before_mcl": len(pairs),
        })

    track_records = {}
    for track in ("GTR",) + LOWER_TRACKS:
        cds_path = args.tracks / (track + ".cds")
        pep_path = args.tracks / (track + ".pep")
        if not cds_path.is_file() or not pep_path.is_file():
            raise SystemExit("ERROR: missing track FASTA for {}".format(track))
        cds_records = read_fasta(cds_path)
        pep_records = read_fasta(pep_path)
        valid, invalid_frame, empty, cds_only, pep_only = select_valid(cds_records, pep_records)
        track_records[track] = valid
        manifest_rows.append({
            "dataset": track,
            "role": "lower_GTR_anchor",
            "cds_records_input": len(cds_records),
            "pep_records_input": len(pep_records),
            "valid_records": len(valid),
            "cds_without_peptide": len(cds_only),
            "peptide_without_cds": len(pep_only),
            "nontriplet_cds_excluded": len(invalid_frame),
            "empty_records_excluded": len(empty),
            "candidate_pairs_before_mcl": "NA",
        })

    quota_by_track = {track: OrderedDict() for track in LOWER_TRACKS}
    for row in quota_rows:
        track = row.get("track_id", "")
        if track not in quota_by_track:
            continue
        query_gene, gtr_gene = row["query_gene"], row["gtr_gene"]
        if query_gene not in track_records[track] or gtr_gene not in track_records["GTR"]:
            continue
        key = (query_gene, gtr_gene)
        current = quota_by_track[track].get(key)
        score = float(row.get("quota_score", "0") or 0)
        if current is None or score > float(current["quota_score"]):
            quota_by_track[track][key] = row

    for track in LOWER_TRACKS:
        rows = []
        families = []
        for index, row in enumerate(quota_by_track[track].values(), 1):
            family_id = "{}_GTR_{:06d}".format(track, index)
            rows.append({
                "family_id": family_id,
                "track_id": track,
                "query_gene": row["query_gene"],
                "gtr_gene": row["gtr_gene"],
                "block_id": row["block_id"],
                "quota_score": row["quota_score"],
            })
            families.append({"family_id": family_id, "GTR": row["gtr_gene"], track: row["query_gene"]})
        if not rows:
            raise SystemExit("ERROR: no valid GTR anchors retained for {}".format(track))
        write_tsv(lower_dir / (track + ".quota_pairs.tsv"), list(rows[0].keys()), rows)
        write_tsv(lower_dir / (track + ".families.tsv"), ["family_id", "GTR", track], families)

    write_tsv(
        args.outdir / "input_filtering_manifest.tsv",
        [
            "dataset", "role", "cds_records_input", "pep_records_input", "valid_records",
            "cds_without_peptide", "peptide_without_cds", "nontriplet_cds_excluded",
            "empty_records_excluded", "candidate_pairs_before_mcl",
        ],
        manifest_rows,
    )
    (args.outdir / "method_contract.txt").write_text(
        "Upper data: PZI A1-B1, PSP A1-B1 and PPA A-B i-ADHoRe anchor pairs.\n"
        "Lower data: GTR-centered quota-filtered anchors for the eight displayed tracks.\n"
        "CDS without matching peptide, empty CDS/proteins and non-triplet CDS are excluded.\n"
        "MCL membership is applied to upper A-B anchors in the subsequent step.\n",
        encoding="ascii",
    )
    (args.outdir / "COMPLETE").write_text("Figure 2a reference input preparation completed.\n", encoding="ascii")
    print("[OK] Prepared Figure 2a inputs: {}".format(args.outdir))


if __name__ == "__main__":
    main()
