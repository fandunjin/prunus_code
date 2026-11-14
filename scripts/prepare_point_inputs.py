#!/usr/bin/env python3
"""Prepare native POInT WGD inputs from the formal GTR-anchored pillars.

The GTR locus is used only to define the order of homologous pillars.  It is
not an extant POInT terminal.  PZI, PSP and PPA are represented by two copies
each; the A/B labels are assembly subgenome labels and are not parental
assignments.
"""

from __future__ import print_function

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


SPECIES = ("PZI", "PSP", "PPA")
TRACKS = {
    "PZI": ("PZI_A", "PZI_B"),
    "PSP": ("PSP_A", "PSP_B"),
    "PPA": ("PPA_A", "PPA_B"),
}


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_text(path, text):
    Path(path).write_text(text, encoding="ascii")


def natural_key(value):
    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", value)]


def parse_bed(path):
    rows = []
    seen = set()
    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError(
                    "{}:{} has fewer than four BED fields".format(path, line_number)
                )
            seqid, start, end, gene = fields[:4]
            if gene in seen:
                raise ValueError("{} repeats gene {}".format(path, gene))
            seen.add(gene)
            rows.append((seqid, int(start), int(end), gene))
    rows.sort(key=lambda row: (natural_key(row[0]), row[1], row[2], row[3]))
    return rows


def read_pillars(path):
    rows = read_tsv(path)
    if not rows:
        raise ValueError("empty pillar table: {}".format(path))
    required = ["gtr_global_index"]
    for species in SPECIES:
        required.extend([species + "_A_gene", species + "_B_gene"])
    missing = [field for field in required if field not in rows[0]]
    if missing:
        raise ValueError("pillar table lacks columns: {}".format(", ".join(missing)))
    rows.sort(key=lambda row: int(row["gtr_global_index"]))
    return rows


def write_gene_order(species, track_rows, used_by_track, out_path):
    """Write the two-copy POInT gene order format.

    POInT expects one genome name followed by ``contig_number gene_name``
    records.  Contig numbers are monotonically assigned across the A and B
    tracks; unselected genes are omitted because POInT only needs the genes
    represented in the pillar file.
    """
    lines = [species]
    contig_number = 0
    written = set()
    for track in TRACKS[species]:
        grouped = defaultdict(list)
        for seqid, start, end, gene in track_rows[track]:
            if gene in used_by_track[track]:
                grouped[seqid].append((start, end, gene))
        for seqid in sorted(grouped, key=natural_key):
            genes = sorted(grouped[seqid], key=lambda item: (item[0], item[1], item[2]))
            for start, end, gene in genes:
                if gene in written:
                    raise ValueError(
                        "{} gene reused in gene-order file: {}".format(species, gene)
                    )
                written.add(gene)
                lines.append("{}\t{}".format(contig_number, gene))
            contig_number += 1
    expected = set()
    for track in TRACKS[species]:
        expected.update(used_by_track[track])
    if written != expected:
        missing = sorted(expected - written)
        extra = sorted(written - expected)
        raise ValueError(
            "{} gene-order mismatch; missing={} extra={}".format(
                species, missing[:5], extra[:5]
            )
        )
    write_text(out_path, "\n".join(lines) + "\n")
    return len(written), contig_number


def write_pillar_file(rows, out_path):
    lines = ["\t".join(SPECIES)]
    for row in rows:
        tokens = []
        for species in SPECIES:
            tokens.append(row[species + "_A_gene"] or "NONE")
        # Read_WGX_Data::get_data() consumes copy 0 for all taxa, then copy 1.
        for species in SPECIES:
            tokens.append(row[species + "_B_gene"] or "NONE")
        if len(tokens) != 6:
            raise AssertionError("POInT WGD pillar must contain six copy tokens")
        lines.append("\t".join(tokens))
    write_text(out_path, "\n".join(lines) + "\n")


def write_pillar_map(rows, out_path):
    fields = [
        "optimized_index",
        "pillar_id",
        "gtr_global_index",
        "gtr_seqid",
        "gtr_gene",
    ]
    with Path(out_path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for index, row in enumerate(rows):
            writer.writerow(
                {
                    "optimized_index": str(index),
                    "pillar_id": row.get("pillar_id", ""),
                    "gtr_global_index": row.get("gtr_global_index", ""),
                    "gtr_seqid": row.get("gtr_seqid", ""),
                    "gtr_gene": row.get("gtr_gene", ""),
                }
            )


def write_nexus(out_path):
    write_text(
        out_path,
        "#NEXUS\n\n"
        "Begin trees;\n"
        "Translate\n"
        "  1 PPA,\n"
        "  2 PSP,\n"
        "  3 PZI\n"
        ";\n"
        "tree fixed = [&R] (1:0.1,(2:0.1,3:0.1):0.1):0.1;\n"
        "End;\n",
    )


def write_model(path, biased):
    if biased:
        name = "PZI_PSP_PPA_WGD_bias"
        parameter_lines = [
            "FixRate ZERO_TO_INF 0.0612",
            "SwitchProb ZERO_TO_ONE 0.0064",
            "Copy2Bias ZERO_TO_ONE 0.7299",
        ]
        matrix = [
            ["Null", "FixRate", "Default", "Copy2Bias"],
            ["Zero", "Null", "Zero", "Zero"],
            ["Zero", "Zero", "Null", "Zero"],
            ["Zero", "Zero", "Zero", "Null"],
        ]
    else:
        name = "PZI_PSP_PPA_WGD_symmetric"
        parameter_lines = [
            "FixRate ZERO_TO_INF 0.0612",
            "SwitchProb ZERO_TO_ONE 0.0064",
        ]
        matrix = [
            ["Null", "FixRate", "Default", "Default"],
            ["Zero", "Null", "Zero", "Zero"],
            ["Zero", "Zero", "Null", "Zero"],
            ["Zero", "Zero", "Zero", "Null"],
        ]
    lines = [
        name + " hierarchical 2",
        "#NumStates",
        "4",
        "#NumParams",
        str(len(parameter_lines)),
        "#StateNames",
        "Dupl 1 1",
        "DuplFix 1 1",
        "Copy1 1 0",
        "Copy2 0 1",
        "#ParameterDescripts",
    ]
    lines.extend(parameter_lines)
    lines.append("#Matrix")
    lines.extend("\t".join(row) for row in matrix)
    # The symmetric comparator must carry the same root-state and redundancy
    # constraints as the published WGD_bias_fix model.  Otherwise its lnL is
    # inflated by a different root model and the two fits are not comparable.
    lines.extend(
        [
            "#OtherLimitations",
            "Redundancy\tDupl\tDuplFix",
            "rootstate\tDupl",
        ]
    )
    lines.append("")
    write_text(path, "\n".join(lines))


def validate_inputs(outdir):
    pillar_path = Path(outdir) / "point_pillars_initial.txt"
    if not pillar_path.is_file() or pillar_path.stat().st_size == 0:
        raise ValueError("missing written POInT pillar file")
    lines = pillar_path.read_text(encoding="ascii").splitlines()
    if not lines or lines[0].split("\t") != list(SPECIES):
        raise ValueError("invalid POInT pillar header")
    expected_tracks = ("PZI_A", "PSP_A", "PPA_A", "PZI_B", "PSP_B", "PPA_B")
    for line_number, line in enumerate(lines[1:], 2):
        tokens = line.split("\t")
        if len(tokens) != 6:
            raise ValueError(
                "pillar line {} does not have six fields".format(line_number)
            )
        # Read_WGX_Data uses a copy-major layout: PZI_A, PSP_A, PPA_A,
        # followed by PZI_B, PSP_B, PPA_B.  Validate the formal prefixes so
        # a sequential species-pair layout cannot silently corrupt a fit.
        for token, track in zip(tokens, expected_tracks):
            if token != "NONE" and not token.startswith(track + "__"):
                raise ValueError(
                    "pillar line {} has {} in the {} column; expected {}__ or NONE".format(
                        line_number, token, track, track
                    )
                )
    for species in SPECIES:
        path = Path(outdir) / (species + "_POInT_geneorder.txt")
        tokens = path.read_text(encoding="ascii").split()
        if not tokens or tokens[0] != species:
            raise ValueError("invalid gene order for {}".format(species))
    return len(lines) - 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pillars", required=True)
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit("ERROR: output directory is not empty: {}".format(outdir))
    outdir.mkdir(parents=True, exist_ok=True)

    rows = read_pillars(args.pillars)
    catalog = Path(args.catalog)
    track_rows = {}
    used_by_track = defaultdict(set)
    for row in rows:
        for species in SPECIES:
            for copy in ("A", "B"):
                gene = row[species + "_" + copy + "_gene"].strip()
                if gene:
                    used_by_track[species + "_" + copy].add(gene)
    for track in sorted(used_by_track):
        bed_path = catalog / "tracks" / (track + ".bed")
        if not bed_path.is_file():
            raise SystemExit("ERROR: missing formal track BED: {}".format(bed_path))
        track_rows[track] = parse_bed(bed_path)
        available = {row[3] for row in track_rows[track]}
        missing = sorted(used_by_track[track] - available)
        if missing:
            raise SystemExit(
                "ERROR: pillar genes missing from {} BED: {}".format(
                    track, ",".join(missing[:5])
                )
            )

    for species in SPECIES:
        count, contigs = write_gene_order(
            species,
            track_rows,
            used_by_track,
            outdir / (species + "_POInT_geneorder.txt"),
        )
        write_text(
            outdir / (species + "_geneorder_stats.tsv"),
            "field\tvalue\n" "genes\t{}\ncontigs\t{}\n".format(count, contigs),
        )

    write_pillar_file(rows, outdir / "point_pillars_initial.txt")
    write_pillar_map(rows, outdir / "initial_pillar_map.tsv")
    write_nexus(outdir / "point_species_tree.nex")
    write_model(outdir / "WGD_symmetric_model.txt", biased=False)
    write_model(outdir / "WGD_biased_model.txt", biased=True)
    n = validate_inputs(outdir)
    write_text(
        outdir / "input_contract.tsv",
        "field\tvalue\n"
        "species_order\tPZI,PSP,PPA\n"
        "copy_order\tall species copy 0 (A), then all species copy 1 (B); A/B are assembly subgenome labels\n"
        "external_anchor\tGTR only; not a POInT terminal\n"
        "pillar_count\t{}\n"
        "topology\t(PPA,(PSP,PZI))\n"
        "model_pair\tsymmetric versus copy-specific loss rate\n".format(n),
    )
    write_text(outdir / "PREPARED", "Native POInT inputs prepared and validated.\n")
    print("[OK] Prepared native POInT inputs: {} pillars".format(n))


if __name__ == "__main__":
    main()
