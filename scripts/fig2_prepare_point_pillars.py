#!/usr/bin/env python3
"""Build GTR-centered POInT pillar candidates from formal quota anchors.

This is an input-preparation step only.  It records observable A/B states and
does not assign a subgenome, estimate branch losses, or calculate a posterior.
Those quantities must come from the POInT likelihood run.
"""

from __future__ import print_function

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


TRACKS = ("PZI_A", "PZI_B", "PSP_A", "PSP_B", "PPA_A", "PPA_B")


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quota-anchors", required=True)
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    anchors_path = Path(args.quota_anchors)
    catalog = Path(args.catalog)
    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit("ERROR: output directory is not empty: {}".format(outdir))
    if not anchors_path.is_file() or anchors_path.stat().st_size == 0:
        raise SystemExit("ERROR: missing formal quota anchors: {}".format(anchors_path))
    if "alltaxa_fractionation_states" in str(anchors_path):
        raise SystemExit("ERROR: candidate JCVI states cannot be POInT input")

    gtr_order = read_tsv(catalog / "gtr_gene_order.tsv")
    anchors = defaultdict(lambda: defaultdict(set))
    for row in read_tsv(anchors_path):
        track = row.get("track_id", "")
        gene = row.get("gtr_gene", "")
        query = row.get("query_gene", "")
        if track in TRACKS and gene and query:
            anchors[gene][track].add(query)
    if not anchors:
        raise SystemExit("ERROR: no PZI/PSP/PPA A/B anchors found")

    outdir.mkdir(parents=True, exist_ok=True)
    fields = ["pillar_id", "gtr_global_index", "gtr_seqid", "gtr_gene"]
    for track in TRACKS:
        fields.extend([track + "_state", track + "_gene"])
    observed_rows = []
    excluded_rows = []
    for index, gtr in enumerate(gtr_order):
        gene = gtr["gtr_gene"]
        row = {
            "pillar_id": "GTR_PILLAR_{:08d}".format(index + 1),
            "gtr_global_index": gtr["gtr_global_index"],
            "gtr_seqid": gtr["gtr_seqid"],
            "gtr_gene": gene,
        }
        for track in TRACKS:
            genes = sorted(anchors[gene].get(track, set()))
            row[track + "_state"] = "1" if genes else "0"
            row[track + "_gene"] = ",".join(genes)
        species_present = []
        for species in ("PZI", "PSP", "PPA"):
            a = row[species + "_A_state"] == "1"
            b = row[species + "_B_state"] == "1"
            species_present.append(a or b)
        if all(species_present):
            observed_rows.append(row)
        else:
            audit = dict(row)
            audit["exclusion_reason"] = "missing_all_copies_in_a_species"
            excluded_rows.append(audit)
    if not observed_rows:
        raise SystemExit("ERROR: no POInT-eligible pillars remain after collinear-block observability checks")
    unambiguous_rows = []
    for row in observed_rows:
        if any("," in row[track + "_gene"] for track in TRACKS):
            audit = dict(row)
            audit["exclusion_reason"] = "multiple_query_genes_in_a_track"
            excluded_rows.append(audit)
        else:
            unambiguous_rows.append(row)
    reuse_counts = {
        track: Counter(row[track + "_gene"] for row in unambiguous_rows if row[track + "_gene"])
        for track in TRACKS
    }
    rows = []
    for row in unambiguous_rows:
        if any(row[track + "_gene"] and reuse_counts[track][row[track + "_gene"]] > 1 for track in TRACKS):
            audit = dict(row)
            audit["exclusion_reason"] = "query_gene_reused_across_pillars"
            excluded_rows.append(audit)
        else:
            rows.append(row)
    if not rows:
        raise SystemExit("ERROR: no unambiguous one-gene-per-track pillars remain")
    write_tsv(outdir / "gtr_centered_pillars.tsv", fields, rows)
    write_tsv(outdir / "excluded_pillar_candidates.tsv", fields + ["exclusion_reason"], excluded_rows)
    # POInT models duplicated genomes only. GTR defines pillar homology and
    # orientation, but is not a duplicated terminal in the likelihood tree.
    (outdir / "point_species_tree.nwk").write_text("(PPA,(PSP,PZI));\n", encoding="ascii")
    (outdir / "point_input_contract.tsv").write_text(
        "field\tdefinition\n"
        "pillar_id\tUnique GTR-centered homologous locus\n"
        "*_state\tObservable quota-filtered anchor presence, 1 or 0\n"
        "*_gene\tZero or one stable prefixed query gene identifier observed at the locus\n"
        "eligibility\tEach PZI/PSP/PPA taxon has at least one copy; no multi-anchor track or reused query gene\n"
        "point_species_tree.nwk\tFrozen duplicated-taxon scaffold; GTR is excluded from the likelihood tree\n",
        encoding="ascii",
    )
    (outdir / "method_status.txt").write_text(
        "Unambiguous pillar candidates were prepared from formal quota-filtered anchors. Each A/B track contributes at most one gene and each query gene occurs in at most one pillar. No POInT likelihood, subgenome assignment, posterior, or branch-loss count was calculated by this script.\n"
        "GTR is an external locus anchor and is not a duplicated POInT terminal.\n"
        "Run the official POInT implementation using its documented pillar/collinear-block format, then normalize its branch output before plotting.\n",
        encoding="ascii",
    )
    (outdir / "COMPLETE").write_text("POInT pillar preparation completed; likelihood inference remains pending.\n", encoding="ascii")
    print("[OK] Wrote {} eligible and {} excluded GTR-centered pillar candidates: {}".format(len(rows), len(excluded_rows), outdir))


if __name__ == "__main__":
    main()
