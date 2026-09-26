#!/usr/bin/env python3
"""Combine formal WGD v2 Ks audits into Figure 2a source tables."""

from __future__ import print_function

import argparse
import csv
import math
import statistics
from pathlib import Path


UPPER = ("PZI", "PSP", "PPA")
LOWER_TRACKS = ("PZI_A", "PZI_B", "PSP_A", "PSP_B", "PPA_A", "PPA_B", "PCAM", "PPE")
LOWER_SPECIES = {"PZI_A": "PZI", "PZI_B": "PZI", "PSP_A": "PSP", "PSP_B": "PSP", "PPA_A": "PPA", "PPA_B": "PPA", "PCAM": "PCAM", "PPE": "PPE"}


def read_tsv(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def finite_ks(row, path):
    try:
        value = float(row["dS"])
    except (KeyError, TypeError, ValueError):
        raise SystemExit("[ERROR] Invalid dS in {}: {}".format(path, row))
    if not math.isfinite(value) or not 0.0 < value < 5.0:
        raise SystemExit("[ERROR] Non-filtered dS in {}: {}".format(path, value))
    return value


def require_complete(path):
    if not (path / "COMPLETE").is_file():
        raise SystemExit("[ERROR] Missing completion marker: {}".format(path / "COMPLETE"))


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upper-root", type=Path, required=True)
    parser.add_argument("--lower-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.outdir.exists() and any(args.outdir.iterdir()):
        raise SystemExit("[ERROR] Output directory is not empty: {}".format(args.outdir))
    args.outdir.mkdir(parents=True, exist_ok=True)

    upper_rows = []
    lower_rows = []
    all_source = []
    summary = []

    for species in UPPER:
        audit_dir = args.upper_root / species / "audit"
        require_complete(args.upper_root / species)
        path = audit_dir / (species + "_AB.ks_0_5.tsv")
        rows = read_tsv(path)
        values = []
        for row in rows:
            value = finite_ks(row, path)
            values.append(value)
            record = {
                "comparison": species,
                "species": species,
                "source_track": species + "_AB",
                "family_id": row["family_id"],
                "gene_a": row.get("A", row.get("gene1", "")),
                "gene_b": row.get("B", row.get("gene2", "")),
                "Ks": "{:.8f}".format(value),
            }
            upper_rows.append(record)
            all_source.append({"panel": "upper", **record})
        summary.append({
            "panel": "upper", "comparison": species, "source_tracks": species + "_AB",
            "n": len(values), "Ks_median": "{:.8f}".format(statistics.median(values)),
            "Ks_mean": "{:.8f}".format(statistics.mean(values)), "filter": "0<Ks<5",
        })

    for track in LOWER_TRACKS:
        species = LOWER_SPECIES[track]
        audit_dir = args.lower_root / track / "audit"
        require_complete(args.lower_root / track)
        path = audit_dir / (track + "_GTR.ks_0_5.tsv")
        rows = read_tsv(path)
        values = []
        for row in rows:
            value = finite_ks(row, path)
            values.append(value)
            record = {
                "comparison": species,
                "species": species,
                "source_track": track + "_GTR",
                "family_id": row["family_id"],
                "gene_a": row.get("gene1", ""),
                "gene_b": row.get("gene2", ""),
                "Ks": "{:.8f}".format(value),
            }
            lower_rows.append(record)
            all_source.append({"panel": "lower", **record})
        summary.append({
            "panel": "lower", "comparison": species, "source_tracks": track + "_GTR",
            "n": len(values), "Ks_median": "{:.8f}".format(statistics.median(values)),
            "Ks_mean": "{:.8f}".format(statistics.mean(values)), "filter": "0<Ks<5",
        })

    write_tsv(args.outdir / "panel_a_upper_ks.tsv", ["comparison", "species", "source_track", "family_id", "gene_a", "gene_b", "Ks"], upper_rows)
    write_tsv(args.outdir / "panel_a_lower_ks.tsv", ["comparison", "species", "source_track", "family_id", "gene_a", "gene_b", "Ks"], lower_rows)
    write_tsv(args.outdir / "panel_a_all_source_data.tsv", ["panel", "comparison", "species", "source_track", "family_id", "gene_a", "gene_b", "Ks"], all_source)
    write_tsv(args.outdir / "panel_a_summary.tsv", ["panel", "comparison", "source_tracks", "n", "Ks_median", "Ks_mean", "filter"], summary)
    (args.outdir / "method_contract.txt").write_text(
        "Formal WGD v2 pairwise CodeML dS values only; all inputs were generated with BLASTP E<=1e-10, MCL I=2.0, MUSCLE 3.8.1551 and PAL2NAL.\n"
        "Upper curves: PZI_AB, PSP_AB and PPA_AB i-ADHoRe anchors after MCL and CDS validation.\n"
        "Lower curves: GTR quota-filtered anchors; A/B tracks are combined within PZI, PSP and PPA.\n"
        "Values retained for plotting: 0<Ks<5; density is normalized within each comparison.\n",
        encoding="ascii",
    )
    (args.outdir / "COMPLETE").write_text("Figure 2a formal Ks source tables completed.\n", encoding="ascii")
    print("[OK] Combined Figure 2a source tables: {}".format(args.outdir))
    print("[OK] Upper rows: {}; lower rows: {}".format(len(upper_rows), len(lower_rows)))


if __name__ == "__main__":
    main()
