#!/usr/bin/env python3
"""Freeze the gene-order and protein catalogs used by Fig. 2c/d.

The reference workflow aligns one combined query catalog to a low-copy target.
Every query track receives a stable prefix so that A/B and species labels cannot
collide during LAST, Quota Align, or pillar construction.  This script only
prepares inputs; it does not infer anchors, gene loss, or subgenome identity.
"""

from __future__ import print_function

import argparse
import csv
import re
from collections import OrderedDict
from pathlib import Path


TRACKS = [
    ("PZI_A", "pzi_hapA1", "PZI", "A"),
    ("PZI_B", "pzi_hapB1", "PZI", "B"),
    ("PSP_A", "psp_hapA1", "PSP", "A"),
    ("PSP_B", "psp_hapB1", "PSP", "B"),
    ("PPA_A", "ppa_hapA", "PPA", "A"),
    ("PPA_B", "ppa_hapB", "PPA", "B"),
    ("PCAM", "pcam", "PCAM", "single_copy_control"),
    ("PPE", "ppe", "PPE", "single_copy_control"),
]

PZI_DISPLAY_PAIRS = [
    ("chr01+chr13", "C1A", "C6B"),
    ("chr02+chr15", "C2A", "C7B"),
    ("chr03+chr12", "C7A", "C1B"),
]

GTR_PRIMARY_CHROMOSOME_RE = re.compile(r"^Chr[0-9]+$")


def read_bed(path):
    rows = []
    seen = set()
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError("BED has fewer than four columns at {}:{}".format(path, line_number))
            if fields[3] in seen:
                raise ValueError("Duplicate BED gene {} in {}".format(fields[3], path))
            seen.add(fields[3])
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                raise ValueError("Non-integer BED coordinates at {}:{}".format(path, line_number))
            rows.append({
                "seqid": fields[0], "start": start, "end": end,
                "gene": fields[3], "score": fields[4] if len(fields) > 4 else "0",
                "strand": fields[5] if len(fields) > 5 else ".",
            })
    if not rows:
        raise ValueError("Empty BED: {}".format(path))
    return rows


def read_fasta(path):
    records = OrderedDict()
    current = None
    chunks = []
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    records[current] = "".join(chunks)
                current = line[1:].split()[0]
                if current in records:
                    raise ValueError("Duplicate FASTA record {} in {}".format(current, path))
                chunks = []
            else:
                chunks.append(line)
    if current is not None:
        records[current] = "".join(chunks)
    if not records:
        raise ValueError("Empty FASTA: {}".format(path))
    return records


def write_fasta(path, records):
    with Path(path).open("w", encoding="utf-8", newline="\n") as handle:
        for name, sequence in records.items():
            handle.write(">{}\n".format(name))
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_track(label, source, outdir, is_reference=False):
    bed_path = source / (source.name + ".bed")
    pep_path = source / (source.name + ".pep")
    cds_path = source / (source.name + ".cds")
    if not bed_path.is_file() or bed_path.stat().st_size == 0:
        raise ValueError("Missing BED: {}".format(bed_path))
    if not pep_path.is_file() or pep_path.stat().st_size == 0:
        raise ValueError("Missing protein FASTA: {}".format(pep_path))
    if not cds_path.is_file() or cds_path.stat().st_size == 0:
        raise ValueError("Missing CDS FASTA: {}".format(cds_path))
    bed_rows = read_bed(bed_path)
    peptides = read_fasta(pep_path)
    coding_sequences = read_fasta(cds_path)
    prefix = "GTR__" if is_reference else label + "__"
    output_bed_rows = []
    output_peptides = OrderedDict()
    output_cds = OrderedDict()
    missing = []
    chromosome_order = OrderedDict()
    for order, row in enumerate(bed_rows):
        if is_reference and not GTR_PRIMARY_CHROMOSOME_RE.match(row["seqid"]):
            continue
        if row["gene"] not in peptides or row["gene"] not in coding_sequences:
            missing.append(row["gene"])
            continue
        chromosome_order.setdefault(row["seqid"], len(chromosome_order))
        output_bed_rows.append({
            "seqid": prefix + row["seqid"],
            "start": row["start"],
            "end": row["end"],
            "gene": prefix + row["gene"],
            "score": row["score"],
            "strand": row["strand"],
            "order": order,
        })
        output_peptides[prefix + row["gene"]] = peptides[row["gene"]]
        output_cds[prefix + row["gene"]] = coding_sequences[row["gene"]]
    if missing:
        preview = ", ".join(missing[:5])
        raise ValueError("{} BED genes lack matching protein records in {} (e.g. {})".format(len(missing), pep_path, preview))
    output_bed = outdir / (label + ".bed")
    output_pep = outdir / (label + ".pep")
    output_cds_path = outdir / (label + ".cds")
    with output_bed.open("w", encoding="utf-8", newline="\n") as handle:
        for row in output_bed_rows:
            handle.write("\t".join(str(row[key]) for key in ("seqid", "start", "end", "gene", "score", "strand")) + "\n")
    write_fasta(output_pep, output_peptides)
    write_fasta(output_cds_path, output_cds)
    seqids = outdir / (label + ".seqids")
    with seqids.open("w", encoding="utf-8", newline="\n") as handle:
        for seqid in chromosome_order:
            handle.write(prefix + seqid + "\n")
    return output_bed_rows, output_peptides, output_cds, prefix


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs-root", required=True, help="downstream_analysis/ks/jcvi/inputs")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    inputs = Path(args.inputs_root)
    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit("ERROR: output directory is not empty: {}".format(outdir))
    tracks_dir = outdir / "tracks"
    tracks_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows = []
    all_query_bed = []
    all_query_peptides = OrderedDict()
    all_query_cds = OrderedDict()
    all_track_rows = []
    for track_id, source_label, species, copy_state in TRACKS:
        rows, peptides, coding_sequences, prefix = write_track(
            track_id, inputs / source_label, tracks_dir, False
        )
        all_query_bed.extend(rows)
        all_query_peptides.update(peptides)
        all_query_cds.update(coding_sequences)
        for order, row in enumerate(rows):
            all_track_rows.append({
                "track_id": track_id,
                "species": species,
                "copy_state": copy_state,
                "query_gene": row["gene"],
                "query_seqid": row["seqid"],
                "query_order": order,
            })
        manifest_rows.append({
            "track_id": track_id,
            "source_label": source_label,
            "species": species,
            "copy_state": copy_state,
            "bed": str((tracks_dir / (track_id + ".bed")).resolve()),
            "pep": str((tracks_dir / (track_id + ".pep")).resolve()),
            "cds": str((tracks_dir / (track_id + ".cds")).resolve()),
            "n_bed_genes": len(rows),
            "n_peptides": len(peptides),
            "n_cds": len(coding_sequences),
        })

    gtr_rows, gtr_peptides, gtr_cds, gtr_prefix = write_track(
        "GTR", inputs / "gtr", tracks_dir, True
    )
    all_query_bed_path = outdir / "query8.bed"
    with all_query_bed_path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in all_query_bed:
            handle.write("\t".join(str(row[key]) for key in ("seqid", "start", "end", "gene", "score", "strand")) + "\n")
    write_fasta(outdir / "query8.pep", all_query_peptides)
    write_fasta(outdir / "query8.cds", all_query_cds)
    with (outdir / "query8.seqids").open("w", encoding="utf-8", newline="\n") as handle:
        for track_id, _, _, _ in TRACKS:
            seqids = tracks_dir / (track_id + ".seqids")
            for line in seqids.read_text(encoding="utf-8").splitlines():
                handle.write(line + "\n")
    # GTR global order is the x-axis for panel c; chromosome boundaries remain explicit.
    gtr_order_rows = []
    index = 0
    last_seqid = None
    chromosome_number = 0
    for row in gtr_rows:
        if row["seqid"] != last_seqid:
            chromosome_number += 1
            last_seqid = row["seqid"]
        gtr_order_rows.append({
            "gtr_global_index": index,
            "gtr_seqid": row["seqid"],
            "gtr_gene": row["gene"],
            "gtr_start": row["start"],
            "gtr_end": row["end"],
            "gtr_chromosome_number": chromosome_number,
        })
        index += 1
    write_tsv(outdir / "gtr_gene_order.tsv", list(gtr_order_rows[0].keys()), gtr_order_rows)
    write_tsv(outdir / "query_gene_track_map.tsv", list(all_track_rows[0].keys()), all_track_rows)
    write_tsv(outdir / "track_manifest.tsv", list(manifest_rows[0].keys()), manifest_rows)
    write_tsv(
        outdir / "gtr_primary_chromosomes.tsv",
        ["gtr_seqid", "n_genes"],
        [
            {
                "gtr_seqid": seqid,
                "n_genes": sum(1 for row in gtr_rows if row["seqid"] == seqid),
            }
            for seqid in sorted({row["seqid"] for row in gtr_rows})
        ],
    )
    write_tsv(outdir / "pzi_display_chromosome_pairs.tsv",
              ["pair_label", "pzi_A_chromosome", "pzi_B_chromosome"],
              [{"pair_label": p, "pzi_A_chromosome": a, "pzi_B_chromosome": b} for p, a, b in PZI_DISPLAY_PAIRS])
    (outdir / "method_contract.txt").write_text(
        "Reference-adapted Fig. 2c/d input catalog. Query tracks are PZI_A, PZI_B, PSP_A, PSP_B, PPA_A, PPA_B, PCAM and PPE; GTR is the low-copy target.\n"
        "Stable prefixes are added before LAST so identifiers cannot collide. GTR is restricted to its nine chromosome-scale sequences Chr01-Chr09; unplaced contigs are excluded. The catalog is not an anchor result and contains no inferred gene-loss calls.\n"
        "The formal next steps are LAST, Quota Align Merge (maximum inter-block distance 20 genes), Quota Align (8:1), FractBias (100-gene window), then POInT on quota-filtered pillars.\n",
        encoding="ascii",
    )
    (outdir / "COMPLETE").write_text("Catalog preparation completed; no formal fractionation inference was performed.\n", encoding="ascii")
    print("[OK] Wrote strict Fig. 2c/d catalogs: {}".format(outdir))


if __name__ == "__main__":
    main()
