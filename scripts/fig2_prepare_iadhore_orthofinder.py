#!/usr/bin/env python3
"""Prepare an i-ADHoRe 3 input from frozen gene order and orthogroups.

The reference workflow uses gene-family homology assignments for i-ADHoRe,
then reads the resulting cross-copy anchor points before the PUG family screen.
This script only converts identifiers and writes i-ADHoRe lists/configuration;
it does not infer collinearity.
"""

from __future__ import print_function

import argparse
import csv
from collections import OrderedDict
from pathlib import Path


TRACKS = OrderedDict([
    ("PZI_A", ("pzi_A", "PZI")),
    ("PZI_B", ("pzi_B", "PZI")),
    ("PSP_A", ("psp_A", "PSP")),
    ("PSP_B", ("psp_B", "PSP")),
    ("PPA_A", ("ppa_A", "PPA")),
    ("PPA_B", ("ppa_B", "PPA")),
])


def read_id_map(path):
    original_to_stable = {}
    stable_to_original = {}
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            original = row["original_id"]
            stable = row["stable_id"]
            for alias in {original, *[part for part in original.split("|") if part]}:
                if alias in original_to_stable and original_to_stable[alias] != stable:
                    raise ValueError("Non-unique original alias {} in {}".format(alias, path))
                original_to_stable[alias] = stable
            if stable in stable_to_original and stable_to_original[stable] != original:
                raise ValueError("Non-unique stable ID {} in {}".format(stable, path))
            stable_to_original[stable] = original
    return original_to_stable, stable_to_original


def read_bed(path):
    rows = []
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise ValueError("Malformed BED at {}:{}".format(path, line_number))
            rows.append({"seqid": fields[0], "gene": fields[3], "strand": fields[5]})
    if not rows:
        raise ValueError("Empty BED: {}".format(path))
    return rows


def read_orthogroups(path, required_labels):
    gene_to_family = {}
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not required_labels.issubset(reader.fieldnames):
            raise ValueError("Orthogroups lacks required columns: {}".format(
                ",".join(sorted(required_labels - set(reader.fieldnames or [])))))
        for row in reader:
            family = row["Orthogroup"]
            for label in required_labels:
                for member in [x.strip() for x in (row.get(label) or "").split(",") if x.strip()]:
                    if member in gene_to_family and gene_to_family[member] != family:
                        raise ValueError("Gene occurs in multiple orthogroups: {}".format(member))
                    gene_to_family[member] = family
    return gene_to_family


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", required=True, type=Path)
    parser.add_argument("--orthogroups", required=True, type=Path)
    parser.add_argument("--id-maps", required=True, type=Path,
                        help="Frozen id_maps directory containing pzi_A.id_map.tsv, etc.")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=12)
    args = parser.parse_args()
    if args.threads < 1:
        raise SystemExit("ERROR: --threads must be positive")
    if args.outdir.exists() and any(args.outdir.iterdir()):
        raise SystemExit("ERROR: output directory is not empty: {}".format(args.outdir))
    args.outdir.mkdir(parents=True, exist_ok=True)
    lists_dir = args.outdir / "lists"
    lists_dir.mkdir()

    maps = {}
    stable_to_original = {}
    for track, (label, _species) in TRACKS.items():
        maps[label], stable_to_original[label] = read_id_map(args.id_maps / (label + ".id_map.tsv"))
    families = read_orthogroups(args.orthogroups, {label for label, _ in TRACKS.values()})

    family_rows = []
    for stable_id, family in sorted(families.items()):
        if any(stable_id in stable_to_original[label] for label, _ in TRACKS.values()):
            family_rows.append((stable_id, family))
    with (args.outdir / "families.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        for stable_id, family in family_rows:
            handle.write("{}\t{}\n".format(stable_id, family))

    manifest = []
    stable_map_rows = []
    for track, (label, species) in TRACKS.items():
        bed_rows = read_bed(args.catalog / "tracks" / (track + ".bed"))
        list_by_seqid = OrderedDict()
        for row in bed_rows:
            raw_gene = row["gene"].split("__", 1)[1] if "__" in row["gene"] else row["gene"]
            stable_id = maps[label].get(raw_gene)
            if stable_id is None:
                raise ValueError("Catalog gene absent from {} id map: {}".format(label, raw_gene))
            if stable_id not in families:
                continue
            list_by_seqid.setdefault(row["seqid"].split("__", 1)[1], []).append((stable_id, row["strand"]))
            stable_map_rows.append({
                "track": track, "species": species, "stable_id": stable_id,
                "original_id": stable_to_original[label][stable_id],
            })
        for seqid, genes in list_by_seqid.items():
            list_name = "{}_{}".format("A" if track.endswith("_A") else "B", seqid)
            species_lists = lists_dir / species
            species_lists.mkdir(exist_ok=True)
            path = species_lists / (list_name + ".lst")
            with path.open("w", encoding="utf-8", newline="\n") as handle:
                for stable_id, strand in genes:
                    handle.write("{}{}\n".format(stable_id, strand if strand in ("+", "-") else "+"))
            manifest.append({"species": species, "track": track, "list_name": list_name,
                             "seqid": seqid, "path": str(path.resolve()), "genes": len(genes)})

    with (args.outdir / "iadhore.ini").open("w", encoding="utf-8", newline="\n") as handle:
        for species in ("PZI", "PSP", "PPA"):
            handle.write("genome= {}\n".format(species))
            for row in manifest:
                if row["species"] == species:
                    handle.write("{} {}\n".format(row["list_name"], row["path"]))
            handle.write("\n")
        handle.write("blast_table= {}\n".format((args.outdir / "families.tsv").resolve()))
        handle.write("table_type= family\n")
        handle.write("output_path= {}\n".format((args.outdir / "iadhore_output").resolve()))
        handle.write("gap_size= 30\ncluster_gap= 35\nq_value= 0.75\nprob_cutoff= 0.01\n")
        handle.write("anchor_points= 3\nalignment_method= gg2\nlevel_2_only= false\n")
        handle.write("multiple_hypothesis_correction= FDR\nvisualizeGHM= false\nnumber_of_threads= {}\n".format(args.threads))

    fields = ["species", "track", "list_name", "seqid", "path", "genes"]
    with (args.outdir / "list_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(manifest)
    fields = ["track", "species", "stable_id", "original_id"]
    with (args.outdir / "stable_gene_map.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(stable_map_rows)
    (args.outdir / "method_contract.txt").write_text(
        "Reference Fig. 3b collinearity input preparation. Orthogroup family assignments are used as i-ADHoRe family labels; i-ADHoRe 3.0 performs the gene-order homology search.\n"
        "The output anchor pairs are the sole collinearity input to the subsequent family screen. No JCVI/DIAMOND anchor table is read.\n",
        encoding="ascii")
    print("[OK] Prepared i-ADHoRe input: {} lists, {} family memberships".format(len(manifest), len(family_rows)))


if __name__ == "__main__":
    main()
