#!/usr/bin/env python3
"""Prepare a species-level 4D alignment for an MCMCtree fossil-calibrated run.

The input protein alignments are the existing trimmed single-copy alignments
with taxon labels.  Orthogroup rows and frozen ID maps recover the original
gene IDs needed to extract the matching CDS.  PZI/PSP/PPA A and B copies are
not treated as separate species here; only their A representative is used.
"""

import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path


FOURFOLD_PREFIXES = {"GC", "GG", "CC", "AC", "GT", "TC", "CG", "CT"}
VALID_BASES = set("ACGT")


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path):
    order = []
    records = {}
    current = None
    chunks = []
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if current is not None:
                    records[current] = "".join(chunks).upper()
                current = line[1:].strip().split()[0]
                order.append(current)
                chunks = []
            else:
                chunks.append(line.strip())
        if current is not None:
            records[current] = "".join(chunks).upper()
    return order, records


def write_fasta(path, order, records):
    with Path(path).open("w", encoding="utf-8") as handle:
        for label in order:
            handle.write(">%s\n" % label)
            seq = records[label]
            for start in range(0, len(seq), 80):
                handle.write(seq[start:start + 80] + "\n")


def parse_cell(value):
    return [x for x in re.split(r"[ ,]+", value.strip()) if x]


def read_id_maps(map_dir, species):
    maps = {}
    for sp in species:
        path = Path(map_dir) / (sp + ".id_map.tsv")
        if not path.is_file():
            raise SystemExit("Missing ID map: %s" % path)
        result = {}
        for row in read_tsv(path):
            stable = row.get("stable_id", "")
            original = row.get("original_id", "")
            if stable and original:
                result[stable] = original
        maps[sp] = result
    return maps


def index_fasta(path):
    _, records = read_fasta(path)
    index = {}
    for key, seq in records.items():
        candidates = {key}
        candidates.add(key.split()[0])
        candidates.add(key.rsplit("|", 1)[-1])
        for candidate in candidates:
            index.setdefault(candidate, seq)
    return index


def find_cds(index, original_id):
    candidates = [original_id, original_id.split()[0]]
    if "|" in original_id:
        candidates.append(original_id.split("|", 1)[0])
    if original_id.endswith(".p"):
        candidates.append(original_id[:-2])
    if ".p" in original_id:
        candidates.append(original_id.rsplit(".p", 1)[0])
    for candidate in candidates:
        if candidate in index:
            return index[candidate]
    return None


def valid_cds(seq):
    seq = seq.upper().replace("-", "")
    return bool(seq) and len(seq) % 3 == 0 and re.fullmatch(r"[ACGTN]+", seq) is not None


def extract_4d(records, order):
    lengths = {len(records[x]) for x in order}
    if len(lengths) != 1:
        return None, 0, "unequal_codon_alignment_lengths"
    length = next(iter(lengths))
    if length == 0 or length % 3:
        return None, 0, "codon_alignment_not_multiple_of_three"
    result = {x: [] for x in order}
    count = 0
    for pos in range(0, length, 3):
        codons = [records[x][pos:pos + 3] for x in order]
        if any(len(c) != 3 or c[:2] not in FOURFOLD_PREFIXES or any(b not in VALID_BASES for b in c) for c in codons):
            continue
        count += 1
        for x in order:
            result[x].append(records[x][pos + 2])
    return {x: "".join(v) for x, v in result.items()}, count, "ok"


def write_phylip(path, order, records):
    length = len(records[order[0]]) if order else 0
    with Path(path).open("w", encoding="utf-8") as handle:
        handle.write("%d %d\n" % (len(order), length))
        for label in order:
            handle.write("%-12s %s\n" % (label, records[label]))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trimmed-dir", required=True)
    parser.add_argument("--orthogroups", required=True)
    parser.add_argument("--single-copy", required=True)
    parser.add_argument("--id-map-dir", required=True)
    parser.add_argument("--cds-map", action="append", required=True, help="species=path; repeat for every species")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--pal2nal", default="pal2nal.pl")
    parser.add_argument("--min-4d-sites", type=int, default=12)
    args = parser.parse_args()

    species = ["gtr", "rch", "ppa", "psp", "pzi", "pmal", "pav", "pspe", "pcam", "psal", "par", "pmu", "pdu", "ppe"]
    original_labels = {"ppa": "ppa_A", "psp": "psp_A", "pzi": "pzi_A"}
    of_species = ["gtr", "rch", "ppa_A", "psp_A", "pzi_A", "pmal", "pav", "pspe", "pcam", "psal", "par", "pmu", "pdu", "ppe"]
    out = Path(args.out_dir)
    for sub in ("codon_alignments", "aa_alignments", "cds_inputs"):
        (out / sub).mkdir(parents=True, exist_ok=True)

    cds_paths = {}
    for item in args.cds_map:
        if "=" not in item:
            raise SystemExit("Invalid --cds-map (expected species=path): %s" % item)
        key, value = item.split("=", 1)
        cds_paths[key] = Path(value)
    for sp in of_species:
        label = original_labels.get(sp, sp)
        if label not in cds_paths and sp not in cds_paths:
            raise SystemExit("Missing CDS map for %s" % sp)

    maps = read_id_maps(args.id_map_dir, of_species)
    cds_indexes = {}
    for sp in of_species:
        label = original_labels.get(sp, sp)
        path = cds_paths.get(label, cds_paths.get(sp))
        if not path.is_file():
            raise SystemExit("Missing CDS FASTA for %s: %s" % (sp, path))
        cds_indexes[sp] = index_fasta(path)

    rows = {row["Orthogroup"]: row for row in read_tsv(args.orthogroups)}
    selected = [x.strip() for x in Path(args.single_copy).read_text(encoding="utf-8").splitlines() if x.strip()]
    trimmed_dir = Path(args.trimmed_dir)
    stats = []
    mapping_rows = []
    used = 0
    for og in selected:
        row = rows.get(og)
        aa_path = trimmed_dir / (og + ".fa")
        if row is None or not aa_path.is_file():
            stats.append((og, "missing_orthogroup_or_alignment", 0, 0))
            continue
        aa_order, aa_records_all = read_fasta(aa_path)
        if any(label not in aa_records_all for label in of_species):
            stats.append((og, "missing_taxon_in_protein_alignment", 0, 0))
            continue
        aa_records = {}
        cds_records = {}
        failed = None
        for sp in of_species:
            stable_values = parse_cell(row.get(sp, ""))
            if len(stable_values) != 1:
                failed = "not_single_copy_%s" % sp
                break
            stable = stable_values[0]
            original = maps[sp].get(stable)
            if not original:
                failed = "stable_id_not_mapped_%s" % sp
                break
            cds = find_cds(cds_indexes[sp], original)
            if cds is None or not valid_cds(cds):
                failed = "cds_not_found_or_invalid_%s" % sp
                break
            aa_records[sp] = aa_records_all[sp]
            cds_records[sp] = cds.rstrip("*")
            mapping_rows.append((og, sp, stable, original))
        if failed:
            stats.append((og, failed, 0, 0))
            continue
        aa_out = out / "aa_alignments" / (og + ".fa")
        cds_in = out / "cds_inputs" / (og + ".fa")
        codon_out = out / "codon_alignments" / (og + ".fa")
        write_fasta(aa_out, of_species, aa_records)
        write_fasta(cds_in, of_species, cds_records)
        cmd = [args.pal2nal, str(aa_out), str(cds_in), "-output", "fasta", "-nogap"]
        proc = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if proc.returncode != 0 or not proc.stdout.startswith(">"):
            stats.append((og, "pal2nal_failed", 0, 0))
            (out / "pal2nal_errors.log").open("a", encoding="utf-8").write("%s\n%s\n" % (og, proc.stderr[-1000:]))
            continue
        codon_out.write_text(proc.stdout, encoding="utf-8")
        _, codon_records = read_fasta(codon_out)
        lengths = {len(codon_records.get(x, "")) for x in of_species}
        if len(lengths) != 1 or not lengths or next(iter(lengths)) % 3:
            stats.append((og, "invalid_codon_alignment", 0, 0))
            codon_out.unlink(missing_ok=True)
            continue
        fourd, n4d, status = extract_4d(codon_records, of_species)
        if status != "ok" or n4d < args.min_4d_sites:
            stats.append((og, "skip_%s" % (status if status != "ok" else "too_few_4d_sites"), len(codon_records[of_species[0]]) // 3, n4d))
            continue
        stats.append((og, "used", len(codon_records[of_species[0]]) // 3, n4d))
        used += 1

    concat = {sp: [] for sp in species}
    used_ogs = []
    for og, status, codons, n4d in stats:
        if status != "used":
            continue
        _, codon_records = read_fasta(out / "codon_alignments" / (og + ".fa"))
        fourd, _, _ = extract_4d(codon_records, of_species)
        for sp in of_species:
            label = original_labels.get(sp, sp)
            concat[label].append(fourd[sp])
        used_ogs.append(og)
    concat_records = {sp: "".join(parts) for sp, parts in concat.items()}
    write_fasta(out / "single_copy_4d_sites.fa", species, concat_records)
    write_phylip(out / "single_copy_4d_sites.phy", species, concat_records)
    with (out / "single_copy_family_stats.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["orthogroup", "status", "codons", "fourfold_sites"])
        writer.writerows(stats)
    with (out / "single_copy_mapping.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["orthogroup", "orthofinder_species", "stable_id", "original_id"])
        writer.writerows(mapping_rows)
    (out / "used_single_copy_orthogroups.txt").write_text("".join(x + "\n" for x in used_ogs), encoding="utf-8")
    (out / "species_order.txt").write_text("".join(x + "\n" for x in species), encoding="utf-8")
    if used < 100:
        raise SystemExit("Too few usable 4D families: %d" % used)
    print("used_families=%d taxa=%d fourfold_sites=%d" % (used, len(species), len(concat_records[species[0]])))


if __name__ == "__main__":
    main()
