#!/usr/bin/env python3
"""Prepare the exact-species, plant-outgroup protein alignment for MCMCtree.

The source alignments contain A/B representatives for PPA, PSP and PZI.
This collapses those homeolog copies to one A representative per species and
retains the nine diploid Prunus taxa plus Gillenia and Rosa as plant outgroups.
"""

import argparse
import csv
import re
from pathlib import Path

SOURCE_TIPS = [
    "ppa_A", "psp_A", "pzi_A",
    "pmal", "pav", "pspe", "pcam", "psal", "par", "pmu", "pdu", "ppe",
    "gtr", "rch",
]
SPECIES_TIPS = [
    "ppa", "psp", "pzi",
    "pmal", "pav", "pspe", "pcam", "psal", "par", "pmu", "pdu", "ppe",
    "gtr", "rch",
]
RENAME = {"ppa_A": "ppa", "psp_A": "psp", "pzi_A": "pzi"}


def read_fasta(path):
    order, records = [], {}
    current, chunks = None, []
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if current is not None:
                    records[current] = "".join(chunks)
                current = line[1:].strip().split()[0]
                order.append(current)
                chunks = []
            else:
                chunks.append(line.strip())
        if current is not None:
            records[current] = "".join(chunks)
    return order, records


def write_fasta(path, order, records):
    with Path(path).open("w", encoding="utf-8") as handle:
        for label in order:
            handle.write(">%s\n" % label)
            sequence = records[label]
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def write_phylip(path, order, records):
    length = len(records[order[0]])
    with Path(path).open("w", encoding="utf-8") as handle:
        handle.write("%d %d\n" % (len(order), length))
        for label in order:
            handle.write("%-10s %s\n" % (label, records[label]))


def allowed_families(partition_file):
    text = Path(partition_file).read_text(encoding="utf-8", errors="replace")
    names = re.findall(r"^\s*charset\s+(\S+)\s*=", text, flags=re.MULTILINE | re.IGNORECASE)
    if not names:
        raise SystemExit("No charset entries found: %s" % partition_file)
    return set(names)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trimmed-dir", required=True)
    parser.add_argument("--partition-file", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--min-families", type=int, default=100)
    args = parser.parse_args()

    selected = allowed_families(args.partition_file)
    out = Path(args.out_dir)
    family_dir = out / "protein_alignments_species14"
    family_dir.mkdir(parents=True, exist_ok=True)
    concat = {tip: [] for tip in SPECIES_TIPS}
    stats = []
    used = 0

    for path in sorted(Path(args.trimmed_dir).glob("*.fa")):
        og = path.stem
        if og not in selected:
            continue
        _, records = read_fasta(path)
        missing = [tip for tip in SOURCE_TIPS if tip not in records]
        if missing:
            stats.append((og, "missing:" + ",".join(missing), 0))
            continue
        lengths = {len(records[tip]) for tip in SOURCE_TIPS}
        if len(lengths) != 1:
            stats.append((og, "unequal_alignment_lengths", 0))
            continue
        species_records = {RENAME.get(tip, tip): records[tip] for tip in SOURCE_TIPS}
        write_fasta(family_dir / (og + ".fa"), SPECIES_TIPS, species_records)
        length = next(iter(lengths))
        for tip in SPECIES_TIPS:
            concat[tip].append(species_records[tip])
        stats.append((og, "used", length))
        used += 1

    if used < args.min_families:
        raise SystemExit("Too few usable families: %d" % used)
    concatenated = {tip: "".join(parts) for tip, parts in concat.items()}
    write_fasta(out / "single_copy_protein_species14.fa", SPECIES_TIPS, concatenated)
    write_phylip(out / "single_copy_protein_species14.phy", SPECIES_TIPS, concatenated)
    with (out / "family_stats.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["orthogroup", "status", "alignment_length"])
        writer.writerows(stats)
    (out / "species_order.txt").write_text("".join(tip + "\n" for tip in SPECIES_TIPS), encoding="utf-8")
    with (out / "species_tip_policy.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["source_tip", "dating_tip", "policy"])
        for source in SOURCE_TIPS:
            policy = "A_representative_of_allopolyploid" if source.endswith("_A") else (
                "plant_outgroup" if source in {"gtr", "rch"} else "diploid_species"
            )
            writer.writerow([source, RENAME.get(source, source), policy])
    print("used_families=%d taxa=%d protein_sites=%d" % (used, len(SPECIES_TIPS), len(concatenated[SPECIES_TIPS[0]])))


if __name__ == "__main__":
    main()
