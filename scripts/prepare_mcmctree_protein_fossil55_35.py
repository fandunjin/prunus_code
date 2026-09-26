#!/usr/bin/env python3
"""Build a Prunus species-level concatenated protein alignment for MCMCtree.

The existing trimmed alignments contain 17 tips (A/B copies for PZI, PSP and
PPA). This keeps one A representative for each species, removes B copies, and
omits the Gillenia and Rosa outgroups because the calibrated Prunus crown is
the root of the dating tree.
"""

import argparse
import csv
import re
from pathlib import Path


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
            seq = records[label]
            for start in range(0, len(seq), 80):
                handle.write(seq[start:start + 80] + "\n")


def write_phylip(path, order, records):
    length = len(records[order[0]]) if order else 0
    with Path(path).open("w", encoding="utf-8") as handle:
        handle.write("%d %d\n" % (len(order), length))
        for label in order:
            handle.write("%-10s %s\n" % (label, records[label]))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trimmed-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--partition-file", default="", help="Existing NEXUS partition file; only listed families are used")
    parser.add_argument("--min-families", type=int, default=100)
    args = parser.parse_args()

    source_order = [
        "ppa_A", "psp_A", "pzi_A", "pmal", "pav", "pspe", "pcam",
        "psal", "par", "pmu", "pdu", "ppe",
    ]
    species_order = [
        "ppa", "psp", "pzi", "pmal", "pav", "pspe", "pcam", "psal",
        "par", "pmu", "pdu", "ppe",
    ]
    rename = {"ppa_A": "ppa", "psp_A": "psp", "pzi_A": "pzi"}
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    family_dir = out / "protein_alignments_species_level"
    family_dir.mkdir(exist_ok=True)

    allowed = None
    if args.partition_file:
        partition_path = Path(args.partition_file)
        if not partition_path.is_file():
            raise SystemExit("Missing partition file: %s" % partition_path)
        allowed = []
        text = partition_path.read_text(encoding="utf-8", errors="replace")
        for match in re.finditer(r"^\s*charset\s+(\S+)\s*=", text, flags=re.MULTILINE | re.IGNORECASE):
            allowed.append(match.group(1))
        if not allowed:
            raise SystemExit("No charset entries found in partition file: %s" % partition_path)
        allowed = set(allowed)

    concat = {label: [] for label in species_order}
    stats = []
    used = 0
    for path in sorted(Path(args.trimmed_dir).glob("*.fa")):
        og = path.stem
        if allowed is not None and og not in allowed:
            continue
        _, records = read_fasta(path)
        if any(label not in records for label in source_order):
            stats.append((og, "missing_required_tip", 0))
            continue
        lengths = {len(records[label]) for label in source_order}
        if len(lengths) != 1:
            stats.append((og, "unequal_alignment_lengths", 0))
            continue
        species_records = {rename.get(label, label): records[label] for label in source_order}
        write_fasta(family_dir / (og + ".fa"), species_order, species_records)
        for label in species_order:
            concat[label].append(species_records[label])
        used += 1
        stats.append((og, "used", next(iter(lengths))))

    concat_records = {label: "".join(parts) for label, parts in concat.items()}
    write_fasta(out / "single_copy_protein.fa", species_order, concat_records)
    write_phylip(out / "single_copy_protein.phy", species_order, concat_records)
    with (out / "family_stats.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["orthogroup", "status", "alignment_length"])
        writer.writerows(stats)
    (out / "species_order.txt").write_text("".join(x + "\n" for x in species_order), encoding="utf-8")
    (out / "species_tip_policy.tsv").write_text(
        "source_tip\tdating_tip\tpolicy\n" +
        "".join("%s\t%s\t%s\n" % (src, rename.get(src, src), "A_representative" if src.endswith("_A") else "species_tip") for src in source_order),
        encoding="utf-8",
    )
    if used < args.min_families:
        raise SystemExit("Too few usable protein families: %d" % used)
    print("used_families=%d taxa=%d protein_sites=%d" % (used, len(species_order), len(concat_records[species_order[0]])))


if __name__ == "__main__":
    main()
