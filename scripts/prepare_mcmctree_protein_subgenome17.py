#!/usr/bin/env python3
"""Prepare the 17-tip subgenome-resolved protein alignment for MCMCtree."""

import argparse
import csv
import re
from pathlib import Path


TIP_ORDER = [
    "ppa_A", "psp_A", "pzi_A", "ppa_B", "psp_B", "pzi_B",
    "pmal", "pav", "pspe", "pcam", "psal", "par", "pmu", "pdu",
    "ppe", "gtr", "rch",
]


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
            handle.write("%-12s %s\n" % (label, records[label]))


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
    args = parser.parse_args()

    selected = allowed_families(args.partition_file)
    out = Path(args.out_dir)
    family_dir = out / "protein_alignments_subgenome17"
    family_dir.mkdir(parents=True, exist_ok=True)
    concat = {tip: [] for tip in TIP_ORDER}
    stats = []
    used = 0

    for path in sorted(Path(args.trimmed_dir).glob("*.fa")):
        og = path.stem
        if og not in selected:
            continue
        _, records = read_fasta(path)
        missing = [tip for tip in TIP_ORDER if tip not in records]
        if missing:
            stats.append((og, "missing:" + ",".join(missing), 0))
            continue
        lengths = {len(records[tip]) for tip in TIP_ORDER}
        if len(lengths) != 1:
            stats.append((og, "unequal_alignment_lengths", 0))
            continue
        write_fasta(family_dir / (og + ".fa"), TIP_ORDER, {tip: records[tip] for tip in TIP_ORDER})
        length = next(iter(lengths))
        for tip in TIP_ORDER:
            concat[tip].append(records[tip])
        stats.append((og, "used", length))
        used += 1

    if used < 100:
        raise SystemExit("Too few usable families: %d" % used)
    concatenated = {tip: "".join(parts) for tip, parts in concat.items()}
    write_fasta(out / "single_copy_protein_subgenome17.fa", TIP_ORDER, concatenated)
    write_phylip(out / "single_copy_protein_subgenome17.phy", TIP_ORDER, concatenated)
    with (out / "family_stats.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["orthogroup", "status", "alignment_length"])
        writer.writerows(stats)
    (out / "species_order.txt").write_text("".join(tip + "\n" for tip in TIP_ORDER), encoding="utf-8")
    with (out / "species_tip_policy.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["tip", "role", "note"])
        for tip in TIP_ORDER:
            role = "subgenome_copy" if tip.endswith("_A") or tip.endswith("_B") else ("plant_outgroup" if tip in {"gtr", "rch"} else "diploid_species")
            writer.writerow([tip, role, "A/B copies are lineage tips, not independent species" if role == "subgenome_copy" else ""])
    print("used_families=%d taxa=%d protein_sites=%d" % (used, len(TIP_ORDER), len(concatenated[TIP_ORDER[0]])))


if __name__ == "__main__":
    main()
