#!/usr/bin/env python3
"""Filter trimmed protein alignments and build a concatenated supermatrix."""

import argparse
import csv
from pathlib import Path


def read_fasta(path):
    records = []
    name = None
    chunks = []
    with path.open("r", encoding="ascii") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.upper())
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def alignment_stats(records):
    if not records:
        return 0, 0, 0, 0
    lengths = set(len(seq) for _, seq in records)
    if len(lengths) != 1:
        raise ValueError("sequences have unequal lengths")
    length = lengths.pop()
    variable = 0
    informative = 0
    for col in range(length):
        states = {}
        for _, seq in records:
            aa = seq[col]
            if aa in "-?X":
                continue
            states[aa] = states.get(aa, 0) + 1
        if len(states) >= 2:
            variable += 1
        if sum(1 for count in states.values() if count >= 2) >= 2:
            informative += 1
    missing = sum(seq.count("-") + seq.count("?") + seq.count("X") for _, seq in records)
    cells = max(1, len(records) * length)
    return length, variable, informative, float(missing) / cells


def write_fasta(path, records):
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for name, sequence in records:
            handle.write(">%s\n%s\n" % (name, sequence))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment-dir", required=True)
    parser.add_argument("--species-file", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--min-length", type=int, default=100)
    parser.add_argument("--min-informative", type=int, default=10)
    parser.add_argument("--max-missing", type=float, default=0.20)
    parser.add_argument("--min-alignments", type=int, default=200)
    args = parser.parse_args()

    source = Path(args.alignment_dir)
    out_dir = Path(args.out_dir)
    retained_dir = out_dir / "retained"
    retained_dir.mkdir(parents=True, exist_ok=False)
    species = [line.strip() for line in Path(args.species_file).read_text(encoding="utf-8").splitlines() if line.strip()]

    rows = []
    retained = []
    for path in sorted(source.glob("*.fa")):
        og = path.stem
        status = "PASS"
        reason = "retained"
        try:
            records = read_fasta(path)
            names = [name for name, _ in records]
            if len(names) != len(set(names)):
                raise ValueError("duplicate taxon")
            unknown = sorted(set(names) - set(species))
            if unknown:
                raise ValueError("unknown taxa: %s" % ",".join(unknown))
            length, variable, informative, missing = alignment_stats(records)
            if length < args.min_length:
                status, reason = "FAIL", "alignment_shorter_than_%d" % args.min_length
            elif informative < args.min_informative:
                status, reason = "FAIL", "fewer_than_%d_informative_sites" % args.min_informative
            elif missing > args.max_missing:
                status, reason = "FAIL", "missing_fraction_above_%.3f" % args.max_missing
        except ValueError as error:
            records = []
            length = variable = informative = 0
            missing = 1.0
            status, reason = "FAIL", str(error).replace("\t", " ")

        rows.append({
            "orthogroup": og,
            "taxa": str(len(records)),
            "aligned_sites": str(length),
            "variable_sites": str(variable),
            "informative_sites": str(informative),
            "missing_fraction": "%.6f" % missing,
            "status": status,
            "reason": reason,
        })
        if status == "PASS":
            destination = retained_dir / path.name
            write_fasta(destination, records)
            retained.append((og, dict(records), length))

    fields = ["orthogroup", "taxa", "aligned_sites", "variable_sites", "informative_sites", "missing_fraction", "status", "reason"]
    with (out_dir / "alignment_qc.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    if len(retained) < args.min_alignments:
        raise SystemExit("ERROR: retained %d alignments; minimum is %d" % (len(retained), args.min_alignments))

    concat = dict((sp, []) for sp in species)
    partitions = []
    start = 1
    for og, records, length in retained:
        for sp in species:
            concat[sp].append(records.get(sp, "-" * length))
        end = start + length - 1
        partitions.append((og, start, end))
        start = end + 1
    write_fasta(out_dir / "concat_protein.fa", [(sp, "".join(concat[sp])) for sp in species])
    with (out_dir / "concat_partitions.nex").open("w", encoding="ascii", newline="\n") as handle:
        handle.write("#nexus\nbegin sets;\n")
        for og, first, last in partitions:
            handle.write("  charset %s = %d-%d;\n" % (og, first, last))
        handle.write("end;\n")
    with (out_dir / "retained_orthogroups.txt").open("w", encoding="ascii", newline="\n") as handle:
        for og, _, _ in retained:
            handle.write(og + "\n")
    (out_dir / "COMPLETE").write_text("Alignment QC and concatenation completed.\n", encoding="ascii")


if __name__ == "__main__":
    main()
