#!/usr/bin/env python3
"""Filter existing JCVI inputs to chromosome-scale external references.

The script preserves source gene IDs and only removes non-chromosomal
sequences.  It deliberately does not rename any PZI sequence to an A/B label.
"""

from __future__ import print_function

import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", action="append", nargs=3, required=True,
                        metavar=("LABEL", "INPUT_DIR", "SEQIDS"),
                        help="Repeat: label, existing JCVI input directory, comma-separated chromosome seqids")
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def read_fasta(path):
    header = None
    seq = []
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header = line[1:].rstrip("\n")
                seq = []
            else:
                seq.append(line.strip())
        if header is not None:
            yield header, "".join(seq)


def write_fasta(path, records):
    with Path(path).open("w", encoding="ascii", newline="\n") as handle:
        for header, seq in records:
            handle.write(">{}\n".format(header))
            for start in range(0, len(seq), 80):
                handle.write(seq[start:start + 80] + "\n")


def read_bed(path, allowed):
    rows = []
    ids = set()
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4 or fields[0] not in allowed:
                continue
            rows.append(fields)
            ids.add(fields[3])
    return rows, ids


def resolve_input(source_dir, label, suffix):
    candidate = Path(source_dir) / (label + suffix)
    if not candidate.is_file():
        raise SystemExit("Missing input: {}".format(candidate))
    return candidate


def main():
    args = parse_args()
    outroot = Path(args.outdir)
    outroot.mkdir(parents=True, exist_ok=True)
    summary = []
    for label, source_dir, sequence_text in args.source:
        allowed = [value for value in sequence_text.split(",") if value]
        if not allowed:
            raise SystemExit("No chromosome seqids supplied for {}".format(label))
        source = Path(source_dir)
        bed_in = resolve_input(source, label, ".bed")
        pep_in = resolve_input(source, label, ".pep")
        cds_in = resolve_input(source, label, ".cds")
        bed_rows, gene_ids = read_bed(bed_in, set(allowed))
        if not bed_rows:
            raise SystemExit("No retained BED records for {}".format(label))
        pep = [(header, seq) for header, seq in read_fasta(pep_in) if header.split()[0] in gene_ids]
        cds = [(header, seq) for header, seq in read_fasta(cds_in) if header.split()[0] in gene_ids]
        pep_ids = set(header.split()[0] for header, _seq in pep)
        cds_ids = set(header.split()[0] for header, _seq in cds)
        missing_pep = sorted(gene_ids - pep_ids)
        missing_cds = sorted(gene_ids - cds_ids)
        if missing_pep or missing_cds:
            raise SystemExit(
                "{} has BED records without peptide/CDS: pep_missing={}, cds_missing={}".format(
                    label, len(missing_pep), len(missing_cds)))
        target = outroot / label
        target.mkdir(parents=True, exist_ok=True)
        with (target / (label + ".bed")).open("w", encoding="ascii", newline="\n") as handle:
            for row in bed_rows:
                handle.write("\t".join(row) + "\n")
        write_fasta(target / (label + ".pep"), pep)
        write_fasta(target / (label + ".cds"), cds)
        (target / (label + ".seqids")).write_text(",".join(allowed) + "\n", encoding="ascii")
        summary.append({
            "label": label,
            "source_dir": str(source),
            "retained_seqids": ",".join(allowed),
            "bed_genes": len(gene_ids),
            "pep_records": len(pep),
            "cds_records": len(cds),
            "missing_pep": len(missing_pep),
            "missing_cds": len(missing_cds),
        })
    with (outroot / "external_input_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0].keys()), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary)
    print("[OK] Prepared chromosome-scale external inputs: {}".format(outroot))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
