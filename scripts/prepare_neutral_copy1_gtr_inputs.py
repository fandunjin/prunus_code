#!/usr/bin/env python3
"""Prepare neutral PZI copy1 and chromosome-only GTR JCVI inputs.

This helper deliberately keeps the PZI labels as ``chr01``--``chr16`` and
does not infer or introduce A/B subgenome labels.  GTR is restricted to its
nine assembled chromosomes (Chr01--Chr09); unplaced sequences are excluded
from the model comparison.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gtr-input", required=True, help="Existing JCVI GTR input directory")
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def fasta_records(path: Path):
    name = None
    seq = []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"): 
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].split()[0]
                seq = []
            elif name is not None:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def write_fasta(path: Path, records):
    with path.open("w", encoding="utf-8") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def main() -> int:
    a = args()
    src = Path(a.gtr_input)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    for required in ("gtr.bed", "gtr.pep", "gtr.cds"):
        path = src / required
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Missing GTR JCVI input: {path}")

    keep = {f"Chr{i:02d}" for i in range(1, 10)}
    ids = []
    bed_lines = []
    with (src / "gtr.bed").open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4 or fields[0] not in keep:
                continue
            bed_lines.append(line.rstrip("\n"))
            ids.append(fields[3])
    if not bed_lines:
        raise SystemExit("No GTR Chr01--Chr09 genes were found in gtr.bed")
    idset = set(ids)
    (out / "gtr.bed").write_text("\n".join(bed_lines) + "\n", encoding="utf-8")
    for ext in ("pep", "cds"):
        selected = [(name, seq) for name, seq in fasta_records(src / f"gtr.{ext}") if name in idset]
        write_fasta(out / f"gtr.{ext}", selected)
        if len(selected) != len(idset) and ext == "pep":
            missing = len(idset) - len(selected)
            raise SystemExit(f"GTR protein records missing for {missing} bed genes")
    (out / "gtr.seqids").write_text(",".join(sorted(keep, key=lambda x: int(x[3:]))) + "\n", encoding="utf-8")
    with (out / "filter_summary.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["field", "value"])
        w.writerow(["gtr_chromosomes", ",".join(sorted(keep, key=lambda x: int(x[3:])))])
        w.writerow(["gtr_bed_records", len(bed_lines)])
        w.writerow(["gtr_protein_records", sum(1 for _ in fasta_records(out / "gtr.pep"))])
        w.writerow(["gtr_cds_records", sum(1 for _ in fasta_records(out / "gtr.cds"))])
    print(f"[OK] Neutral GTR inputs: {out} ({len(bed_lines)} genes; 9 chromosomes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
