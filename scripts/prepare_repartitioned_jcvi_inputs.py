#!/usr/bin/env python3
"""Stage BED/protein inputs for the revised PZI/PPA/PSP JCVI run.

The revised haplotype package is the authoritative source.  Do not fall back
to the legacy annotation/repartition tree: those files were generated before
the chromosome reassignment and would silently mix coordinate systems.
"""

from __future__ import annotations

import argparse
import os
import shutil
from pathlib import Path


LABELS = {
    "pzi": ("A1", "A2", "B1", "B2"),
    "psp": ("A1", "A2", "B1", "B2"),
    "ppa": ("A", "B"),
}


def attrs(text: str) -> dict[str, str]:
    out = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = value
        elif " " in item:
            key, value = item.split(None, 1)
            out[key] = value.strip('"')
    return out


def stage_one(root: Path, out: Path, species: str, label: str) -> tuple[Path, Path, int]:
    # These GFF3/pep files were extracted after the final haplotype
    # reassignment and use the same chromosome vocabulary as the revised FASTA.
    input_root = root / "haploid_repartition" / "inputs" / species
    gff = input_root / "gff" / f"{species}_hap{label}.gff3"
    pep = input_root / "pep" / f"{species}_hap{label}.pep.fa"
    if not gff.is_file() or gff.stat().st_size == 0:
        raise SystemExit(f"Missing postprocessed GFF3: {gff}")
    if not pep.is_file() or pep.stat().st_size == 0:
        raise SystemExit(f"Missing postprocessed protein FASTA: {pep}")
    bed = out / f"{species}_{label}.bed"
    records = []
    seen = set()
    with gff.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"mRNA", "transcript"}:
                continue
            a = attrs(fields[8])
            feature_id = a.get("ID") or a.get("transcript_id") or a.get("Name")
            if not feature_id or feature_id in seen:
                continue
            seen.add(feature_id)
            records.append((fields[0], int(fields[3]) - 1, int(fields[4]), feature_id, fields[6]))
    records.sort(key=lambda row: (row[0], row[1], row[2], row[3]))
    if not records:
        raise SystemExit(f"No mRNA/transcript records found in {gff}")
    with bed.open("w", encoding="utf-8") as handle:
        for chrom, start, end, feature_id, strand in records:
            handle.write(f"{chrom}\t{start}\t{end}\t{feature_id}\t0\t{strand}\n")
    pep_link = out / f"{species}_{label}.pep"
    if pep_link.exists() or pep_link.is_symlink():
        pep_link.unlink()
    try:
        os.symlink(pep.resolve(), pep_link)
    except OSError:
        shutil.copy2(pep, pep_link)
    return bed, pep_link, len(records)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out-root", type=Path, required=True)
    args = parser.parse_args()
    input_dir = args.out_root / "work" / "jcvi_inputs"
    input_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for species, labels in LABELS.items():
        for label in labels:
            bed, pep, count = stage_one(args.root, input_dir, species, label)
            rows.append({"name": f"{species}_{label}", "bed": str(bed), "protein": str(pep), "genes": count})
    with (input_dir / "manifest.tsv").open("w", encoding="utf-8") as handle:
        handle.write("name\tbed\tprotein\tgenes\n")
        for row in rows:
            handle.write("{name}\t{bed}\t{protein}\t{genes}\n".format(**row))
    (input_dir / ".staged.ok").write_text("repartitioned JCVI inputs staged\n", encoding="ascii")
    print(f"[OK] Staged {len(rows)} revised JCVI inputs in {input_dir}")


if __name__ == "__main__":
    main()
