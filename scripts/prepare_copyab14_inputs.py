#!/usr/bin/env python3
"""Prepare a 14-terminal OrthoFinder protein set with A1+B1 copyAB taxa."""

from __future__ import annotations

import argparse
from pathlib import Path


def rewrite_fasta(src: Path, dst: Path, prefix: str) -> int:
    count = 0
    with src.open("r", encoding="utf-8", errors="replace") as inp, dst.open("a", encoding="utf-8") as out:
        for line in inp:
            if line.startswith(">"):
                token = line[1:].strip().split()[0]
                rest = line[1:].strip()[len(token) :]
                out.write(f">{prefix}__{token}{rest}\n")
                count += 1
            else:
                out.write(line)
    return count


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", type=Path, required=True)
    ap.add_argument("--base", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--manifest", type=Path, required=True)
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    for old in args.out.glob("*.fa"):
        old.unlink()
    args.manifest.parent.mkdir(parents=True, exist_ok=True)

    # The 11 non-polyploid terminals use the already frozen 0903 proteins.
    ordinary = ["gtr", "par", "pav", "pcam", "pdu", "pmal", "pmu", "ppe", "psal", "pspe", "rch"]
    records: list[tuple[str, str, int]] = []
    for label in ordinary:
        src = args.base / f"{label}.fa"
        if not src.is_file():
            raise SystemExit(f"missing ordinary input: {src}")
        dst = args.out / f"{label}.fa"
        n = rewrite_fasta(src, dst, label)
        records.append((label, str(src), n))

    copyab = {
        "pzi_copyAB": [args.source / "pzi" / "pep" / "pzi_hapA1.pep.fa", args.source / "pzi" / "pep" / "pzi_hapB1.pep.fa"],
        "psp_copyAB": [args.source / "psp" / "pep" / "psp_hapA1.pep.fa", args.source / "psp" / "pep" / "psp_hapB1.pep.fa"],
        "ppa_copyAB": [args.source / "ppa" / "pep" / "ppa_hapA.pep.fa", args.source / "ppa" / "pep" / "ppa_hapB.pep.fa"],
    }
    for label, sources in copyab.items():
        dst = args.out / f"{label}.fa"
        total = 0
        for src in sources:
            if not src.is_file():
                raise SystemExit(f"missing copyAB input: {src}")
            hap = "A1" if "hapA1" in src.name or "hapA." in src.name else "B1"
            total += rewrite_fasta(src, dst, f"{label}_{hap}")
        records.append((label, ";".join(str(x) for x in sources), total))

    with args.manifest.open("w", encoding="utf-8") as out:
        out.write("label\tsource\tprotein_records\tfile\n")
        for label, source, n in records:
            out.write(f"{label}\t{source}\t{n}\t{args.out / (label + '.fa')}\n")
    print(f"[OK] Prepared {len(records)} OrthoFinder terminals in {args.out}")
    print(f"[OK] Manifest: {args.manifest}")


if __name__ == "__main__":
    main()
