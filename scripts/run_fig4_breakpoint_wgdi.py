#!/usr/bin/env python3
"""Map the nine-unit Fig. 4 reference to resolved Prunus genomes with WGDI.

Run this script in the existing WSL analysis environment. It consumes only
the staged reference and target annotations beneath 0903 and writes raw,
auditable DIAMOND, collinearity and blockinfo outputs for junction discovery.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


def run(command: list[str], log: Path) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("a", encoding="utf-8") as handle:
        handle.write("$ " + " ".join(command) + "\n")
        subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)


def write_conf(path: Path, section: str, values: dict[str, object]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f"[{section}]\n")
        for key, value in values.items():
            handle.write(f"{key} = {value}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True, help="0903 directory")
    parser.add_argument("--wgdi", required=True)
    parser.add_argument("--diamond", required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--targets", default="", help="Comma-separated target labels; default maps every staged target")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    root = args.root.resolve()
    work = root / "karyotype" / "fig4_breakpoint_candidates"
    reference = work / "reference"
    reference_gff = reference / "prunus9_reference.gff"
    reference_lens = reference / "prunus9_reference.lens"
    reference_pep = reference / "prunus9_reference.pep.fa"
    target_manifest = work / "target_manifest.tsv"
    for path in (reference_gff, reference_lens, reference_pep, target_manifest):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    output = work / "mapping"
    logs = work / "logs"
    output.mkdir(parents=True, exist_ok=True)

    with target_manifest.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    requested = {item.strip() for item in args.targets.split(",") if item.strip()}
    if requested:
        missing = requested - {row["target_label"] for row in rows}
        if missing:
            raise ValueError(f"Unknown targets: {', '.join(sorted(missing))}")
        rows = [row for row in rows if row["target_label"] in requested]

    for row in rows:
        label = row["target_label"]
        target_gff = root / row["gff"]
        target_lens = root / row["lens"]
        target_pep = root / row["pep"]
        for path in (target_gff, target_lens, target_pep):
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(path)
        pair_dir = output / label
        pair_dir.mkdir(parents=True, exist_ok=True)
        database = pair_dir / f"{label}.dmnd"
        blast = pair_dir / f"prunus9_vs_{label}.blast.tsv"
        collinearity = pair_dir / f"prunus9_vs_{label}.collinearity"
        blockinfo = pair_dir / f"prunus9_vs_{label}.blockinfo.csv"
        if args.force or not database.is_file():
            run([args.diamond, "makedb", "--in", str(target_pep), "--db", str(database.with_suffix(""))], logs / f"{label}.makedb.log")
        if args.force or not blast.is_file():
            run([
                args.diamond, "blastp", "--query", str(reference_pep), "--db", str(database), "--out", str(blast),
                "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                "--evalue", "1e-5", "--max-target-seqs", "20", "--ultra-sensitive", "--threads", str(args.threads),
            ], logs / f"{label}.blast.log")
        if args.force or not collinearity.is_file():
            config = pair_dir / "icl.conf"
            write_conf(config, "collinearity", {
                "gff1": reference_gff, "gff2": target_gff, "lens1": reference_lens, "lens2": target_lens,
                "blast": blast, "blast_reverse": "false", "multiple": 1, "process": args.threads,
                "evalue": "1e-5", "score": 100, "grading": "50,40,25", "mg": "25,25", "pvalue": 1,
                "repeat_number": 10, "position": "order", "savefile": collinearity,
            })
            run([args.wgdi, "-icl", str(config)], logs / f"{label}.icl.log")
        if args.force or not blockinfo.is_file():
            config = pair_dir / "blockinfo.conf"
            write_conf(config, "blockinfo", {
                "blast": blast, "gff1": reference_gff, "gff2": target_gff, "lens1": reference_lens,
                "lens2": target_lens, "collinearity": collinearity, "score": 100, "evalue": "1e-5",
                "repeat_number": 10, "position": "order", "ks": "none", "ks_col": "ks_NG86", "savefile": blockinfo,
            })
            run([args.wgdi, "-bi", str(config)], logs / f"{label}.blockinfo.log")
        print(f"[OK] {label}: {blockinfo}")


if __name__ == "__main__":
    main()
