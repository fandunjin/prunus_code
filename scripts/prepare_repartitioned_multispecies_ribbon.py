#!/usr/bin/env python3
"""Prepare a native JCVI three-row ribbon for revised PZI, PSP and PPA."""

from __future__ import annotations

import argparse
from pathlib import Path


COMPONENTS = {
    "pzi": ("pzi_A1", "pzi_B1"),
    "psp": ("psp_A1", "psp_B1"),
    "ppa": ("ppa_A", "ppa_B"),
}
PAIRS = {
    "pzi_psp": (("pzi_A1", "psp_A1"), ("pzi_B1", "psp_B1")),
    "psp_ppa": (("psp_A1", "ppa_A"), ("psp_B1", "ppa_B")),
}
SEQIDS = [f"C{i}A" for i in range(1, 9)] + [f"C{i}B" for i in range(1, 9)]


def require(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f"Missing or empty {label}: {path}")


def read_bed(path: Path) -> list[str]:
    rows = []
    seen = set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 4 or fields[3] in seen:
                continue
            seen.add(fields[3])
            rows.append(raw.rstrip("\r\n"))
    return rows


def block_rows(path: Path, valid: set[str]) -> list[str]:
    rows = []
    seen = set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#") or line.startswith("StartGeneA"):
                continue
            if "*" in line:
                line = line.split("*", 1)[1].strip()
            fields = line.split()
            if len(fields) != 6 or not all(g in valid for g in fields[:4]):
                continue
            key = tuple(fields)
            if key in seen:
                continue
            seen.add(key)
            rows.append("\t".join(fields))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jcvi-root", type=Path, required=True)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--out-root", type=Path, required=True)
    args = parser.parse_args()

    bed_out = args.out_root / "bed"
    anchor_out = args.out_root / "anchors"
    bed_out.mkdir(parents=True, exist_ok=True)
    anchor_out.mkdir(parents=True, exist_ok=True)

    valid: dict[str, set[str]] = {}
    for species, components in COMPONENTS.items():
        rows = []
        seen = set()
        for component in components:
            source = args.input_root / f"{component}.bed"
            require(source, f"JCVI BED {component}")
            for row in read_bed(source):
                gene = row.split("\t", 4)[3]
                if gene not in seen:
                    seen.add(gene)
                    rows.append(row)
        target = bed_out / f"{species}.bed"
        target.write_text("\n".join(rows) + "\n", encoding="utf-8")
        valid[species] = seen

    counts = {}
    for name, pair_list in PAIRS.items():
        target = anchor_out / f"{name}.blocks.simple"
        rows = []
        for left, right in pair_list:
            left_species = left.split("_", 1)[0]
            right_species = right.split("_", 1)[0]
            source = args.jcvi_root / f"{left}_{right}" / "blocks.simple"
            require(source, f"JCVI blocks {left}-{right}")
            rows.extend(block_rows(source, valid[left_species] | valid[right_species]))
        target.write_text("\n".join(rows) + ("\n" if rows else ""), encoding="utf-8")
        counts[name] = len(rows)

    seq_line = ",".join(SEQIDS)
    (args.out_root / "three_species_A1B1.seqids").write_text(
        "\n".join([seq_line] * 3) + "\n", encoding="ascii"
    )
    layout = ["# y, xstart, xend, rotation, color, label, va, bed"]
    for y, species in ((0.90, "pzi"), (0.50, "psp"), (0.10, "ppa")):
        layout.append(f" {y:.2f}, .08, .92, 0, #4C5660, {species.upper()}, center, bed/{species}.bed")
    layout.extend(
        [
            "# A1/B1 direct JCVI comparisons",
            " e, 0, 1, anchors/pzi_psp.blocks.simple",
            " e, 1, 2, anchors/psp_ppa.blocks.simple",
        ]
    )
    (args.out_root / "three_species_A1B1.layout").write_text("\n".join(layout) + "\n", encoding="ascii")
    with (args.out_root / "block_summary.tsv").open("w", encoding="ascii") as handle:
        handle.write("pair\tblocks\tfile\n")
        for name, count in counts.items():
            handle.write(f"{name}\t{count}\tanchors/{name}.blocks.simple\n")
    print(f"[OK] Prepared three-species ribbon inputs: {args.out_root}")
    for name, count in counts.items():
        print(f"[OK] {name}: {count} blocks")


if __name__ == "__main__":
    main()
