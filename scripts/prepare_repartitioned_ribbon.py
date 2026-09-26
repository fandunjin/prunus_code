#!/usr/bin/env python3
"""Create native JCVI karyotype inputs for the revised PZI four-haplotype ribbon."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


ROWS = ("pzi_A1", "pzi_A2", "pzi_B1", "pzi_B2")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jcvi-root", type=Path, required=True)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--out-root", type=Path, required=True)
    args = parser.parse_args()
    bed_dir = args.out_root / "bed"
    anchor_dir = args.out_root / "anchors"
    bed_dir.mkdir(parents=True, exist_ok=True)
    anchor_dir.mkdir(parents=True, exist_ok=True)
    for row in ROWS:
        source = args.input_root / f"{row}.bed"
        if not source.is_file() or source.stat().st_size == 0:
            raise SystemExit(f"Missing JCVI BED: {source}")
        target = bed_dir / f"{row}.bed"
        if target.exists() or target.is_symlink():
            target.unlink()
        try:
            target.symlink_to(source.resolve())
        except OSError:
            shutil.copy2(source, target)

    pairs = (("pzi_A1", "pzi_A2"), ("pzi_B1", "pzi_B2"), ("pzi_A1", "pzi_B1"), ("pzi_A2", "pzi_B2"))
    for left, right in pairs:
        source = args.jcvi_root / f"{left}_{right}" / "blocks.simple"
        if not source.is_file() or source.stat().st_size == 0:
            raise SystemExit(f"Missing JCVI blocks: {source}")
        target = anchor_dir / f"{left}_{right}.blocks.simple"
        if target.exists() or target.is_symlink():
            target.unlink()
        try:
            target.symlink_to(source.resolve())
        except OSError:
            shutil.copy2(source, target)

    chroms = [f"C{i}A" for i in range(1, 9)]
    bchroms = [f"C{i}B" for i in range(1, 9)]
    (args.out_root / "repartitioned.seqids").write_text(
        ",".join(chroms) + "\n" + ",".join(chroms) + "\n" + ",".join(bchroms) + "\n" + ",".join(bchroms) + "\n",
        encoding="ascii",
    )
    lines = ["# y, xstart, xend, rotation, color, label, va, bed"]
    labels = ("PZI A1", "PZI A2", "PZI B1", "PZI B2")
    for i, (row, label) in enumerate(zip(ROWS, labels)):
        y = 0.92 - i * 0.27
        lines.append(f" {y:.3f}, .08, .92, 0, #4C5660, {label}, center, bed/{row}.bed")
    lines.extend(
        [
            "# revised partition: gray within-phase and colored A/B comparisons",
            " e, 0, 1, anchors/pzi_A1_pzi_A2.blocks.simple",
            " e, 2, 3, anchors/pzi_B1_pzi_B2.blocks.simple",
            " e, 0, 2, anchors/pzi_A1_pzi_B1.blocks.simple",
            " e, 1, 3, anchors/pzi_A2_pzi_B2.blocks.simple",
        ]
    )
    (args.out_root / "repartitioned.layout").write_text("\n".join(lines) + "\n", encoding="ascii")
    (args.out_root / "README.txt").write_text(
        "Native JCVI ribbon for revised A/B partition. Chromosome labels are C1A-C8A and C1B-C8B.\n"
        "Edges use the same JCVI blocks as the Circos links.\n",
        encoding="ascii",
    )
    print(f"[OK] Native ribbon inputs: {args.out_root}")


if __name__ == "__main__":
    main()
