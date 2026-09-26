#!/usr/bin/env python3
"""Increase native WGDI point opacity without changing point coordinates or count."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


POINT_STYLE = re.compile(r'(<use\b[^>]*?style="[^"]*?fill-opacity:\s*)0\.5(\s*"[^>]*>)')


def darken_file(source: Path, destination: Path, opacity: float) -> tuple[int, int]:
    text = source.read_text(encoding="utf-8")
    before = len(re.findall(r"<use\b", text))
    text, changed = POINT_STYLE.subn(rf"\g<1>{opacity:.2f}\g<2>", text)
    after = len(re.findall(r"<use\b", text))
    if before != after:
        raise ValueError(f"point count changed for {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(text, encoding="utf-8")
    return changed, after


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--opacity", type=float, default=0.85)
    args = parser.parse_args()
    if not 0 < args.opacity <= 1:
        raise SystemExit("--opacity must be between 0 and 1")

    rows = []
    for source in sorted(args.input_dir.glob("*.svg")):
        destination = args.output_dir / source.name
        changed, points = darken_file(source, destination, args.opacity)
        rows.append((source.name, changed, points, destination.stat().st_size))
        print(f"{source.name}\tpoints={points}\tchanged={changed}\topacity={args.opacity:.2f}")
    if not rows:
        raise SystemExit(f"No SVG files found in {args.input_dir}")
    print(f"Processed {len(rows)} SVG files")


if __name__ == "__main__":
    main()
