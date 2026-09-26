#!/usr/bin/env python3
"""Create seqids/layout files for native JCVI karyotype plotting."""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query", required=True)
    parser.add_argument("--ref", required=True)
    parser.add_argument("--query-seqids", required=True)
    parser.add_argument("--ref-seqids", required=True)
    parser.add_argument("--query-bed", required=True)
    parser.add_argument("--ref-bed", required=True)
    parser.add_argument("--anchors-simple", required=True)
    parser.add_argument("--out-seqids", required=True)
    parser.add_argument("--out-layout", required=True)
    return parser.parse_args()


def read_seqids(path: Path) -> str:
    if path.is_file():
        text = path.read_text(encoding="utf-8", errors="replace").strip()
        if text:
            return text
    raise SystemExit(f"Missing or empty seqids file: {path}")


def main() -> int:
    args = parse_args()
    query_seqids = read_seqids(Path(args.query_seqids))
    ref_seqids = read_seqids(Path(args.ref_seqids))
    out_seqids = Path(args.out_seqids)
    out_layout = Path(args.out_layout)
    out_seqids.parent.mkdir(parents=True, exist_ok=True)
    out_layout.parent.mkdir(parents=True, exist_ok=True)

    out_seqids.write_text(query_seqids + "\n" + ref_seqids + "\n", encoding="utf-8")
    layout = "\n".join(
        [
            "# y, xstart, xend, rotation, color, label, va, bed",
            f" .65, .08, .92, 0, #4477AA, {args.query}, top, {Path(args.query_bed).name}",
            f" .35, .08, .92, 0, #CC6677, {args.ref}, bottom, {Path(args.ref_bed).name}",
            "# edges",
            f" e, 0, 1, {Path(args.anchors_simple).name}",
            "",
        ]
    )
    out_layout.write_text(layout, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
