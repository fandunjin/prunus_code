#!/usr/bin/env python3
"""Build the confirmed-only A1/A2/B1/B2 copy-fate matrix.

This is a lightweight post-processing entry point. It reads the confirmed
candidate, anchor, and synteny-component outputs already produced by
run_inflorescence_family_pipeline.py and writes
copy_fate/confirmed_ab_copy_fate.tsv.
"""
from __future__ import annotations

import argparse
from pathlib import Path

from run_inflorescence_family_pipeline import DEFAULT_OUT, build_confirmed_ab_copy_fate


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    build_confirmed_ab_copy_fate(args.outdir.resolve())
    print(args.outdir.resolve() / "copy_fate" / "confirmed_ab_copy_fate.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
