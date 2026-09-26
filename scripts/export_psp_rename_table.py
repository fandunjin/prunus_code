#!/usr/bin/env python3
"""Convert the PSP constrained-assignment map to a compact naming table."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


HAP_COLUMNS = ("hapA1", "hapA2", "hapB1", "hapB2")


def natural_key(value: str) -> tuple:
    """Sort chr01 ... chr32 numerically while retaining arbitrary source IDs."""
    pieces = re.split(r"(\d+)", value)
    return tuple(int(piece) if piece.isdigit() else piece.lower() for piece in pieces)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Export a chr32/copy1/copy2/A1/A2/B1/B2 PSP rename table."
    )
    parser.add_argument("--rename-map", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    with args.rename_map.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    required = {"old_seqid", "new_seqid", "hap_set", "phase_suffix"}
    missing = required - set(rows[0] if rows else [])
    if missing:
        raise SystemExit(
            "ERROR: rename map lacks required column(s): " + ", ".join(sorted(missing))
        )
    if len(rows) != 32:
        raise SystemExit(f"ERROR: expected 32 rows in rename map, found {len(rows)}")

    table = []
    for row in rows:
        hap_set = row["hap_set"]
        if hap_set not in HAP_COLUMNS:
            raise SystemExit(f"ERROR: unexpected hap_set for {row['old_seqid']}: {hap_set}")

        phase_suffix = row["phase_suffix"]
        if phase_suffix not in {"1", "2"}:
            raise SystemExit(
                f"ERROR: unexpected phase suffix for {row['old_seqid']}: {phase_suffix}"
            )

        output_row = {
            "chr32": row["old_seqid"],
            "copy1": row["old_seqid"] if phase_suffix == "1" else "",
            "copy2": row["old_seqid"] if phase_suffix == "2" else "",
            **{column: "" for column in HAP_COLUMNS},
        }
        output_row[hap_set] = row["new_seqid"]
        table.append(output_row)

    seen_names = [row[column] for row in table for column in HAP_COLUMNS if row[column]]
    if len(set(seen_names)) != 32:
        raise SystemExit("ERROR: final haplotype names are not unique")

    table.sort(key=lambda row: natural_key(row["chr32"]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=("chr32", "copy1", "copy2", *HAP_COLUMNS),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(table)


if __name__ == "__main__":
    main()
