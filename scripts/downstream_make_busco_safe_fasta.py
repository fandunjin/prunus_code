#!/usr/bin/env python3
"""Write a BUSCO-safe protein FASTA and a header mapping table.

BUSCO 5 can fail when FASTA headers contain path-like characters such as "/".
This helper strips descriptions and replaces unsafe characters in the record ID
for BUSCO QC only. Canonical downstream FASTA files are not modified.
"""

import argparse
import csv
import re
from pathlib import Path


SAFE_ID_RE = re.compile(r"[^A-Za-z0-9_.-]+")


def fasta_iter(path):
    current = None
    desc = ""
    chunks = []
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if current is not None:
                    yield current, desc, "".join(chunks)
                header = line[1:].strip()
                current = header.split()[0] if header else ""
                desc = header[len(current):].strip() if current else ""
                chunks = []
            else:
                chunks.append(line.strip())
        if current is not None:
            yield current, desc, "".join(chunks)


def safe_id(label, original_id, index):
    clean = SAFE_ID_RE.sub("_", original_id).strip("_.-")
    if not clean:
        clean = "record_%06d" % index
    return "%s__%s" % (label, clean)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--map", required=True)
    parser.add_argument("--label", required=True)
    args = parser.parse_args()

    out_path = Path(args.output)
    map_path = Path(args.map)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    map_path.parent.mkdir(parents=True, exist_ok=True)

    seen = {}
    records = 0
    with out_path.open("w", encoding="utf-8") as out, map_path.open("w", encoding="utf-8", newline="") as map_handle:
        writer = csv.DictWriter(
            map_handle,
            delimiter="\t",
            fieldnames=["label", "safe_id", "original_id", "original_description"],
            lineterminator="\n",
        )
        writer.writeheader()
        for records, (original_id, desc, seq) in enumerate(fasta_iter(args.input), start=1):
            base_id = safe_id(args.label, original_id, records)
            used = seen.get(base_id, 0)
            seen[base_id] = used + 1
            final_id = base_id if used == 0 else "%s_dup%d" % (base_id, used + 1)

            out.write(">%s\n" % final_id)
            for idx in range(0, len(seq), 60):
                out.write(seq[idx:idx + 60] + "\n")
            writer.writerow({
                "label": args.label,
                "safe_id": final_id,
                "original_id": original_id,
                "original_description": desc,
            })

    print("busco_safe_records=%d output=%s map=%s" % (records, out_path, map_path))


if __name__ == "__main__":
    main()
