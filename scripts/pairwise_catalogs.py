#!/usr/bin/env python3
"""Pool final A and B locus catalogs into the three Fig. 2 pairwise inputs."""
import argparse
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--catalog", type=Path, required=True, help="fractionation/01_catalog/tracks"
    )
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    outputs = {}
    for species in ["pzi", "psp", "ppa"]:
        for suffix in ["bed", "pep"]:
            target = args.out / f"{species}.{suffix}"
            if target.exists():
                raise SystemExit(f"Output exists: {target}")
            outputs[target] = b"".join(
                (args.catalog / f"{species.upper()}_{subgenome}.{suffix}").read_bytes()
                for subgenome in ["A", "B"]
            )
    args.out.mkdir(parents=True, exist_ok=True)
    checksums = {}
    for target, data in outputs.items():
        target.write_bytes(data)
        checksums[target.name] = hashlib.sha256(data).hexdigest()
    (args.out / "pair_input_sha256.json").write_text(
        json.dumps(checksums, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
