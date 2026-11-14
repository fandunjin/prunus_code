#!/usr/bin/env python3
"""Contract low-UFBoot branches and combine best gene trees."""

import argparse
from pathlib import Path
from phylo_newick import contract_low_support, parse_newick, serialize


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree-dir", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--threshold", type=float, default=10.0)
    args = parser.parse_args()
    paths = sorted(Path(args.tree_dir).glob("*.treefile"))
    if not paths:
        raise SystemExit("ERROR: no IQ-TREE gene trees found")
    with Path(args.out).open("w", encoding="ascii", newline="\n") as output:
        for path in paths:
            root = parse_newick(path.read_text(encoding="ascii"))
            root = contract_low_support(root, args.threshold)
            output.write(serialize(root) + ";\n")


if __name__ == "__main__":
    main()
