#!/usr/bin/env python3
"""Write a POInT model file with parameters from an optimized tree output."""

from __future__ import print_function

import argparse
import re
from pathlib import Path


PARAM_RE = re.compile(r"^>\s*(FixRate|SwitchProb|Copy1Bias|Copy2Bias):\s*([-+0-9.eE]+)\s*$")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--template", required=True)
    ap.add_argument("--tree-output", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--fixed-ab", action="store_true",
                    help="Allow the fixed A/B template to omit the unidentifiable SwitchProb parameter")
    args = ap.parse_args()

    values = {}
    for line in Path(args.tree_output).read_text(encoding="ascii", errors="replace").splitlines():
        match = PARAM_RE.match(line.strip())
        if match:
            values[match.group(1)] = match.group(2)
    required = ("FixRate",) if args.fixed_ab else ("FixRate", "SwitchProb")
    if args.fixed_ab and "Copy1Bias" not in values and "Copy2Bias" not in values:
        raise SystemExit("ERROR: fixed A/B optimized tree lacks Copy1Bias or Copy2Bias")
    missing = [name for name in required if name not in values]
    if missing:
        raise SystemExit("ERROR: optimized POInT tree lacks parameters: {}".format(", ".join(missing)))

    template_parameters = set()
    for line in Path(args.template).read_text(encoding="ascii").splitlines():
        fields = line.split()
        if fields and fields[0] in values and len(fields) >= 3:
            template_parameters.add(fields[0])

    out_lines = []
    for line in Path(args.template).read_text(encoding="ascii").splitlines():
        fields = line.split()
        if fields and fields[0] in template_parameters and len(fields) >= 3:
            fields[2] = values[fields[0]]
            line = "\t".join(fields)
        out_lines.append(line)
    Path(args.out).write_text("\n".join(out_lines) + "\n", encoding="ascii")
    print("[OK] Wrote optimized POInT model: {}".format(args.out))


if __name__ == "__main__":
    main()
