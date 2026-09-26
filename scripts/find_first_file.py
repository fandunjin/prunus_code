#!/usr/bin/env python3
import argparse
import fnmatch
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Print the first matching file in a directory tree.")
    parser.add_argument("--dir", required=True)
    parser.add_argument("--pattern", default="*")
    return parser.parse_args()


def main():
    args = parse_args()
    root = Path(args.dir)
    if not root.is_dir():
        raise SystemExit(f"ERROR: directory does not exist: {root}")
    for path in sorted(root.rglob("*")):
        if path.is_file() and fnmatch.fnmatch(path.name, args.pattern):
            print(path)
            return
    raise SystemExit(f"ERROR: no file matching {args.pattern!r} under {root}")


if __name__ == "__main__":
    main()
