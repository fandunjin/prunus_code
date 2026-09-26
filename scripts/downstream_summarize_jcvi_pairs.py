#!/usr/bin/env python3
"""Summarize JCVI syntenic gene pairs with BED coordinates."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", required=True)
    parser.add_argument("--query", required=True)
    parser.add_argument("--ref", required=True)
    parser.add_argument("--comparison-type", default="external_reference")
    parser.add_argument("--aligns", required=True, help="Output from jcvi.compara.synteny toaligns")
    parser.add_argument("--anchors", required=True, help="Screened JCVI anchors file")
    parser.add_argument("--qbed", required=True)
    parser.add_argument("--sbed", required=True)
    parser.add_argument("--out-pairs", required=True)
    parser.add_argument("--out-blocks", required=True)
    parser.add_argument("--out-status", required=True)
    return parser.parse_args()


def read_bed(path: Path) -> dict[str, tuple[str, str, str, str]]:
    features: dict[str, tuple[str, str, str, str]] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            strand = fields[5] if len(fields) > 5 else "."
            features[fields[3]] = (fields[0], fields[1], fields[2], strand)
    return features


def read_aligns(path: Path) -> list[tuple[str, str, str]]:
    pairs: list[tuple[str, str, str]] = []
    if not path.is_file() or path.stat().st_size == 0:
        return pairs
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                fields = line.split()
            if len(fields) < 3:
                continue
            pairs.append((fields[0], fields[1], fields[2]))
    return pairs


def write_outputs(args: argparse.Namespace) -> None:
    qbed = read_bed(Path(args.qbed))
    sbed = read_bed(Path(args.sbed))
    aligns = read_aligns(Path(args.aligns))

    out_pairs = Path(args.out_pairs)
    out_blocks = Path(args.out_blocks)
    out_status = Path(args.out_status)
    for path in (out_pairs, out_blocks, out_status):
        path.parent.mkdir(parents=True, exist_ok=True)

    block_rows: dict[str, list[tuple[str, str]]] = defaultdict(list)
    with out_pairs.open("w", encoding="utf-8") as out:
        print(
            "\t".join(
                [
                    "comparison",
                    "comparison_type",
                    "block_id",
                    "query_label",
                    "ref_label",
                    "query_gene",
                    "ref_gene",
                    "query_seqid",
                    "query_start0",
                    "query_end",
                    "query_strand",
                    "ref_seqid",
                    "ref_start0",
                    "ref_end",
                    "ref_strand",
                ]
            ),
            file=out,
        )
        for block_id, qgene, sgene in aligns:
            qcoord = qbed.get(qgene, ("NA", "NA", "NA", "NA"))
            scoord = sbed.get(sgene, ("NA", "NA", "NA", "NA"))
            block_rows[block_id].append((qgene, sgene))
            print(
                "\t".join(
                    [
                        args.comparison,
                        args.comparison_type,
                        block_id,
                        args.query,
                        args.ref,
                        qgene,
                        sgene,
                        qcoord[0],
                        qcoord[1],
                        qcoord[2],
                        qcoord[3],
                        scoord[0],
                        scoord[1],
                        scoord[2],
                        scoord[3],
                    ]
                ),
                file=out,
            )

    with out_blocks.open("w", encoding="utf-8") as out:
        print(
            "\t".join(
                [
                    "comparison",
                    "comparison_type",
                    "block_id",
                    "query_label",
                    "ref_label",
                    "n_pairs",
                    "query_seqids",
                    "query_start0_min",
                    "query_end_max",
                    "ref_seqids",
                    "ref_start0_min",
                    "ref_end_max",
                ]
            ),
            file=out,
        )
        for block_id in sorted(block_rows):
            qgenes = [x[0] for x in block_rows[block_id]]
            sgenes = [x[1] for x in block_rows[block_id]]
            qcoords = [qbed[g] for g in qgenes if g in qbed]
            scoords = [sbed[g] for g in sgenes if g in sbed]
            qseqids = sorted({x[0] for x in qcoords})
            sseqids = sorted({x[0] for x in scoords})
            qstarts = [int(x[1]) for x in qcoords if x[1].isdigit()]
            qends = [int(x[2]) for x in qcoords if x[2].isdigit()]
            sstarts = [int(x[1]) for x in scoords if x[1].isdigit()]
            sends = [int(x[2]) for x in scoords if x[2].isdigit()]
            print(
                "\t".join(
                    [
                        args.comparison,
                        args.comparison_type,
                        block_id,
                        args.query,
                        args.ref,
                        str(len(block_rows[block_id])),
                        ",".join(qseqids) if qseqids else "NA",
                        str(min(qstarts)) if qstarts else "NA",
                        str(max(qends)) if qends else "NA",
                        ",".join(sseqids) if sseqids else "NA",
                        str(min(sstarts)) if sstarts else "NA",
                        str(max(sends)) if sends else "NA",
                    ]
                ),
                file=out,
            )

    status = "ok" if aligns else "no_screened_pairs"
    anchors = Path(args.anchors)
    with out_status.open("w", encoding="utf-8") as out:
        print(
            "\t".join(
                [
                    "comparison",
                    "query_label",
                    "ref_label",
                    "comparison_type",
                    "status",
                    "n_blocks",
                    "n_syntenic_pairs",
                    "screened_anchors",
                    "jcvi_toaligns",
                ]
            ),
            file=out,
        )
        print(
            "\t".join(
                [
                    args.comparison,
                    args.query,
                    args.ref,
                    args.comparison_type,
                    status,
                    str(len(block_rows)),
                    str(len(aligns)),
                    str(anchors),
                    str(Path(args.aligns)),
                ]
            ),
            file=out,
        )


def main() -> int:
    args = parse_args()
    write_outputs(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
