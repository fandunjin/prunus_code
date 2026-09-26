#!/usr/bin/env python3
"""Finalize provenance markers for an already completed formal FraBias run."""

from pathlib import Path
import argparse


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    args = ap.parse_args()
    required = [
        args.root / "01_catalog" / "gtr_gene_order.tsv",
        args.root / "04_quota_align" / "normalized" / "quota_filtered_anchors.tsv",
        args.root / "05_fractbias" / "fractbias_windows.tsv",
        args.root / "06_panel_c" / "figure2_panel_c_fractbias_reference_exact.svg",
    ]
    missing = [str(p) for p in required if not p.is_file() or p.stat().st_size == 0]
    if missing:
        raise SystemExit("ERROR: formal FraBias outputs missing: {}".format(", ".join(missing)))
    (args.root / "software_and_parameters.txt").write_text(
        "method=LAST (protein-feature mode), tandem filter, DAGChainer, Quota Align Merge, Quota Align, FraBias\n"
        "last=lastdb -p -c; lastal -p BL62 -P 2 -i256M -f BlastTab\n"
        "tandem_Nmax=10\n"
        "dagchainer=D20,g10,A5; effective maximum E-value 1e-5\n"
        "quota_align_merge=Dm20\n"
        "quota_align=logical query:GTR 8:1; file-axis determined from merged QA; Nm40; GLPK\n"
        "fractbias=overlapping 100-gene windows on all nine GTR chromosomes\n"
        "tracks=PZI_A,PZI_B,PSP_A,PSP_B,PPA_A,PPA_B,PCAM,PPE\n"
        "display_pairs=C06,C02,C03\n"
        "boundary=quota-filtered anchor retention, not validated gene loss\n"
        "JCVI_or_Dollo_inputs=forbidden\n", encoding="ascii")
    (args.root / "tool_versions.txt").write_text(
        "LAST=1595 (recorded during run)\n"
        "DAGChainer=run_DAG_chainer.pl (version not emitted by executable)\n"
        "QuotaAlign=local quota-alignment checkout; GLPK 5.0\n"
        "FraBias=local I/O port of CoGe fractionation_bias.py; algorithm unchanged\n", encoding="ascii")
    (args.root / "COMPLETE").write_text("Formal reference-method FraBias panel completed.\n", encoding="ascii")
    print("[OK] Formal FraBias provenance finalized: {}".format(args.root))


if __name__ == "__main__":
    main()
