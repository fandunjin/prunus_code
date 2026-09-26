#!/usr/bin/env python3
"""Create chromosome-renamed PZI haplotype inputs using C01-C08 aliases.

The original FASTA/GFF/BED/JCVI inputs are never modified.  This produces a
parallel, coordinate-compatible package in which chromosome sequence IDs are
named C01_A1 ... C08_B2.  Gene IDs and all feature coordinates are preserved.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

from prepare_pzi_four_haplotype_circos_jcvi import (
    HAPS,
    PAIR_LABELS,
    COPY1_PAIR_ROWS,
    annotation_paths,
    genome_path,
    load_quartet_assignment,
    open_text,
    require_file,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args()


def rewrite_fasta(source: Path, target: Path, aliases: dict[str, str]) -> int:
    count = 0
    with open_text(source) as inp, target.open("w", encoding="utf-8") as out:
        for raw in inp:
            if raw.startswith(">"):
                header = raw.rstrip("\n")
                fields = header[1:].split(None, 1)
                original = fields[0]
                renamed = aliases.get(original, original)
                suffix = f" {fields[1]}" if len(fields) == 2 else ""
                out.write(f">{renamed}{suffix}\n")
                count += 1
            else:
                out.write(raw)
    return count


def rewrite_gff(source: Path, target: Path, aliases: dict[str, str]) -> tuple[int, int]:
    records = 0
    renamed = 0
    with open_text(source) as inp, target.open("w", encoding="utf-8") as out:
        for raw in inp:
            if raw.startswith("#") or not raw.strip():
                out.write(raw)
                continue
            fields = raw.rstrip("\n").split("\t")
            if fields and fields[0] in aliases:
                fields[0] = aliases[fields[0]]
                renamed += 1
            out.write("\t".join(fields) + "\n")
            records += 1
    return records, renamed


def rewrite_bed(source: Path, target: Path, aliases: dict[str, str]) -> int:
    records = 0
    with open_text(source) as inp, target.open("w", encoding="utf-8") as out:
        for raw in inp:
            fields = raw.rstrip("\n").split("\t")
            if fields and fields[0] in aliases:
                fields[0] = aliases[fields[0]]
            out.write("\t".join(fields) + "\n")
            if raw.strip():
                records += 1
    return records


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    out_dir = args.out_dir.resolve()
    genome_dir = out_dir / "genome"
    gff_dir = out_dir / "gff"
    cds_dir = out_dir / "cds"
    pep_dir = out_dir / "pep"
    bed_dir = out_dir / "bed"
    for directory in (genome_dir, gff_dir, cds_dir, pep_dir, bed_dir):
        directory.mkdir(parents=True, exist_ok=True)

    aliases, reverse_aliases, quartet_meta = load_quartet_assignment(root)
    with (out_dir / "chromosome_rename_map.tsv").open("w", encoding="ascii") as handle:
        handle.write(
            "display_id\toriginal_id\thaplotype\tdisplay_group\tsource_group\t"
            "source_group_confidence\tmapping_note\n"
        )
        for display in sorted(reverse_aliases):
            original = reverse_aliases[display]
            hap = display.rsplit("_", 1)[1]
            display_group = display.rsplit("_", 1)[0]
            source_number = int(original[3:5])
            source_group = f"C{source_number:02d}"
            confidence = quartet_meta[source_number][1]
            mapping_note = (
                "manual_A_C01_C06_swap"
                if display_group != source_group
                else "quartet_default"
            )
            handle.write(
                f"{display}\t{original}\t{hap}\t{display_group}\t{source_group}\t"
                f"{confidence}\t{mapping_note}\n"
            )

    with (out_dir / "copy1_pair_map.tsv").open("w", encoding="ascii") as handle:
        handle.write("group\tcopy1_A1_slot\tcopy1_B1_slot\tcopy2_A2_slot\tcopy2_B2_slot\tstatus\n")
        for group, a1_slot, b1_slot, a2_slot, b2_slot, status in COPY1_PAIR_ROWS:
            handle.write(
                f"{group}\t{a1_slot}\t{b1_slot}\t{a2_slot}\t{b2_slot}\t{status}\n"
            )

    manifest_rows: list[tuple[str, str, str, str, str, str]] = []
    for hap in HAPS:
        label = f"pzi_hap{hap}"
        genome = genome_path(root, hap)
        gff, pep = annotation_paths(root, hap)
        genome_out = genome_dir / f"{label}.C_alias.fa"
        gff_out = gff_dir / f"{label}.C_alias.gff3"
        cds_out = cds_dir / f"{label}.C_alias.cds.fa"
        pep_out = pep_dir / f"{label}.C_alias.pep.fa"
        nseq = rewrite_fasta(genome, genome_out, aliases)
        ngff, nrename = rewrite_gff(gff, gff_out, aliases)
        cds_candidates = (
            root / "downstream_inputs" / "subgenomes" / f"{label}.primary.cds.fa",
            root / "final_haplotype_genomes_eviann" / f"{label}.cds.fa",
        )
        cds_source = next((candidate for candidate in cds_candidates if candidate.is_file() and candidate.stat().st_size > 0), None)
        if cds_source is not None:
            shutil.copy2(cds_source, cds_out)
        else:
            cds_out.write_text("", encoding="ascii")
        shutil.copy2(pep, pep_out)

        bed_source = root / "downstream_analysis" / "circos" / "work" / "jcvi_inputs" / f"hap{hap}.bed"
        bed_out = bed_dir / f"hap{hap}.C_alias.bed"
        nbed = 0
        if bed_source.is_file() and bed_source.stat().st_size > 0:
            nbed = rewrite_bed(bed_source, bed_out, aliases)
        manifest_rows.append((hap, str(genome_out), str(gff_out), str(cds_out), str(pep_out), f"{nseq};{ngff};{nrename};{nbed}"))

    with (out_dir / "manifest.tsv").open("w", encoding="utf-8") as handle:
        handle.write("haplotype\tgenome\tgff3\tcds\tprotein\tcounts(seq,gff,renamed_gff,bed)\n")
        for row in manifest_rows:
            handle.write("\t".join(row) + "\n")
    (out_dir / "README.txt").write_text(
        "C01-C08 chromosome aliases derived from the final four-haplotype assignment.\n"
        "Original FASTA/GFF/JCVI files are unchanged; gene IDs and coordinates are preserved.\n"
        "The C01/C06 aliases are manually exchanged in both A1 and A2 following JCVI review.\n"
        "C04, C05, C07 and C08 are stable copy1 self-synteny backbones; C01, C02, C03 and C06\n"
        "retain rearrangement/provisional status.\n",
        encoding="ascii",
    )
    print(f"[DONE] C-aliased inputs: {out_dir}")


if __name__ == "__main__":
    main()
