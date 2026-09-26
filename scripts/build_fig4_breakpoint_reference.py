#!/usr/bin/env python3
"""Stage the nine-unit Prunus reference used for Fig. 4 breakpoint discovery.

The reference is deliberately restricted to current files under 0903:
PZI hapB1 C2B and C5B are retained as complete source chromosomes, and
GTR Chr03--Chr09 provide the other seven units.  It is a coordinate anchor
for detecting structural junctions, not a reconstructed DNA ancestor.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


UNITS = (
    ("C2B", "pzi_hapB1", "C2B"),
    ("C5B", "pzi_hapB1", "C5B"),
    ("GTR03", "gtr", "Chr03"),
    ("GTR04", "gtr", "Chr04"),
    ("GTR05", "gtr", "Chr05"),
    ("GTR06", "gtr", "Chr06"),
    ("GTR07", "gtr", "Chr07"),
    ("GTR08", "gtr", "Chr08"),
    ("GTR09", "gtr", "Chr09"),
)

TARGETS = (
    ("pzi_hapA1", "PZI", "A"),
    ("pzi_hapA2", "PZI", "A"),
    ("pzi_hapB1", "PZI", "B"),
    ("pzi_hapB2", "PZI", "B"),
    ("psp_hapA1", "PSP", "A"),
    ("psp_hapA2", "PSP", "A"),
    ("psp_hapB1", "PSP", "B"),
    ("psp_hapB2", "PSP", "B"),
    ("ppa_hapA", "PPA", "A"),
    ("ppa_hapB", "PPA", "B"),
    ("pcam", "PCAM", "diploid"),
    ("ppe", "PPE", "diploid"),
    ("gtr", "GTR", "outgroup"),
    ("rch", "RCH", "outgroup"),
)


def read_gff(path: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 6:
                continue
            try:
                int(float(fields[2]))
                int(float(fields[3]))
                int(float(fields[5]))
            except ValueError:
                continue
            rows.append(fields)
    if not rows:
        raise ValueError(f"No usable WGDI GFF rows: {path}")
    return rows


def read_lens(path: Path) -> dict[str, tuple[int, int]]:
    rows: dict[str, tuple[int, int]] = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                rows[fields[0]] = (int(float(fields[1])), int(float(fields[2])))
            except ValueError:
                continue
    return rows


def read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    ident: str | None = None
    header = ""
    chunks: list[str] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line.startswith(">"):
                if ident is not None:
                    records[ident] = (header, "".join(chunks))
                header = line[1:].strip()
                ident = header.split()[0]
                chunks = []
            elif ident is not None:
                chunks.append(line.strip())
    if ident is not None:
        records[ident] = (header, "".join(chunks))
    return records


def write_wrapped(handle, sequence: str, width: int = 80) -> None:
    for start in range(0, len(sequence), width):
        handle.write(sequence[start : start + width] + "\n")


def source_paths(prepared: Path, label: str) -> tuple[Path, Path, Path]:
    base = prepared / label
    return base / f"{label}.gff", base / f"{label}.lens", base / f"{label}.pep.fa"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True, help="0903 directory")
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    root = args.root.resolve()
    prepared = root / "karyotype" / "prepared"
    output = (args.output or root / "karyotype" / "fig4_breakpoint_candidates").resolve()
    reference = output / "reference"
    reference.mkdir(parents=True, exist_ok=True)

    source_cache: dict[str, tuple[list[list[str]], dict[str, tuple[int, int]], dict[str, tuple[str, str]]]] = {}
    for _unit, label, _chrom in UNITS:
        if label in source_cache:
            continue
        gff, lens, pep = source_paths(prepared, label)
        for path in (gff, lens, pep):
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(path)
        source_cache[label] = (read_gff(gff), read_lens(lens), read_fasta(pep))

    gff_out = reference / "prunus9_reference.gff"
    lens_out = reference / "prunus9_reference.lens"
    pep_out = reference / "prunus9_reference.pep.fa"
    unit_manifest = reference / "prunus9_reference_units.tsv"
    manifest_rows: list[dict[str, object]] = []
    with gff_out.open("w", encoding="utf-8") as gff_handle, lens_out.open("w", encoding="utf-8") as lens_handle, pep_out.open("w", encoding="utf-8") as pep_handle:
        for unit, label, chromosome in UNITS:
            gff_rows, lens, peptides = source_cache[label]
            selected = [row[:] for row in gff_rows if row[0] == chromosome]
            if not selected:
                raise ValueError(f"{label} has no GFF rows for {chromosome}")
            selected.sort(key=lambda row: int(float(row[5])))
            length, lens_genes = lens.get(chromosome, (0, 0))
            if length <= 0:
                length = max(int(float(row[3])) for row in selected)
            selected_ids: list[str] = []
            for new_order, row in enumerate(selected, 1):
                old_id = row[1]
                if old_id not in peptides:
                    raise ValueError(f"Missing peptide {old_id} in {label}")
                new_id = f"{unit}__{old_id}"
                row[0] = unit
                row[1] = new_id
                row[5] = str(new_order)
                gff_handle.write("\t".join(row[:6]) + "\n")
                selected_ids.append(old_id)
                header, sequence = peptides[old_id]
                pep_handle.write(f">{new_id} source={label}:{old_id} {header}\n")
                write_wrapped(pep_handle, sequence)
            lens_handle.write(f"{unit}\t{length}\t{len(selected_ids)}\n")
            manifest_rows.append(
                {
                    "reference_unit": unit,
                    "source_label": label,
                    "source_chromosome": chromosome,
                    "source_length_bp": length,
                    "source_gene_count": len(selected_ids),
                    "source_lens_gene_count": lens_genes,
                    "definition": (
                        "ancestral-karyotype reference unit represented by complete extant PZI hapB1 chromosome"
                        if unit in {"C2B", "C5B"}
                        else "ancestral-karyotype reference unit represented by complete GTR chromosome"
                    ),
                }
            )

    with unit_manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest_rows)

    target_rows: list[dict[str, object]] = []
    for label, species, subgenome in TARGETS:
        gff, lens, pep = source_paths(prepared, label)
        for path in (gff, lens, pep):
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(path)
        target_rows.append(
            {
                "target_label": label,
                "species": species,
                "subgenome": subgenome,
                # Relative paths keep the manifest valid from Windows Python
                # and the WSL WGDI environment.
                "gff": gff.relative_to(root).as_posix(),
                "lens": lens.relative_to(root).as_posix(),
                "pep": pep.relative_to(root).as_posix(),
            }
        )
    target_manifest = output / "target_manifest.tsv"
    with target_manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(target_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(target_rows)

    (output / "README.txt").write_text(
        "Fig. 4 breakpoint-discovery reference\n"
        "All inputs are under 0903. The nine ancestral-karyotype reference units are complete PZI hapB1 C2B/C5B and complete GTR Chr03--Chr09.\n"
        "C2B and C5B are used as complete extant-chromosome proxies for the corresponding ancestral units; the reference is not a reconstructed ancestral DNA sequence.\n"
        "WGDI-derived source transitions require cross-species and outgroup review before being labelled as RCT, EEJ or inversion.\n",
        encoding="utf-8",
    )
    print(f"[OK] reference: {reference}")
    print(f"[OK] targets: {target_manifest}")


if __name__ == "__main__":
    main()
