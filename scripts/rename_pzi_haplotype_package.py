#!/usr/bin/env python3
"""Create a self-contained chromosome-renamed PZI haplotype package.

The source genomes are the final four haplotype FASTA files in
``final_haplotype_genomes``.  Matching primary EviAnn annotations are taken
from ``downstream_inputs/subgenomes``.  Only chromosome/sequence IDs are
changed; gene, transcript, CDS and protein IDs and all coordinates remain
unchanged.

The seven-column table ``pzi_chr32_copy_haplotype_naming.tsv`` is treated as
the authoritative display-name map.  Each final FASTA header contains an
``old_chr=chrNN`` tag, which connects the compact source ID (``chr01_A1``)
to the original 32-chromosome row in that table.  This avoids accidentally
using the older pre-C01/C06-swap mapping.
"""

from __future__ import annotations

import argparse
import hashlib
import re
from pathlib import Path
from typing import Iterable, TextIO


HAPS = ("A1", "A2", "B1", "B2")
TABLE_COLUMNS = ("hapA1", "hapA2", "hapB1", "hapB2")
SOURCE_ID_RE = re.compile(r"^chr\d{2}_(?:A1|A2|B1|B2)$")
OLD_CHR_RE = re.compile(r"(?:^|\s)old_chr=(chr\d{2})(?:\s|$)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument(
        "--map",
        type=Path,
        help="Seven-column naming table (default: final_haplotype_genomes/pzi_chr32_copy_haplotype_naming.tsv)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Allow replacing files already present in --out-dir",
    )
    return parser.parse_args()


def open_text(path: Path) -> TextIO:
    # The current inputs are plain text, but accepting gzip makes the helper
    # safe to reuse with archived copies without changing output semantics.
    import gzip

    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def require_file(path: Path, label: str) -> Path:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty {label}: {path}")
    return path


def parse_naming_table(path: Path) -> tuple[dict[str, dict[str, str]], list[str]]:
    """Return ``old_chr -> column -> display ID`` and preserve row order."""

    rows: dict[str, dict[str, str]] = {}
    order: list[str] = []
    with open_text(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        required = {"original_chr32", "copy1", "copy2", *TABLE_COLUMNS}
        missing = required - set(header)
        if missing:
            raise ValueError(f"Naming table lacks columns {sorted(missing)}: {path}")
        index = {name: i for i, name in enumerate(header)}
        for line_no, raw in enumerate(handle, 2):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < len(header):
                raise ValueError(f"Malformed naming-table row {line_no}: {raw.rstrip()}")
            old = fields[index["original_chr32"]]
            if not re.fullmatch(r"chr\d{2}", old):
                raise ValueError(f"Invalid original_chr32 {old!r} at row {line_no}")
            if old in rows:
                raise ValueError(f"Duplicate original_chr32 {old} at row {line_no}")
            rows[old] = {name: fields[index[name]] for name in header}
            order.append(old)
    if len(rows) != 32:
        raise ValueError(f"Expected 32 chromosome rows in {path}, found {len(rows)}")
    return rows, order


def fasta_records(path: Path) -> Iterable[tuple[str, str, str]]:
    """Yield ``(first_token, complete_header, sequence)`` records."""

    name: str | None = None
    header = ""
    chunks: list[str] = []
    with open_text(path) as handle:
        for raw in handle:
            if raw.startswith(">"):
                if name is not None:
                    yield name, header, "".join(chunks)
                header = raw.rstrip("\r\n")[1:]
                name = header.split(None, 1)[0]
                chunks = []
            else:
                chunks.append(raw.strip())
    if name is not None:
        yield name, header, "".join(chunks)


def read_genome_aliases(
    path: Path,
    hap: str,
    table: dict[str, dict[str, str]],
) -> tuple[dict[str, str], dict[str, str], dict[str, str]]:
    """Build source-ID -> display-ID aliases from final FASTA headers."""

    aliases: dict[str, str] = {}
    display_to_source: dict[str, str] = {}
    source_headers: dict[str, str] = {}
    column = f"hap{hap}"
    for source_id, header, _sequence in fasta_records(path):
        if not SOURCE_ID_RE.fullmatch(source_id):
            raise ValueError(f"Unexpected chromosome ID {source_id!r} in {path}")
        match = OLD_CHR_RE.search(header)
        if match is None:
            raise ValueError(f"FASTA header lacks old_chr tag: >{header}")
        old_chr = match.group(1)
        display = table.get(old_chr, {}).get(column, "")
        if not display:
            raise ValueError(f"No {column} display ID for {old_chr} in naming table")
        if not re.fullmatch(r"C\d{2}_(?:A1|A2|B1|B2)", display):
            raise ValueError(f"Invalid display ID {display!r} for {old_chr}")
        if source_id in aliases and aliases[source_id] != display:
            raise ValueError(f"Conflicting alias for {source_id}")
        if display in display_to_source and display_to_source[display] != source_id:
            raise ValueError(f"Display ID {display} maps to multiple source IDs")
        aliases[source_id] = display
        display_to_source[display] = source_id
        source_headers[source_id] = header
    if len(aliases) != 8:
        raise ValueError(f"Expected 8 chromosomes for hap{hap}, found {len(aliases)}")
    return aliases, display_to_source, source_headers


def replace_ids(text: str, aliases: dict[str, str]) -> str:
    # IDs are delimited by non-word characters in FASTA descriptions and GFF
    # attributes.  Longest-first replacement prevents partial collisions.
    for source, display in sorted(aliases.items(), key=lambda item: -len(item[0])):
        text = re.sub(rf"(?<![A-Za-z0-9_]){re.escape(source)}(?![A-Za-z0-9_])", display, text)
    return text


def rewrite_fasta(path: Path, target: Path, aliases: dict[str, str]) -> tuple[int, int, int]:
    records = changed_headers = bases = 0
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="\n") as out:
        for source_id, header, sequence in fasta_records(path):
            new_header = replace_ids(header, aliases)
            # The first token must be the new chromosome ID for genome FASTA;
            # for CDS/PEP, no chromosome token is normally present.
            if source_id in aliases:
                new_header = aliases[source_id] + (" " + header.split(None, 1)[1] if " " in header else "")
            if new_header != header:
                changed_headers += 1
            out.write(f">{new_header}\n")
            for start in range(0, len(sequence), 80):
                out.write(sequence[start : start + 80] + "\n")
            records += 1
            bases += len(sequence)
    return records, changed_headers, bases


def rewrite_gff(path: Path, target: Path, aliases: dict[str, str]) -> tuple[int, int, int, int]:
    records = renamed_seqids = renamed_regions = transcripts = 0
    target.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path) as inp, target.open("w", encoding="utf-8", newline="\n") as out:
        for raw in inp:
            line = raw.rstrip("\r\n")
            if line.startswith("##sequence-region"):
                fields = line.split()
                if len(fields) >= 2 and fields[1] in aliases:
                    fields[1] = aliases[fields[1]]
                    renamed_regions += 1
                out.write(" ".join(fields) + "\n")
                continue
            if line.startswith("#") or not line.strip():
                out.write(line + "\n")
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                raise ValueError(f"Malformed GFF row in {path}: {line}")
            if fields[0] in aliases:
                fields[0] = aliases[fields[0]]
                renamed_seqids += 1
            if fields[2] in {"mRNA", "transcript"}:
                transcripts += 1
            records += 1
            out.write("\t".join(fields) + "\n")
    return records, renamed_seqids, renamed_regions, transcripts


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def refuse_existing_outputs(out: Path, force: bool) -> None:
    """Protect an existing package from accidental partial overwrites."""

    if force or not out.exists():
        return
    existing = [
        path for path in out.rglob("*")
        if path.is_file() and path.stat().st_size > 0
    ]
    if existing:
        preview = "\n  ".join(str(path) for path in existing[:10])
        suffix = "\n  ..." if len(existing) > 10 else ""
        raise FileExistsError(
            f"Output directory already contains files; refusing to overwrite: {out}\n  "
            f"{preview}{suffix}\nUse --force only after checking the existing package."
        )


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    out = args.out_dir.resolve()
    refuse_existing_outputs(out, args.force)
    if args.map:
        map_path = args.map.resolve()
    else:
        map_candidates = (
            root / "final_haplotype_genomes" / "pzi_chr32_copy_haplotype_naming.tsv",
            root / "haploid" / "final_haplotypes" / "pzi_chr32_copy_haplotype_naming.tsv",
        )
        map_path = next((path.resolve() for path in map_candidates if path.is_file() and path.stat().st_size > 0), map_candidates[0].resolve())
    require_file(map_path, "seven-column chromosome naming table")
    table, table_order = parse_naming_table(map_path)

    dirs = {name: out / name for name in ("genome", "gff", "cds", "pep")}
    for directory in dirs.values():
        directory.mkdir(parents=True, exist_ok=True)

    mapping_rows: list[tuple[str, str, str, str, str, str, str]] = []
    manifest_rows: list[dict[str, object]] = []
    all_outputs: list[Path] = []
    for hap in HAPS:
        label = f"pzi_hap{hap}"
        genome_candidates = (
            root / "final_haplotype_genomes" / f"{label}.fa",
            root / "haploid" / "final_haplotypes" / f"{label}.fa",
        )
        genome_src = next((path for path in genome_candidates if path.is_file() and path.stat().st_size > 0), genome_candidates[0])
        require_file(genome_src, f"{label} genome")
        gff_src = require_file(root / "downstream_inputs" / "subgenomes" / f"{label}.primary.gff3", f"{label} primary GFF3")
        cds_src = require_file(root / "downstream_inputs" / "subgenomes" / f"{label}.primary.cds.fa", f"{label} primary CDS")
        pep_src = require_file(root / "downstream_inputs" / "subgenomes" / f"{label}.primary.pep.fa", f"{label} primary protein")
        aliases, _display_to_source, source_headers = read_genome_aliases(genome_src, hap, table)

        genome_out = dirs["genome"] / f"{label}.fa"
        gff_out = dirs["gff"] / f"{label}.gff3"
        cds_out = dirs["cds"] / f"{label}.cds.fa"
        pep_out = dirs["pep"] / f"{label}.pep.fa"
        genome_records, genome_changed, genome_bp = rewrite_fasta(genome_src, genome_out, aliases)
        gff_records, gff_seqids, gff_regions, gff_transcripts = rewrite_gff(gff_src, gff_out, aliases)
        cds_records, cds_changed, cds_bp = rewrite_fasta(cds_src, cds_out, aliases)
        pep_records, pep_changed, pep_aa = rewrite_fasta(pep_src, pep_out, aliases)
        all_outputs.extend((genome_out, gff_out, cds_out, pep_out))

        # Record each mapping once, including the copy columns for downstream
        # interpretation and a direct link to the original 32-chromosome ID.
        for source_id, display in sorted(aliases.items()):
            source_header = source_headers[source_id]
            old_match = OLD_CHR_RE.search(source_header)
            assert old_match is not None
            old_chr = old_match.group(1)
            row = table[old_chr]
            mapping_rows.append((old_chr, source_id, display, hap, row["copy1"], row["copy2"], source_header))

        # GFF seqids must all be represented in the genome mapping.  This
        # catches mismatched annotation packages before they reach Circos/JCVI.
        gff_seqid_set: set[str] = set()
        with open_text(gff_src) as handle:
            for raw in handle:
                if raw.startswith("#") or not raw.strip():
                    continue
                gff_seqid_set.add(raw.split("\t", 1)[0])
        unknown_gff = sorted(gff_seqid_set - set(aliases))
        if unknown_gff:
            raise ValueError(f"{label} GFF contains unmapped seqids: {unknown_gff[:5]}")
        manifest_rows.append({
            "haplotype": hap,
            "genome_source": genome_src.relative_to(root).as_posix(),
            "gff3_source": gff_src.relative_to(root).as_posix(),
            "cds_source": cds_src.relative_to(root).as_posix(),
            "protein_source": pep_src.relative_to(root).as_posix(),
            "genome": genome_out.relative_to(out).as_posix(),
            "gff3": gff_out.relative_to(out).as_posix(),
            "cds": cds_out.relative_to(out).as_posix(),
            "protein": pep_out.relative_to(out).as_posix(),
            "genome_records": genome_records,
            "genome_bp": genome_bp,
            "gff_features": gff_records,
            "gff_transcripts": gff_transcripts,
            "cds_records": cds_records,
            "protein_records": pep_records,
            "genome_headers_changed": genome_changed,
            "cds_headers_changed": cds_changed,
            "protein_headers_changed": pep_changed,
            "gff_seqids_renamed": gff_seqids,
            "gff_regions_renamed": gff_regions,
        })

    mapping_rows.sort(key=lambda row: (int(row[0][3:]), row[3]))
    map_out = out / "chromosome_rename_map.tsv"
    with map_out.open("w", encoding="ascii", newline="\n") as handle:
        handle.write("original_chr32\tsource_seqid\tdisplay_seqid\thaplotype\tcopy1\tcopy2\tsource_header\n")
        for row in mapping_rows:
            handle.write("\t".join(row) + "\n")
    all_outputs.append(map_out)

    manifest_out = out / "manifest.tsv"
    fields = list(manifest_rows[0])
    with manifest_out.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(fields) + "\n")
        for row in manifest_rows:
            handle.write("\t".join(str(row[field]) for field in fields) + "\n")
    all_outputs.append(manifest_out)

    validation_out = out / "validation.tsv"
    with validation_out.open("w", encoding="ascii", newline="\n") as handle:
        handle.write("haplotype\tstatus\tgenome_records\tgff_features\tgff_transcripts\tcds_records\tprotein_records\tgff_seqids_renamed\n")
        for row in manifest_rows:
            status = "PASS" if (
                row["genome_records"] == 8
                and row["gff_transcripts"] == row["cds_records"] == row["protein_records"]
                and row["gff_seqids_renamed"] > 0
            ) else "CHECK"
            handle.write("\t".join(map(str, (row["haplotype"], status, row["genome_records"], row["gff_features"], row["gff_transcripts"], row["cds_records"], row["protein_records"], row["gff_seqids_renamed"]))) + "\n")
    all_outputs.append(validation_out)

    checksums_out = out / "checksums.md5"
    with checksums_out.open("w", encoding="ascii", newline="\n") as handle:
        for path in sorted(all_outputs):
            handle.write(f"{md5sum(path)}  {path.relative_to(out).as_posix()}\n")

    readme = out / "README.txt"
    readme.write_text(
        "P. zippeliana chromosome-renamed four-haplotype package\n"
        "========================================================\n"
        "Genome FASTA files are from final_haplotype_genomes. Primary GFF3/CDS/protein\n"
        "files are matched files from downstream_inputs/subgenomes. Chromosome IDs use\n"
        "the C01-C08 aliases in pzi_chr32_copy_haplotype_naming.tsv; gene/transcript/\n"
        "protein IDs and coordinates are unchanged. GFF seqids and sequence-region\n"
        "directives are rewritten to the display IDs. See manifest.tsv, validation.tsv\n"
        "and checksums.md5 for provenance and integrity checks.\n",
        encoding="ascii",
    )
    print(f"[DONE] Renamed haplotype package: {out}")
    for row in manifest_rows:
        print(f"[PASS] hap{row['haplotype']}: genome={row['genome_records']} sequences, GFF transcripts={row['gff_transcripts']}, CDS={row['cds_records']}, PEP={row['protein_records']}")


if __name__ == "__main__":
    main()
