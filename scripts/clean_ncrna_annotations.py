#!/usr/bin/env python3
"""Build a non-redundant, standard GFF3 ncRNA catalogue from tRNAscan-SE and cmscan tblout.

The historical cmscan ``*.gff3`` files in this project are human-readable
reports, not nine-column GFF3.  This script therefore parses the reproducible
cmscan tblout directly and leaves all raw inputs untouched.
"""

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Tuple


def attrs_to_dict(text: str) -> Dict[str, str]:
    result: Dict[str, str] = {}
    for item in text.rstrip(";").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result


def gff_attrs(items: dict[str, object]) -> str:
    def esc(value: object) -> str:
        return str(value).replace("%", "%25").replace(";", "%3B").replace("=", "%3D").replace("&", "%26").replace(",", "%2C")

    return ";".join(f"{key}={esc(value)}" for key, value in items.items())


def overlap(a: Tuple[int, int], b: Tuple[int, int]) -> float:
    left = max(a[0], b[0])
    right = min(a[1], b[1])
    if right < left:
        return 0.0
    inter = right - left + 1
    return inter / min(a[1] - a[0] + 1, b[1] - b[0] + 1)


def classify(model: str, description: str) -> str:
    value = f"{model} {description}".lower()
    if "rrna" in value or "ribosomal rna" in value:
        return "rRNA"
    if re.search(r"(^|[^a-z])mir[-_]?[a-z0-9]", value) or model.lower().startswith("mir"):
        return "miRNA"
    if model.lower().startswith(("sno", "snor", "snord")) or "small nucleolar" in value:
        return "snoRNA"
    if re.search(r"(^|[^a-z])u(?:1|2|3|4|5|6|7|8|11|12)(?:atac)?([^0-9]|$)", model.lower()):
        return "snRNA"
    if "signal recognition particle" in value or "_srp" in model.lower() or model.lower().endswith("srp"):
        return "SRP RNA"
    if model.lower() == "trna":
        return "tRNA"
    return "other ncRNA"


def parse_trna(path: Path):
    rows = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9 or fields[2] not in {"tRNA", "pseudogene"}:
            continue
        start, end = sorted((int(fields[3]), int(fields[4])))
        attrs = attrs_to_dict(fields[8])
        kind = "tRNA" if fields[2] == "tRNA" else "tRNA pseudogene"
        rows.append(
            {
                "seqid": fields[0],
                "start": start,
                "end": end,
                "strand": fields[6] if fields[6] in {"+", "-"} else ".",
                "score": float(fields[5]) if fields[5] not in {".", ""} else 0.0,
                "feature": fields[2],
                "class": kind,
                "source": "tRNAscan-SE",
                "attrs": attrs,
                "raw": fields,
            }
        )
    return rows


def parse_tblout(path: Path):
    rows = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 29:
            continue
        try:
            start, end = sorted((int(fields[9]), int(fields[10])))
            score = float(fields[16])
            evalue = fields[17]
        except (ValueError, IndexError):
            continue
        if fields[18] != "!":
            continue
        model, accession, seqid, strand = fields[1], fields[2], fields[3], fields[11]
        description = " ".join(fields[28:])
        rows.append(
            {
                "seqid": seqid,
                "start": start,
                "end": end,
                "strand": strand if strand in {"+", "-"} else ".",
                "score": score,
                "evalue": evalue,
                "model": model,
                "accession": accession,
                "class": classify(model, description),
                "description": description,
                "olp": fields[19],
            }
        )
    return rows


def select_nonredundant(rows):
    # Prefer the model most appropriate for a plant genome, then use the
    # Infernal score to break ties.  Rfam contains multiple kingdom-specific
    # models for conserved rRNAs, so score alone can retain redundant calls.
    def rank(row):
        model = row["model"].lower()
        value = 0
        if model in {"5s_rrna", "5_8s_rrna"}:
            value = 90
        elif model in {"ssu_rrna_eukarya", "lsu_rrna_eukarya"}:
            value = 100
        elif "rrna" in model and any(x in model for x in ("bacteria", "archaea", "microsporidia")):
            value = 10
        elif row["class"] in {"miRNA", "snRNA", "snoRNA", "SRP RNA"}:
            value = 80
        else:
            value = 50
        return value

    kept = []
    decisions = []
    for row in sorted(rows, key=lambda x: (-rank(x), -x["score"], x["seqid"], x["start"], x["end"], x["model"])):
        duplicate = None
        for old in kept:
            if old["seqid"] != row["seqid"]:
                continue
            # Different-strand calls are retained unless both are rRNA models;
            # antisense ncRNAs can be genuine, whereas kingdom-specific rRNA
            # models are the known source of same-locus redundancy.
            if old["strand"] != row["strand"] and old["class"] != "rRNA" and row["class"] != "rRNA":
                continue
            if overlap((old["start"], old["end"]), (row["start"], row["end"])) >= 0.5:
                duplicate = old
                break
        if duplicate is None:
            kept.append(row)
            decisions.append((row, "kept", ""))
        else:
            decisions.append((row, "dropped_overlap", f"{duplicate['model']}:{duplicate['start']}-{duplicate['end']}"))
    return kept, decisions


def fasta_seqids(path: Path):
    """Return FASTA record IDs (the first token after '>')."""
    ids = set()
    if path is None or not path.is_file():
        return ids
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                fields = line[1:].strip().split(None, 1)
                if fields:
                    ids.add(fields[0])
    return ids


def locate_input(src: Path, haplotype: str, suffix: str) -> Path:
    """Find a non-empty haplotype input, accepting prefixed filenames."""
    candidates = [src / f"{haplotype}{suffix}", src / f"{src.parent.name}_{haplotype}{suffix}"]
    for candidate in candidates:
        if candidate.is_file() and candidate.stat().st_size > 0:
            return candidate
    matches = sorted(p for p in src.glob(f"*{suffix}") if p.is_file() and p.stat().st_size > 0)
    if len(matches) == 1:
        return matches[0]
    raise SystemExit(f"could not uniquely locate {suffix} in {src}: {matches}")


def write_outputs(species: str, root: Path, out_root: Path, haplotype: str = "", genome: Path = None) -> dict:
    src = root / species if not haplotype else root / species / haplotype
    label = species if not haplotype else f"{species}_{haplotype}"
    out = out_root / species / ("final" if not haplotype else Path(haplotype) / "final")
    out.mkdir(parents=True, exist_ok=True)
    trna_path = src / "tRNAscan-SE.gff3" if not haplotype else locate_input(src, haplotype, ".tRNAscan-SE.gff3")
    tbl_path = src / "rfam.cmscan.tblout" if not haplotype else locate_input(src, haplotype, ".rfam.cmscan.tblout")
    if not trna_path.is_file() or trna_path.stat().st_size == 0:
        raise SystemExit(f"missing tRNAscan GFF3: {trna_path}")
    if not tbl_path.is_file() or tbl_path.stat().st_size == 0:
        raise SystemExit(f"missing cmscan tblout: {tbl_path}")

    trnas = parse_trna(trna_path)
    rfam = parse_tblout(tbl_path)
    # Rfam tRNA models are redundant with tRNAscan-SE whenever their loci
    # overlap; retain the dedicated tRNAscan prediction as the authoritative call.
    rfam_nontrna = []
    decisions = []
    for row in rfam:
        if row["class"] == "tRNA":
            hit = next((t for t in trnas if t["seqid"] == row["seqid"] and overlap((t["start"], t["end"]), (row["start"], row["end"])) >= 0.5), None)
            if hit:
                decisions.append((row, "dropped_tRNA_overlap", f"{hit['seqid']}:{hit['start']}-{hit['end']}"))
            else:
                decisions.append((row, "dropped_tRNA_model", "tRNAscan-SE is authoritative for tRNA loci"))
            continue
        rfam_nontrna.append(row)
    rfam_kept, rfam_decisions = select_nonredundant(rfam_nontrna)
    decisions.extend(rfam_decisions)

    gff_path = out / f"{label}.ncRNA.final.gff3"
    with gff_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("##gff-version 3\n")
        handle.write(f"##source ncRNA-cleaner;species={label};inputs=tRNAscan-SE+Infernal-Rfam-tblout\n")
        for idx, row in enumerate(sorted(trnas, key=lambda x: (x["seqid"], x["start"], x["end"])), 1):
            attrs = dict(row["attrs"])
            attrs.setdefault("ID", f"{label}.trna{idx}")
            attrs["ncRNA_class"] = row["class"]
            handle.write("\t".join([row["seqid"], row["source"], row["feature"], str(row["start"]), str(row["end"]), f"{row['score']:.1f}", row["strand"], ".", gff_attrs(attrs)]) + "\n")
        for idx, row in enumerate(sorted(rfam_kept, key=lambda x: (x["seqid"], x["start"], x["end"])), 1):
            attrs = {
                "ID": f"{label}.rfam{idx}",
                "Name": row["model"],
                "rfam_accession": row["accession"],
                "ncRNA_class": row["class"],
                "evalue": row["evalue"],
                "description": row["description"],
            }
            handle.write("\t".join([row["seqid"], "Infernal-Rfam", "ncRNA", str(row["start"]), str(row["end"]), f"{row['score']:.1f}", row["strand"], ".", gff_attrs(attrs)]) + "\n")

    stats = Counter([t["class"] for t in trnas] + [r["class"] for r in rfam_kept])
    lengths = defaultdict(int)
    for row in trnas + rfam_kept:
        lengths[row["class"]] += row["end"] - row["start"] + 1
    with (out / f"{label}.ncRNA.final.stats.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["species", "class", "count", "total_bp"])
        for cls in ["tRNA", "tRNA pseudogene", "rRNA", "snRNA", "snoRNA", "miRNA", "SRP RNA", "other ncRNA"]:
            writer.writerow([label, cls, stats[cls], lengths[cls]])

    by_seqid = defaultdict(Counter)
    by_seqid_bp = defaultdict(Counter)
    for row in trnas + rfam_kept:
        by_seqid[row["seqid"]][row["class"]] += 1
        by_seqid_bp[row["seqid"]][row["class"]] += row["end"] - row["start"] + 1
    classes = ["tRNA", "tRNA pseudogene", "rRNA", "snRNA", "snoRNA", "miRNA", "SRP RNA", "other ncRNA"]
    with (out / f"{label}.ncRNA.final.by_seqid.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["species", "seqid", "class", "count", "total_bp"])
        for seqid in sorted(by_seqid):
            for cls in classes:
                writer.writerow([label, seqid, cls, by_seqid[seqid][cls], by_seqid_bp[seqid][cls]])

    with (out / f"{label}.rfam.audit.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["species", "model", "accession", "seqid", "start", "end", "strand", "class", "score", "evalue", "decision", "reason"])
        for row, decision, reason in decisions:
            writer.writerow([label, row.get("model", ""), row.get("accession", ""), row["seqid"], row["start"], row["end"], row["strand"], row["class"], row["score"], row.get("evalue", ""), decision, reason])

    seqid_check = "not_requested"
    if genome is not None:
        expected = fasta_seqids(genome)
        observed = {row["seqid"] for row in trnas + rfam_kept}
        unknown = sorted(observed - expected)
        seqid_check = "pass" if not unknown else f"fail:{len(unknown)}_unknown_seqids"
        if unknown:
            raise SystemExit(f"{label}: GFF3 contains sequence IDs absent from {genome}: {unknown[:5]}")

    with (out / f"{label}.ncRNA.final.manifest.txt").open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"Species: {species}\n")
        if haplotype:
            handle.write(f"Haplotype: {haplotype}\n")
        handle.write(f"Label: {label}\n")
        handle.write(f"tRNAscan input: {trna_path}\n")
        handle.write(f"Infernal tblout input: {tbl_path}\n")
        handle.write(f"Genome input: {genome if genome else 'not supplied'}\n")
        handle.write(f"Sequence-ID validation: {seqid_check}\n")
        handle.write("Rfam filter: included hits only (inc=!); all Rfam tRNA models excluded because tRNAscan-SE is authoritative; reciprocal overlap >=0.5 retained at highest plant-model priority and score\n")
        handle.write(f"Final GFF3: {gff_path}\n")
        handle.write(f"tRNA loci retained: {len(trnas)}\n")
        handle.write(f"Rfam hits retained: {len(rfam_kept)}\n")
        handle.write("Note: Rfam model names were mapped to broad classes using model names/descriptions; accession-level curation can refine this mapping.\n")
    print(f"{label}: tRNA_loci={len(trnas)} rfam_retained={len(rfam_kept)} -> {gff_path}")
    return {
        "label": label,
        "species": species,
        "haplotype": haplotype,
        "gff3": str(gff_path),
        "stats": stats,
        "total_bp": lengths,
        "trna_input": str(trna_path),
        "rfam_input": str(tbl_path),
        "genome": "" if genome is None else str(genome),
        "seqid_check": seqid_check,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, default=Path("0903/annotation/ncrna"))
    parser.add_argument("--output-root", type=Path, default=Path("0903/annotation/ncrna_clean"))
    parser.add_argument("--species", nargs="+", default=["pzi", "psp", "ppa"])
    parser.add_argument("--haplotypes", action="store_true", help="process species/haplotype directories")
    parser.add_argument("--genome-root", type=Path, default=None, help="optional root containing species/genome/<species>_<haplotype>.fa")
    args = parser.parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)
    results = []
    if args.haplotypes:
        for species in args.species:
            species_dir = args.input_root / species
            if not species_dir.is_dir():
                raise SystemExit(f"missing species directory: {species_dir}")
            hap_dirs = sorted(p for p in species_dir.iterdir() if p.is_dir())
            if not hap_dirs:
                raise SystemExit(f"no haplotype directories under: {species_dir}")
            for hap_dir in hap_dirs:
                genome = None
                if args.genome_root:
                    genome = args.genome_root / species / "genome" / f"{species}_{hap_dir.name}.fa"
                    if not genome.is_file():
                        raise SystemExit(f"missing haplotype genome for {species}/{hap_dir.name}: {genome}")
                results.append(write_outputs(species, args.input_root, args.output_root, hap_dir.name, genome))
    else:
        for species in args.species:
            results.append(write_outputs(species, args.input_root, args.output_root))

    summary_path = args.output_root / ("ncrna.haploid.summary.tsv" if args.haplotypes else "ncrna.final.summary.tsv")
    with summary_path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        classes = ["tRNA", "tRNA pseudogene", "rRNA", "snRNA", "snoRNA", "miRNA", "SRP RNA", "other ncRNA"]
        if args.haplotypes:
            writer.writerow(["label", "species", "haplotype", "class", "count", "total_bp", "seqid_check", "gff3"])
            for result in results:
                for cls in classes:
                    writer.writerow([result["label"], result["species"], result["haplotype"], cls, result["stats"][cls], result["total_bp"][cls], result["seqid_check"], result["gff3"]])
        else:
            writer.writerow(["species", "class", "count", "total_bp"])
            for result in results:
                for cls in classes:
                    writer.writerow([result["species"], cls, result["stats"][cls], result["total_bp"][cls]])

    if args.haplotypes:
        species_summary = args.output_root / "ncrna.haploid.species_summary.tsv"
        with species_summary.open("w", encoding="utf-8", newline="\n") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["species", "class", "haplotypes", "count_sum", "total_bp_sum", "count_min", "count_max"])
            for species in args.species:
                species_results = [r for r in results if r["species"] == species]
                for cls in ["tRNA", "tRNA pseudogene", "rRNA", "snRNA", "snoRNA", "miRNA", "SRP RNA", "other ncRNA"]:
                    counts = [r["stats"][cls] for r in species_results]
                    bps = [r["total_bp"][cls] for r in species_results]
                    writer.writerow([species, cls, len(species_results), sum(counts), sum(bps), min(counts), max(counts)])

        manifest_path = args.output_root / "ncrna.haploid.manifest.tsv"
        with manifest_path.open("w", encoding="utf-8", newline="\n") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["label", "species", "haplotype", "tRNA_input", "rfam_tblout_input", "genome", "gff3", "seqid_check"])
            for result in results:
                writer.writerow([result["label"], result["species"], result["haplotype"], result["trna_input"], result["rfam_input"], result["genome"], result["gff3"], result["seqid_check"]])
    print(f"summary: {summary_path}")


if __name__ == "__main__":
    main()
