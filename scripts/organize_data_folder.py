#!/usr/bin/env python3
import csv
import shutil
from pathlib import Path


ROOT_SPECIES = {
    "crpi", "cydo", "ejap", "fve", "masi", "mbac", "mdom", "msyl", "par", "par_anno",
    "pav", "pbre", "pcam", "pcom", "pdav", "pdu", "pjam", "pmal", "pmir", "pmong",
    "pmu", "ppa", "ppe", "psal", "pserru", "pseu", "pspe", "rch", "rgla"
}

KEEP_EXT = {
    ".fa", ".fna", ".fasta", ".gff", ".gff3", ".pep", ".faa", ".cds",
    ".fa.gz", ".fna.gz", ".fasta.gz", ".gff.gz", ".gff3.gz", ".faa.gz", ".pep.gz", ".cds.gz",
    ".txt", ".tsv", ".md"
}

DETAIL_DIRS = {
    "_ncbi_reports", "_ncbi_reports_expanded", "_ncbi_reports_plabipd",
}

DETAIL_FILES = {
    "assembly_summary_genbank.txt", "assembly_summary_refseq.txt", "comparative_download_plan.tsv",
    "download_report.tsv", "downloaded_genome_audit.tsv", "expanded_candidate_species.tsv",
    "jgi_main.js", "jgi_ppersica_page.html", "ncbi_chromosome_candidates.tsv", "ncbi_direct_download_report.tsv",
    "ncbi_direct_ftp_plan.tsv", "ncbi_direct_ftp_plan_expanded.tsv", "ncbi_direct_ftp_plan_plabipd.tsv",
    "ncbi_direct_ftp_plan_plabipd_genome_first.tsv", "plabipd_candidate_publications.tsv", "plabipd_genomes_angio.json",
    "plabipd_publication_scored_download_targets.tsv", "plabipd_pubplant.js", "plabipd_pubplant_cladogram1.html",
    "plabipd_rosaceae_candidates.tsv", "ppe_test.zip", "prunus_focused_expansion_plan.tsv", "README_COMPARATIVE_GENOMES.md",
    "refseq_rosaceae_chromosome_annotated_candidates.tsv", "selected_chromosome_genomes.tsv", "selected_expanded_genomes.tsv",
    "selected_plabipd_best_ncbi_genomes.tsv", "selected_plabipd_genomes.tsv", "test_prunus_persica_dataset_report.json",
    "tmp_download_pav.tsv", "transcriptome_accessions_template.tsv", "_test_malus_report.json",
}


def is_detail_candidate(path: Path) -> bool:
    if path.name in DETAIL_FILES:
        return True
    if path.name.startswith("plabipd_") and path.suffix in {".json", ".js", ".html"}:
        return True
    if path.name.startswith("ncbi_") and path.suffix in {".tsv", ".txt", ".json"}:
        return True
    if path.name.startswith("selected_") and path.suffix == ".tsv":
        return True
    if path.name.startswith("assembly_summary_") and path.suffix == ".txt":
        return True
    if path.name.startswith("_test_") or path.name.startswith("tmp_"):
        return True
    if path.suffix == ".zip":
        return True
    return False


def safe_move(src: Path, dst: Path):
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        if dst.is_file():
            dst.unlink()
        else:
            shutil.rmtree(dst)
    shutil.move(str(src), str(dst))


def canonical_species_file(path: Path) -> bool:
    name = path.name.lower()
    return any(name.endswith(ext) for ext in [".fa", ".fna", ".fasta", ".gff", ".gff3", ".pep", ".faa", ".cds", ".fa.gz", ".fna.gz", ".fasta.gz", ".gff.gz", ".gff3.gz", ".pep.gz", ".faa.gz", ".cds.gz"])


def main():
    root = Path("data")
    details = root / "details"
    details.mkdir(exist_ok=True)

    # Move detail files at top level.
    for path in list(root.iterdir()):
        if path.name in ROOT_SPECIES or path.name == "details":
            continue
        if path.is_dir() and path.name in DETAIL_DIRS:
            safe_move(path, details / path.name)
            continue
        if path.is_file() and is_detail_candidate(path):
            safe_move(path, details / path.name)

    # Flatten species dirs so only canonical files remain in top-level species dirs.
    for species_dir in [p for p in root.iterdir() if p.is_dir() and p.name in ROOT_SPECIES]:
        species_details = details / species_dir.name
        species_details.mkdir(exist_ok=True)
        for path in list(species_dir.rglob("*")):
            if path.is_dir():
                continue
            if canonical_species_file(path):
                continue
            rel = path.relative_to(species_dir)
            safe_move(path, species_details / rel)
        # Keep only the canonical package files in the species directory.
        for path in list(species_dir.iterdir()):
            if path.is_file() and not canonical_species_file(path):
                safe_move(path, species_details / path.name)
            elif path.is_dir() and path.name not in {"deep", "keep"}:
                safe_move(path, species_details / path.name)

    # Move any non-species directories inside data to details.
    for path in list(root.iterdir()):
        if path.is_dir() and path.name not in ROOT_SPECIES and path.name != "details":
            safe_move(path, details / path.name)

    print(root / "details")


if __name__ == "__main__":
    main()
