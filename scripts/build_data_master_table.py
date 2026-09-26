#!/usr/bin/env python3
import csv
from pathlib import Path


def label_from_path(path: Path):
    return path.name


def file_types(d: Path):
    types = []
    for p in d.iterdir():
        if p.is_file():
            n = p.name.lower()
            if n.endswith((".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")):
                types.append("fa")
            if n.endswith((".gff", ".gff3", ".gff.gz", ".gff3.gz")):
                types.append("gff")
            if n.endswith((".pep", ".faa", ".pep.gz", ".faa.gz")):
                types.append("pep")
            if n.endswith((".cds", ".cds.gz")):
                types.append("cds")
    return sorted(set(types))


def main():
    root = Path("data")
    details = root / "details"
    pub = details / "plabipd_publication_scored_download_targets.tsv"
    pub_map = {}
    if pub.is_file():
        with pub.open(encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                pub_map[row["label"]] = row

    report_map = {}
    report_sources = [
        details / "_ncbi_reports",
        details / "_ncbi_reports_expanded",
        details / "_ncbi_reports_plabipd",
    ]
    for repdir in report_sources:
        if not repdir.is_dir():
            continue
        for rp in repdir.glob("*_dataset_report.json"):
            report_map[rp.name.replace("_dataset_report.json", "")] = str(rp)

    rows = []
    for d in sorted(p for p in root.iterdir() if p.is_dir() and p.name != "details" and not p.name.startswith("_")):
        if d.name == "ppe_pytozome":
            continue
        ft = ",".join(file_types(d))
        pubrow = pub_map.get(d.name, {})
        rows.append({
            "label": d.name,
            "species": pubrow.get("species", ""),
            "source": pubrow.get("source", ""),
            "doi": pubrow.get("doi", ""),
            "publication_year": pubrow.get("year", ""),
            "downloaded_from": pubrow.get("selected_ftp", pubrow.get("selected_ftp", "")),
            "files_present": ft,
            "annotation_status": pubrow.get("selected_annotation", ""),
            "notes": pubrow.get("title", "")[:120],
            "report_file": report_map.get(d.name, ""),
        })

    out = root / "README_DATA_TABLE.tsv"
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["label", "species", "source", "doi", "publication_year", "downloaded_from", "files_present", "annotation_status", "notes", "report_file"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(out)


if __name__ == "__main__":
    main()
