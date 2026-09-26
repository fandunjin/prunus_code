#!/usr/bin/env python3
import json
import sys
from pathlib import Path


def score_report(report):
    info = report.get("assembly_info", {})
    stats = report.get("assembly_stats", {})
    ann = report.get("annotation_info", {})
    accession = report.get("accession", "")
    source = report.get("source_database", "")
    level = info.get("assembly_level", "")
    refcat = info.get("refseq_category", "")
    status = info.get("assembly_status", "")
    score = 0
    if status == "current":
        score += 1000
    if level == "Chromosome":
        score += 500
    elif level == "Complete Genome":
        score += 300
    if accession.startswith("GCF_") or source == "SOURCE_DATABASE_REFSEQ":
        score += 200
    if refcat == "reference genome":
        score += 300
    elif refcat == "representative genome":
        score += 200
    if ann:
        score += 100
    if ann.get("status") == "Full annotation":
        score += 100
    try:
        score += min(int(stats.get("scaffold_n50", 0)) // 1000000, 100)
    except Exception:
        pass
    try:
        score += min(int(stats.get("total_number_of_chromosomes", 0)), 50)
    except Exception:
        pass
    return score


def main():
    if len(sys.argv) < 3:
        raise SystemExit("Usage: select_ncbi_chromosome_genomes.py OUT_TSV REPORT_JSON...")

    out_tsv = Path(sys.argv[1])
    rows = []
    for path_s in sys.argv[2:]:
        path = Path(path_s)
        label = path.stem.replace("_dataset_report", "")
        with path.open(encoding="utf-8") as handle:
            data = json.load(handle)
        reports = data.get("reports", [])
        chrom = [
            r for r in reports
            if r.get("assembly_info", {}).get("assembly_level") in ("Chromosome", "Complete Genome")
            and r.get("assembly_info", {}).get("assembly_status") == "current"
        ]
        candidates = chrom if chrom else reports
        if not candidates:
            rows.append([label, "", "", "", "", "", "", "", "NO_REPORTS"])
            continue
        best = sorted(candidates, key=score_report, reverse=True)[0]
        org = best.get("organism", {})
        info = best.get("assembly_info", {})
        stats = best.get("assembly_stats", {})
        ann = best.get("annotation_info", {})
        cultivar = org.get("infraspecific_names", {}).get("cultivar", "")
        rows.append([
            label,
            best.get("accession", ""),
            org.get("organism_name", ""),
            cultivar,
            info.get("assembly_name", ""),
            info.get("assembly_level", ""),
            info.get("refseq_category", ""),
            str(stats.get("scaffold_n50", "")),
            ann.get("status", "NO_ANNOTATION") if ann else "NO_ANNOTATION",
        ])

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", encoding="utf-8", newline="") as out:
        out.write("label\taccession\tspecies\tcultivar\tassembly_name\tassembly_level\trefseq_category\tscaffold_n50\tannotation_status\n")
        for row in rows:
            out.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
