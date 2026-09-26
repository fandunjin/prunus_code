#!/usr/bin/env python3
import json
import sys
from pathlib import Path


def main():
    if len(sys.argv) < 3:
        raise SystemExit("Usage: summarize_ncbi_reports.py OUT_TSV REPORT_JSON...")
    out = Path(sys.argv[1])
    rows = []
    for path_s in sys.argv[2:]:
        path = Path(path_s)
        label = path.name.replace("_dataset_report.json", "")
        with path.open(encoding="utf-8") as handle:
            data = json.load(handle)
        for report in data.get("reports", []):
            info = report.get("assembly_info", {})
            stats = report.get("assembly_stats", {})
            ann = report.get("annotation_info", {})
            if info.get("assembly_level") not in ("Chromosome", "Complete Genome"):
                continue
            org = report.get("organism", {})
            rows.append([
                label,
                report.get("accession", ""),
                report.get("source_database", ""),
                org.get("organism_name", ""),
                org.get("infraspecific_names", {}).get("cultivar", ""),
                info.get("assembly_name", ""),
                info.get("assembly_level", ""),
                info.get("refseq_category", ""),
                str(stats.get("total_number_of_chromosomes", "")),
                str(stats.get("total_sequence_length", "")),
                str(stats.get("number_of_scaffolds", "")),
                str(stats.get("scaffold_n50", "")),
                ann.get("status", "") if ann else "",
                ann.get("name", "") if ann else "",
            ])
    rows.sort(key=lambda r: (r[0], r[12] == "", r[7] != "reference genome", r[1]))
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as handle:
        handle.write("label\taccession\tsource\tspecies\tcultivar\tassembly_name\tassembly_level\trefseq_category\tchromosomes\ttotal_length\tscaffolds\tscaffold_n50\tannotation_status\tannotation_name\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
