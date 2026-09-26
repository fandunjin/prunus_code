#!/usr/bin/env python3
import csv
import json
import sys
from pathlib import Path


def score_report(report):
    info = report.get("assembly_info", {})
    stats = report.get("assembly_stats", {})
    ann = report.get("annotation_info")
    accession = report.get("accession", "")
    score = 0
    if info.get("assembly_status") == "current":
        score += 1_000_000
    if info.get("assembly_level") == "Complete Genome":
        score += 500_000
    elif info.get("assembly_level") == "Chromosome":
        score += 450_000
    elif info.get("assembly_level") == "Scaffold":
        score += 100_000
    if accession.startswith("GCF_") or report.get("source_database") == "SOURCE_DATABASE_REFSEQ":
        score += 150_000
    if info.get("refseq_category") == "reference genome":
        score += 100_000
    if ann:
        score += 50_000
    if ann and ann.get("status") == "Full annotation":
        score += 30_000
    try:
        score += min(int(stats.get("scaffold_n50", 0)) // 100_000, 5000)
    except Exception:
        pass
    return score


def ftp_path(report):
    info = report.get("assembly_info", {})
    if info.get("ftp_path"):
        return info["ftp_path"].rstrip("/")
    accession = report.get("accession", "")
    asm = info.get("assembly_name", "")
    if not accession or not asm:
        return ""
    digits = accession.split("_", 1)[1].split(".", 1)[0]
    if len(digits) < 9:
        return ""
    return (
        f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{accession[:3]}/"
        f"{digits[:3]}/{digits[3:6]}/{digits[6:9]}/{accession}_{asm}"
    )


def load_candidate_meta(path):
    meta = {}
    if not path or not Path(path).exists():
        return meta
    with open(path, encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            meta[row["label"]] = row
    return meta


def main():
    if len(sys.argv) < 5:
        raise SystemExit(
            "Usage: select_ncbi_genomes_from_reports.py REPORT_DIR CANDIDATE_TSV SELECTED_TSV PLAN_TSV"
        )
    report_dir = Path(sys.argv[1])
    candidate_tsv = Path(sys.argv[2])
    selected_tsv = Path(sys.argv[3])
    plan_tsv = Path(sys.argv[4])
    meta = load_candidate_meta(candidate_tsv)

    rows = []
    for report_path in sorted(report_dir.glob("*_dataset_report.json")):
        label = report_path.name.replace("_dataset_report.json", "")
        data = json.loads(report_path.read_text(encoding="utf-8"))
        reports = data.get("reports", [])
        current = [
            r for r in reports
            if r.get("assembly_info", {}).get("assembly_status") == "current"
            and r.get("assembly_info", {}).get("assembly_level") in ("Complete Genome", "Chromosome", "Scaffold")
        ]
        chosen = sorted(current, key=score_report, reverse=True)[0] if current else None
        m = meta.get(label, {})
        if not chosen:
            rows.append({
                "label": label,
                "species": m.get("species", ""),
                "accession": "",
                "assembly_name": "",
                "assembly_level": "",
                "annotation": "no_ncbi_genome_found",
                "ftp_path": "",
                "group": m.get("group", ""),
                "priority": m.get("priority", ""),
                "notes": m.get("notes", ""),
            })
            continue
        ann = chosen.get("annotation_info") or {}
        rows.append({
            "label": label,
            "species": chosen.get("organism", {}).get("organism_name", m.get("species", "")),
            "accession": chosen.get("accession", ""),
            "assembly_name": chosen.get("assembly_info", {}).get("assembly_name", ""),
            "assembly_level": chosen.get("assembly_info", {}).get("assembly_level", ""),
            "annotation": ann.get("status") or ("submitted_annotation" if ann else "genome_only"),
            "protein_coding_gene_count": (
                ann.get("stats", {}).get("gene_counts", {}).get("protein_coding", "")
                if ann else ""
            ),
            "scaffold_n50": chosen.get("assembly_stats", {}).get("scaffold_n50", ""),
            "ftp_path": ftp_path(chosen),
            "group": m.get("group", ""),
            "priority": m.get("priority", ""),
            "notes": m.get("notes", ""),
        })

    selected_tsv.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "label", "species", "accession", "assembly_name", "assembly_level",
        "annotation", "protein_coding_gene_count", "scaffold_n50", "ftp_path",
        "group", "priority", "notes",
    ]
    rows.sort(key=lambda r: (0 if r["accession"] else 1, int(r.get("priority") or 9), r.get("group", ""), r["species"]))
    with selected_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    plan_fields = ["label", "accession", "assembly_name", "species", "need_annotation", "ftp_path", "source_note"]
    with plan_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=plan_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for r in rows:
            if not r["accession"]:
                continue
            need_annotation = "yes" if r["annotation"] in ("Full annotation", "submitted_annotation") else "no"
            writer.writerow({
                "label": r["label"],
                "accession": r["accession"],
                "assembly_name": r["assembly_name"],
                "species": r["species"],
                "need_annotation": need_annotation,
                "ftp_path": r["ftp_path"],
                "source_note": f"{r['assembly_level']}; {r['annotation']}; {r.get('notes','')}",
            })

    print(selected_tsv)
    print(plan_tsv)


if __name__ == "__main__":
    main()
