#!/usr/bin/env python3
import csv
import re
import sys
from pathlib import Path


JOURNAL_WEIGHTS = [
    (re.compile(r"\bScience\b|Nature Genetics|Nat Genet|Nature Communications|Nat Commun|Genome Biology|Genome Biol", re.I), 100),
    (re.compile(r"Scientific Data|Sci Data|GigaScience|Gigascience|Horticulture Research|Hortic Res", re.I), 90),
    (re.compile(r"Plant Journal|Plant J|J Integr Plant Biol|DNA Research|DNA Res|New Phytologist|New Phytol", re.I), 80),
    (re.compile(r"Current Biology|Curr Biol|BMC Biology|BMC Biol", re.I), 75),
    (re.compile(r"Frontiers|Front Genet|Plants \\(Basel\\)|IJMS|Int J Mol Sci|PeerJ|G3", re.I), 55),
]


def journal_score(source):
    for pattern, score in JOURNAL_WEIGHTS:
        if pattern.search(source or ""):
            return score
    return 50


def title_score(title):
    t = (title or "").lower()
    score = 0
    if "telomere-to-telomere" in t or "t2t" in t:
        score += 20
    if "chromosome" in t:
        score += 15
    if "genome assembly and annotation" in t or "assembly and annotation" in t:
        score += 15
    if "pan-genome" in t or "population genomics" in t:
        score += 5
    if "transcriptomic" in t and "genomic" in t:
        score += 5
    return score


def main():
    if len(sys.argv) < 4:
        raise SystemExit("Usage: score_plabipd_publications.py PUBLICATIONS_TSV SELECTED_GENOMES_TSV OUT_TSV")
    pub_tsv = Path(sys.argv[1])
    selected_tsv = Path(sys.argv[2])
    out_tsv = Path(sys.argv[3])

    selected = {}
    with selected_tsv.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            selected[row["label"]] = row

    rows = []
    with pub_tsv.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            label = row["label"]
            sel = selected.get(label, {})
            try:
                year = int(row.get("year") or 0)
            except ValueError:
                year = 0
            score = journal_score(row.get("source", "")) + title_score(row.get("title", ""))
            if year >= 2024:
                score += 10
            elif year >= 2020:
                score += 5
            if sel.get("assembly_level") == "Complete Genome":
                score += 15
            elif sel.get("assembly_level") == "Chromosome":
                score += 10
            if sel.get("annotation") == "Full annotation":
                score += 15
            elif sel.get("annotation") == "submitted_annotation":
                score += 10
            row.update({
                "selected_accession": sel.get("accession", ""),
                "selected_assembly": sel.get("assembly_name", ""),
                "selected_level": sel.get("assembly_level", ""),
                "selected_annotation": sel.get("annotation", ""),
                "selected_ftp": sel.get("ftp_path", ""),
                "publication_priority_score": str(score),
            })
            rows.append(row)

    fields = [
        "label", "species", "group", "priority", "genome_size_mbp", "ref_index",
        "year", "doi", "source", "title", "authorship",
        "selected_accession", "selected_assembly", "selected_level",
        "selected_annotation", "selected_ftp", "publication_priority_score",
    ]
    rows.sort(key=lambda r: (-int(r["publication_priority_score"]), int(r.get("priority") or 9), r["species"]))
    with out_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(out_tsv)


if __name__ == "__main__":
    main()
