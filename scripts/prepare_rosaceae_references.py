#!/usr/bin/env python3
"""Curate Rosaceae reference proteins for inflorescence-family trees.

The raw UniProt TSV exports are kept in the pipeline reference directory.
This script selects complete, explicitly annotated proteins, removes exact
sequence duplicates within each tree family, and writes a stable FASTA plus
metadata table consumed by ``run_inflorescence_family_pipeline.py``.
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REF_DIR = ROOT / "analysis" / "inflorescence_family_pipeline" / "reference"
OUT_FASTA = REF_DIR / "rosaceae_reference.fa"
OUT_META = REF_DIR / "rosaceae_reference_metadata.tsv"

# Accessions are selected from the downloaded raw UniProt exports.  The
# labels are conservative tree anchors; they are not claims about PZI gene
# function.  DAM proteins are retained as the Rosaceae SVP/DAM branch.
SELECTIONS = [
    # LFY/FLO
    ("A0A2P6P4I3", "C_LFY_FLO", "LFY", "LFY/FLO"),
    ("A0A6P5TNW4", "C_LFY_FLO", "LFY", "LFY/FLO"),
    ("A2I8Z8", "C_LFY_FLO", "LFY", "LFY/FLO"),
    ("Q8W3Y8", "C_LFY_FLO", "FLO", "LFY/FLO"),
    ("Q8W3Y7", "C_LFY_FLO", "FLO", "LFY/FLO"),
    # PEBP
    ("Q75QW8", "PEBP", "FT", "FT-like"),
    ("B9VJ16", "PEBP", "FT2", "FT-like"),
    ("E2D2C4", "PEBP", "FT3", "FT-like"),
    ("A0A1X9QDN9", "PEBP", "FT", "FT-like"),
    ("I1XUU5", "PEBP", "TFL1", "TFL1-like"),
    ("A0A6P5RKF1", "PEBP", "MFT", "MFT-like"),
    # MADS: inflorescence and floral-organ subfamilies
    ("D6PW23", "MADS", "CAL", "AP1/CAL"),
    ("A0A6P5SRA4", "MADS", "SEP1", "SEP"),
    ("A0A6P5TEN7", "MADS", "AG", "AG/SHP/STK"),
    ("Q0GPY8", "MADS", "SHP", "AG/SHP/STK"),
    ("A5Z0S5", "MADS", "STK", "AG/SHP/STK"),
    ("D0VYC0", "MADS", "SOC1", "SOC1"),
    ("A0A6P5TEX7", "MADS", "SOC1", "SOC1"),
    ("A0A6P5SBJ9", "MADS", "SVP", "SVP/DAM"),
    ("A0A8F9WIG2", "MADS", "DAM4", "SVP/DAM"),
    ("A0A0D5CT77", "MADS", "DAM3", "SVP/DAM"),
    ("A0A1B2LP65", "MADS", "DAM2", "SVP/DAM"),
    ("A0A1B2LP93", "MADS", "DAM1", "SVP/DAM"),
    ("A0A6P5TEV8", "MADS", "AGL6", "AGL6"),
    ("A0A251Q9R2", "MADS", "AGL15", "other-MADS"),
    ("A0A6P5TH17", "MADS", "AGL18", "other-MADS"),
    ("A0A0U5A402", "MADS", "AP3", "other-MADS"),
    ("A0A0U4ZUQ7", "MADS", "PI", "other-MADS"),
    ("A0A6P5RB90", "MADS", "AGL62", "other-MADS"),
]


def load_rows() -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    for path in sorted(REF_DIR.glob("rosaceae_*.tsv")):
        with path.open(encoding="utf-8", errors="ignore", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                acc = (row.get("Entry") or "").strip()
                seq = re.sub(r"[^A-Za-z*]", "", row.get("Sequence") or "").replace("*", "")
                if not acc or not seq:
                    continue
                # Multiple exports contain the same accession.  Keep the
                # longest non-empty sequence and retain the first metadata.
                old = rows.get(acc)
                if old is None or len(seq) > len(old["sequence"]):
                    row = dict(row)
                    row["sequence"] = seq
                    rows[acc] = row
    return rows


def species_name(raw: str) -> str:
    raw = (raw or "").strip()
    return raw.split(" (")[0]


def safe_token(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")


def main() -> int:
    rows = load_rows()
    selected: list[dict[str, str]] = []
    missing: list[str] = []
    seen: set[tuple[str, str]] = set()
    for accession, family, symbol, group in SELECTIONS:
        row = rows.get(accession)
        if row is None:
            missing.append(accession)
            continue
        sequence = row["sequence"]
        # Exact isoform duplicates are not independent phylogenetic evidence.
        key = (family, sequence)
        if key in seen:
            continue
        seen.add(key)
        species = species_name(row.get("Organism", ""))
        species_code = {
            "Prunus persica": "Ppers",
            "Prunus avium": "Pavium",
            "Malus domestica": "Mdom",
            "Rosa chinensis": "Rchin",
            "Fragaria vesca": "Fves",
        }.get(species, safe_token(species))
        reference_id = f"ROS_{species_code}_{symbol}_{accession}"
        selected.append({
            "protein_id": reference_id,
            "species": species,
            "accession": accession,
            "gene_symbol": symbol,
            "tree_family": family,
            "subfamily_group": group,
            "curation_source": "UniProt_raw_TSV",
            "source_entry_name": row.get("Entry Name", ""),
            "protein_name": row.get("Protein names", ""),
            "length": str(len(sequence)),
            "sequence": sequence,
        })
    if missing:
        raise SystemExit("Missing selected accessions: " + ",".join(missing))
    selected.sort(key=lambda r: (r["tree_family"], r["subfamily_group"], r["species"], r["gene_symbol"], r["accession"]))
    fields = [
        "protein_id", "species", "accession", "gene_symbol", "tree_family",
        "subfamily_group", "curation_source", "source_entry_name",
        "protein_name", "length",
    ]
    REF_DIR.mkdir(parents=True, exist_ok=True)
    with OUT_META.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in selected:
            writer.writerow({k: row[k] for k in fields})
    with OUT_FASTA.open("w", encoding="utf-8") as fh:
        for row in selected:
            fh.write(f">{row['protein_id']} species={row['species']} accession={row['accession']} gene_symbol={row['gene_symbol']} tree_family={row['tree_family']} subfamily_group={row['subfamily_group']}\n")
            seq = row["sequence"]
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")
    print(f"wrote {len(selected)} unique Rosaceae references")
    for family in ("C_LFY_FLO", "PEBP", "MADS"):
        print(f"{family}: {sum(r['tree_family'] == family for r in selected)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
