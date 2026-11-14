#!/usr/bin/env python3
"""Build gene -> KO -> KEGG pathway mapping tables for PZI, PSP and PPA.

KEGG gene IDs were transferred from Swiss-Prot cross-references
(annotation/functional/<species>/swissprot_kegg.tsv). Local snapshots from the KEGG REST API resolve:

  link/ko/<org>      KEGG gene -> KO                (one cached request per organism)
  link/pathway/ko    KO -> reference pathway map    (cached once)
  list/pathway       pathway map -> name            (cached once)
  get/br:ko00001     pathway map -> BRITE class A/B (cached once)

Outputs one long table per species:
  annotation/functional/<species>/swissprot_kegg_pathway.tsv
with columns protein_id, kegg_id, ko_id, pathway_id, pathway_name,
pathway_class_a, pathway_class_b.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--root", type=Path, required=True)
args = parser.parse_args()
ROOT = args.root
ANN = ROOT / "annotation" / "functional"
CACHE = ANN / "kegg_cache"
CACHE.mkdir(parents=True, exist_ok=True)
SPECIES = ["pzi", "psp", "ppa"]
REST = "https://rest.kegg.jp"


def fetch(url: str, cache_name: str) -> Path:
    path = CACHE / cache_name
    if not path.is_file():
        raise FileNotFoundError(f"Download {url} to {path} before running this script")
    return path


def parse_gene_ko(path: Path) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        gene, ko = parts[0].strip(), parts[1].strip().removeprefix("ko:")
        if gene and ko:
            out.setdefault(gene, []).append(ko)
    return out


def parse_link_pathway_ko(path: Path) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        ko = parts[0].strip().removeprefix("ko:")
        pathway = parts[1].strip().removeprefix("path:")
        if not pathway.startswith("map"):
            continue
        out.setdefault(ko, []).append(pathway)
    return out


def parse_pathway_names(path: Path) -> dict[str, str]:
    out = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split("\t", 1)
        if len(parts) == 2:
            out[parts[0].strip()] = parts[1].strip()
    return out


def parse_brite_classes(path: Path) -> dict[str, tuple[str, str]]:
    out: dict[str, tuple[str, str]] = {}
    class_a = class_b = ""
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("A") and len(line) > 1 and not line[1].isspace():
            rest = line[1:].strip()
            tokens = rest.split(None, 1)
            class_a = tokens[1] if len(tokens) > 1 else rest
        elif line.startswith("B") and len(line) > 1 and line[1].isspace():
            rest = line[1:].strip()
            tokens = rest.split(None, 1)
            class_b = tokens[1] if len(tokens) > 1 else rest
        elif line.startswith("C"):
            paths = re.findall(r"PATH:([A-Za-z]+\d{5})", line)
            for pathway in paths:
                map_id = re.sub(r"^ko", "map", pathway)
                if map_id not in out:
                    out[map_id] = (class_a, class_b)
    return out


def main() -> None:
    kegg = {}
    for species in SPECIES:
        df = pd.read_csv(
            ANN / species / "swissprot_kegg.tsv",
            sep="\t",
            usecols=["protein_id", "kegg_id"],
            dtype=str,
        ).dropna()
        kegg[species] = df
        print(
            f"{species}: {len(df)} rows, {df['protein_id'].nunique()} proteins",
            flush=True,
        )

    orgs = sorted(
        set().union(*(set(df["kegg_id"].str.split(":").str[0]) for df in kegg.values()))
    )
    print(f"organisms to resolve: {len(orgs)}", flush=True)

    gene2ko: dict[str, list[str]] = {}
    for i, org in enumerate(orgs, start=1):
        path = fetch(f"{REST}/link/ko/{org}", f"link_ko_{org}.tsv")
        if path is None:
            continue
        by_gene = parse_gene_ko(path)
        gene2ko.update(by_gene)
        print(f"  [{i}/{len(orgs)}] {org}: {len(by_gene)} genes with KO", flush=True)

    ko2path_path = fetch(f"{REST}/link/pathway/ko", "link_pathway_ko.tsv")
    ko2path = parse_link_pathway_ko(ko2path_path) if ko2path_path else {}
    names_path = fetch(f"{REST}/list/pathway", "list_pathway.tsv")
    names = parse_pathway_names(names_path) if names_path else {}
    brite_path = fetch(f"{REST}/get/br:ko00001", "br_ko00001.txt")
    classes = parse_brite_classes(brite_path) if brite_path else {}
    print(
        f"KO->pathway links: {len(ko2path)} KOs | pathways: {len(names)} | "
        f"BRITE classes: {len(classes)}",
        flush=True,
    )

    summary = []
    for species in SPECIES:
        rows = []
        for protein_id, kegg_id in kegg[species].itertuples(index=False):
            for ko in gene2ko.get(kegg_id, []):
                for pathway in ko2path.get(ko, []):
                    class_a, class_b = classes.get(pathway, ("", ""))
                    rows.append(
                        {
                            "protein_id": protein_id,
                            "kegg_id": kegg_id,
                            "ko_id": ko,
                            "pathway_id": pathway,
                            "pathway_name": names.get(pathway, pathway),
                            "pathway_class_a": class_a,
                            "pathway_class_b": class_b,
                        }
                    )
        table = pd.DataFrame(rows).drop_duplicates()
        out = ANN / species / "swissprot_kegg_pathway.tsv"
        table.to_csv(out, sep="\t", index=False)
        summary.append(
            {
                "species": species,
                "proteins": int(kegg[species]["protein_id"].nunique()),
                "proteins_with_pathway": int(table["protein_id"].nunique()),
                "ko_ids": int(table["ko_id"].nunique()),
                "pathways": int(table["pathway_id"].nunique()),
            }
        )
        print(f"{species}: saved {out} ({len(table)} rows)", flush=True)

    print(pd.DataFrame(summary).to_string(index=False))


if __name__ == "__main__":
    main()
