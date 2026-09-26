#!/usr/bin/env python3
"""Transfer Swiss-Prot metadata to existing DIAMOND best-hit results.

The query proteins and DIAMOND alignments are not changed.  UniProt flat-file
records are indexed by accession, then GO, KEGG, KO, EC, InterPro and selected
functional fields are joined to each query protein through its Swiss-Prot
best hit.  This is homology transfer, not de novo functional prediction.
"""

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def clean_field(value: str) -> str:
    return re.sub(r"\s+", " ", value.strip())


def parse_full_name(de_lines):
    names = []
    for line in de_lines:
        for match in re.finditer(r"Full=([^;]+)", line):
            names.append(clean_field(match.group(1)))
    return "; ".join(dict.fromkeys(names))


def parse_gene_names(gn_lines):
    names = []
    for line in gn_lines:
        for match in re.finditer(r"(?:Name|Synonyms)=([^;]+)", line):
            names.append(clean_field(match.group(1)))
    return "; ".join(dict.fromkeys(names))


def parse_record(lines):
    accessions = []
    de_lines = []
    gn_lines = []
    function_lines = []
    pathway_lines = []
    cc_section = ""
    go = []
    kegg = []
    ko = []
    ec = []
    interpro = []
    for line in lines:
        tag = line[:2]
        payload = line[5:].rstrip("\n") if len(line) > 5 else ""
        if tag == "AC":
            accessions.extend(x.strip().rstrip(";") for x in payload.split(";") if x.strip())
        elif tag == "DE":
            de_lines.append(payload)
        elif tag == "GN":
            gn_lines.append(payload)
        elif tag == "CC":
            if payload.startswith("-!- FUNCTION:"):
                cc_section = "function"
                function_lines.append(payload.split(":", 1)[1].strip())
            elif payload.startswith("-!- PATHWAY:"):
                cc_section = "pathway"
                pathway_lines.append(payload.split(":", 1)[1].strip())
            elif cc_section == "function" and payload.strip() and not payload.startswith("-!- "):
                function_lines.append(payload.strip())
            elif cc_section == "pathway" and payload.strip() and not payload.startswith("-!- "):
                pathway_lines.append(payload.strip())
            elif payload.startswith("-!- "):
                cc_section = ""
        elif tag == "DR":
            fields = [clean_field(x) for x in payload.split(";")]
            if not fields:
                continue
            db = fields[0]
            values = [x for x in fields[1:] if x]
            if db == "GO" and values:
                go_id = values[0]
                aspect = ""
                description = ""
                evidence = ""
                if len(values) > 1 and ":" in values[1]:
                    aspect, description = values[1].split(":", 1)
                if len(values) > 2 and ":" in values[2]:
                    evidence = values[2].split(":", 1)[0]
                go.append((go_id, aspect, description, evidence, "; ".join(values[1:])))
            elif db == "KEGG" and values and values[0] not in {"-", "."}:
                kegg.append(values[0])
            elif db == "KO" and values and values[0] not in {"-", "."}:
                ko.append(values[0])
            elif db == "EC" and values and values[0] not in {"-", "."}:
                ec.append(values[0])
            elif db == "InterPro" and values and values[0] not in {"-", "."}:
                interpro.append((values[0], values[1] if len(values) > 1 else ""))
    if not accessions:
        return None
    record = {
        "accessions": list(dict.fromkeys(accessions)),
        "description": parse_full_name(de_lines),
        "gene_names": parse_gene_names(gn_lines),
        "function": clean_field(" ".join(function_lines)),
        "pathway": clean_field(" ".join(pathway_lines)),
        "go": list(dict.fromkeys(go)),
        "kegg": list(dict.fromkeys(kegg)),
        "ko": list(dict.fromkeys(ko)),
        "ec": list(dict.fromkeys(ec)),
        "interpro": list(dict.fromkeys(interpro)),
    }
    return record


def load_dat(path: Path, wanted=None):
    records = {}
    current = []
    wanted = set(wanted or ())
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("//"):
                record = parse_record(current)
                record_keys = {accession.split("-")[0] for accession in record["accessions"]} if record else set()
                if record and (not wanted or wanted.intersection(record_keys)):
                    for accession in record["accessions"]:
                        key = accession.split("-")[0]
                        if not wanted or key in wanted:
                            records[key] = record
                current = []
            else:
                current.append(line)
    if current:
        record = parse_record(current)
        record_keys = {accession.split("-")[0] for accession in record["accessions"]} if record else set()
        if record and (not wanted or wanted.intersection(record_keys)):
            for accession in record["accessions"]:
                key = accession.split("-")[0]
                if not wanted or key in wanted:
                    records[key] = record
    return records


def parse_accession(value: str) -> str:
    parts = value.split("|")
    if len(parts) >= 3 and parts[0] in {"sp", "tr"}:
        return parts[1].split("-")[0]
    return value.split("-")[0]


def load_hits(path: Path):
    hits = {}
    duplicate_queries = 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 13:
                raise ValueError(f"{path}:{line_number}: expected >=13 DIAMOND columns, found {len(fields)}")
            if fields[0].lower() == "qseqid":
                continue
            query, subject = fields[0], fields[1]
            if query in hits:
                duplicate_queries += 1
                continue
            hits[query] = {
                "query": query,
                "subject": subject,
                "accession": parse_accession(subject),
                "pident": fields[2],
                "alignment_length": fields[3],
                "evalue": fields[10],
                "bitscore": fields[11],
                "title": "\t".join(fields[12:]),
            }
    return hits, duplicate_queries


def join_values(values):
    return ";".join(dict.fromkeys(str(value) for value in values if str(value)))


def write_outputs(out_dir: Path, species: str, hits, metadata, dat_path: Path, diamond_path: Path, duplicate_queries: int):
    out_dir.mkdir(parents=True, exist_ok=True)
    annotation_rows = []
    go_rows = []
    kegg_rows = []
    ko_rows = []
    ec_rows = []
    interpro_rows = []
    missing = 0
    for query in sorted(hits):
        hit = hits[query]
        record = metadata.get(hit["accession"])
        if record is None:
            missing += 1
            record = {"description": "", "gene_names": "", "function": "", "pathway": "", "go": [], "kegg": [], "ko": [], "ec": [], "interpro": []}
        base = [species, hit["query"], hit["subject"], hit["accession"], hit["pident"], hit["alignment_length"], hit["evalue"], hit["bitscore"]]
        for go_id, aspect, description, evidence, raw in record["go"]:
            go_rows.append(base + [go_id, aspect, description, evidence, raw, "Swiss-Prot homology transfer"])
        for value in record["kegg"]:
            kegg_rows.append(base + [value, "Swiss-Prot cross-reference"])
        for value in record["ko"]:
            ko_rows.append(base + [value, "Swiss-Prot cross-reference"])
        for value in record["ec"]:
            ec_rows.append(base + [value, "Swiss-Prot cross-reference"])
        for value, description in record["interpro"]:
            interpro_rows.append(base + [value, description, "Swiss-Prot cross-reference"])
        annotation_rows.append(base + [record["description"], record["gene_names"], record["function"], record["pathway"], join_values(x[0] for x in record["go"]), join_values(record["kegg"]), join_values(record["ko"]), join_values(record["ec"]), join_values(x[0] for x in record["interpro"])])

    def write(path, header, rows):
        with path.open("w", encoding="utf-8", newline="\n") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            writer.writerows(rows)

    base_header = ["species", "protein_id", "swissprot_entry", "swissprot_accession", "pident", "alignment_length", "evalue", "bitscore"]
    write(out_dir / "swissprot_go.tsv", base_header + ["go_id", "go_aspect", "go_description", "go_evidence", "go_raw", "transfer_method"], go_rows)
    write(out_dir / "swissprot_kegg.tsv", base_header + ["kegg_id", "mapping_source"], kegg_rows)
    write(out_dir / "swissprot_ko.tsv", base_header + ["ko_id", "mapping_source"], ko_rows)
    write(out_dir / "swissprot_ec.tsv", base_header + ["ec_number", "mapping_source"], ec_rows)
    write(out_dir / "swissprot_interpro.tsv", base_header + ["interpro_id", "interpro_description", "mapping_source"], interpro_rows)
    write(out_dir / "swissprot_functional_annotation.tsv", base_header + ["protein_description", "gene_names", "function", "pathway", "go_terms", "kegg_ids", "ko_ids", "ec_numbers", "interpro_ids"], annotation_rows)

    with (out_dir / "swissprot_metadata_manifest.txt").open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"Species: {species}\n")
        handle.write(f"DIAMOND input: {diamond_path}\n")
        handle.write(f"Swiss-Prot metadata: {dat_path}\n")
        handle.write("Join key: Swiss-Prot primary accession parsed from sp|ACCESSION|ENTRY_NAME\n")
        handle.write("GO transfer: source terms copied from the Swiss-Prot hit; target proteins are marked as homology transfer\n")
        handle.write("KEGG/KO/EC/InterPro: Swiss-Prot cross-reference fields\n")
        handle.write(f"Query proteins: {len(hits)}\n")
        handle.write(f"Swiss-Prot records indexed: {len(metadata)} accession keys\n")
        handle.write(f"Queries without metadata match: {missing}\n")
        handle.write(f"Duplicate query rows ignored after first best-hit row: {duplicate_queries}\n")
    return len(annotation_rows), len(go_rows), len(kegg_rows), len(ko_rows), len(ec_rows), len(interpro_rows), missing, duplicate_queries


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", required=True)
    parser.add_argument("--diamond", type=Path, required=True)
    parser.add_argument("--swissprot-dat", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    if not args.diamond.is_file() or args.diamond.stat().st_size == 0:
        raise SystemExit(f"missing DIAMOND result: {args.diamond}")
    if not args.swissprot_dat.is_file() or args.swissprot_dat.stat().st_size == 0:
        raise SystemExit(f"missing Swiss-Prot metadata: {args.swissprot_dat}")
    hits, duplicate_queries = load_hits(args.diamond)
    if not hits:
        raise SystemExit(f"no DIAMOND hits found: {args.diamond}")
    metadata = load_dat(args.swissprot_dat, {hit["accession"] for hit in hits.values()})
    counts = write_outputs(args.out_dir, args.species, hits, metadata, args.swissprot_dat, args.diamond, duplicate_queries)
    print(f"{args.species}: proteins={counts[0]} GO={counts[1]} KEGG={counts[2]} KO={counts[3]} EC={counts[4]} InterPro={counts[5]} metadata_missing={counts[6]}")


if __name__ == "__main__":
    main()
