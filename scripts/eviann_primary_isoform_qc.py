#!/usr/bin/env python3
"""Build primary-isoform EviAnn annotation sets and summarize BUSCO results."""

import argparse
import csv
import hashlib
import re
import shutil
from pathlib import Path


HAPS = ["hapA1", "hapA2", "hapB1", "hapB2"]
COMBOS = {
    "A1A2": ["hapA1", "hapA2"],
    "B1B2": ["hapB1", "hapB2"],
    "A1B1": ["hapA1", "hapB1"],
    "A1B2": ["hapA1", "hapB2"],
    "A2B1": ["hapA2", "hapB1"],
    "A2B2": ["hapA2", "hapB2"],
    "A1A2B1B2": ["hapA1", "hapA2", "hapB1", "hapB2"],
}
PACKAGE_FILES = {
    "genome": ("genome", ".fa"),
    "gff": ("gff", ".gff3"),
    "cds": ("cds", ".cds.fa"),
    "pep": ("pep", ".pep.fa"),
}


def parse_attrs(text):
    attrs = []
    for part in text.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
            attrs.append((key, value))
        else:
            attrs.append((part, ""))
    return attrs


def attrs_to_dict(attrs):
    return {key: value for key, value in attrs}


def format_attrs(attrs):
    parts = []
    for key, value in attrs:
        if value == "":
            parts.append(key)
        else:
            parts.append("%s=%s" % (key, value))
    return ";".join(parts)


def prefixed(hap, value):
    prefix = "pzi_%s_" % hap
    if value.startswith(prefix):
        return value
    return prefix + value


def package_input(package_dir, hap, kind):
    """Resolve an input from the renamed subdirectory package or old flat layout."""
    subdir, suffix = PACKAGE_FILES[kind]
    label = "pzi_%s" % hap
    candidates = [
        Path(package_dir) / subdir / (label + suffix),
        Path(package_dir) / (label + suffix),
    ]
    for path in candidates:
        if path.is_file() and path.stat().st_size > 0:
            return path
    raise FileNotFoundError(
        "No usable %s input for %s; checked: %s"
        % (kind, hap, ", ".join(str(path) for path in candidates))
    )


def transform_attrs(attrs, hap, feature, selected_gene_ids, selected_mrna_ids):
    out = []
    original_id = ""
    for key, value in attrs:
        if key == "ID":
            original_id = value
            out.append((key, prefixed(hap, value)))
        elif key == "Parent":
            parents = []
            for parent in value.split(","):
                if parent in selected_gene_ids or parent in selected_mrna_ids:
                    parents.append(prefixed(hap, parent))
                else:
                    parents.append(parent)
            out.append((key, ",".join(parents)))
        elif key == "geneID" and value in selected_gene_ids:
            out.append((key, prefixed(hap, value)))
        else:
            out.append((key, value))
    if feature == "mRNA":
        out.append(("primary_isoform", "true"))
    if original_id:
        out.append(("original_id", original_id))
    return out


def read_fasta(path):
    records = {}
    order = []
    current = None
    desc = ""
    chunks = []
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if current is not None:
                    records[current] = (desc, "".join(chunks))
                header = line[1:].rstrip("\n")
                current = header.split()[0]
                desc = header[len(current):].strip()
                order.append(current)
                chunks = []
            else:
                chunks.append(line.strip())
        if current is not None:
            records[current] = (desc, "".join(chunks))
    return order, records


def write_fasta(records, path):
    with Path(path).open("w", encoding="utf-8") as out:
        for seq_id, desc, seq in records:
            header = ">%s" % seq_id
            if desc:
                header += " " + desc
            out.write(header + "\n")
            for idx in range(0, len(seq), 60):
                out.write(seq[idx:idx + 60] + "\n")


def parse_gff(path):
    records = []
    gene_to_mrnas = {}
    mrna_to_gene = {}
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw.strip():
                continue
            if raw.startswith("#"):
                records.append({"kind": "comment", "raw": raw.rstrip("\n")})
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = parse_attrs(fields[8])
            ad = attrs_to_dict(attrs)
            feature = fields[2]
            rec = {"kind": "feature", "fields": fields, "attrs": attrs, "attr_dict": ad}
            records.append(rec)
            if feature == "gene" and "ID" in ad:
                gene_to_mrnas.setdefault(ad["ID"], [])
            elif feature == "mRNA" and "ID" in ad and "Parent" in ad:
                parent = ad["Parent"].split(",")[0]
                gene_to_mrnas.setdefault(parent, []).append(ad["ID"])
                mrna_to_gene[ad["ID"]] = parent
    return records, gene_to_mrnas, mrna_to_gene


def load_qc(path):
    if not Path(path).is_file():
        return {}
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return {row["mrna_id"]: row for row in csv.DictReader(handle, delimiter="\t")}


def numeric(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def evidence_rank(evidence):
    evidence = (evidence or "").lower()
    if evidence == "complete":
        return 4
    if "transcript" in evidence and "protein" in evidence:
        return 3
    if "transcript" in evidence:
        return 2
    if "protein" in evidence:
        return 1
    return 0


def candidate_score(hap, mrna_id, mrna_attrs, qc_rows, cds_records, pep_records):
    qc = qc_rows.get(mrna_id, {})
    flags = set(filter(None, qc.get("flags", "").split(",")))
    pep_len = int(numeric(qc.get("protein_len"), len(pep_records.get(mrna_id, ("", ""))[1]) if mrna_id in pep_records else 0))
    cds_len = int(numeric(qc.get("cds_len"), len(cds_records.get(mrna_id, ("", ""))[1]) if mrna_id in cds_records else 0))
    has_protein = 1 if mrna_id in pep_records and pep_len > 0 else 0
    pass_flag = 1 if "PASS" in flags else 0
    no_internal_stop = 1 if qc.get("internal_stop", "0") in ("0", "", "False", "false") else 0
    multiple_of_three = 1 if qc.get("cds_not_multiple_of_3", "0") in ("0", "", "False", "false") else 0
    has_start = 1 if qc.get("has_start_M", "0") in ("1", "True", "true") else 0
    has_stop = 1 if qc.get("has_terminal_stop", "0") in ("1", "True", "true") else 0
    not_short = 1 if "short_cds" not in flags else 0
    low_te = 1 if "high_gene_te_overlap" not in flags and "high_cds_te_overlap" not in flags else 0
    ev_rank = evidence_rank(mrna_attrs.get("Evidence", ""))
    # Prefer valid protein-coding models, then QC PASS, then complete ORFs and evidence support, then length.
    return (
        has_protein,
        pass_flag,
        no_internal_stop,
        multiple_of_three,
        has_start + has_stop,
        not_short,
        low_te,
        ev_rank,
        pep_len,
        cds_len,
        -len(mrna_id),
    )


def choose_primary(hap, gene_to_mrnas, mrna_attr_by_id, qc_rows, cds_records, pep_records):
    selected = {}
    rows = []
    skipped = []
    for gene_id in sorted(gene_to_mrnas):
        candidates = gene_to_mrnas[gene_id]
        if not candidates:
            skipped.append({"gene_id": gene_id, "reason": "no_mrna"})
            continue
        ranked = []
        for mrna_id in candidates:
            score = candidate_score(hap, mrna_id, mrna_attr_by_id.get(mrna_id, {}), qc_rows, cds_records, pep_records)
            ranked.append((score, mrna_id))
        ranked.sort(reverse=True)
        best_score, best_mrna = ranked[0]
        if best_score[0] == 0:
            skipped.append({"gene_id": gene_id, "reason": "no_protein_for_any_isoform"})
            continue
        selected[gene_id] = best_mrna
        qc = qc_rows.get(best_mrna, {})
        rows.append({
            "hap_set": hap,
            "gene_id": gene_id,
            "primary_mrna_id": best_mrna,
            "candidate_isoforms": str(len(candidates)),
            "protein_len": str(int(numeric(qc.get("protein_len"), len(pep_records.get(best_mrna, ("", ""))[1])))),
            "cds_len": str(int(numeric(qc.get("cds_len"), len(cds_records.get(best_mrna, ("", ""))[1])))),
            "flags": qc.get("flags", ""),
            "evidence": mrna_attr_by_id.get(best_mrna, {}).get("Evidence", ""),
            "selection_score": ",".join(map(str, best_score)),
        })
    return selected, rows, skipped


def write_tsv(path, rows, fieldnames):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def count_gff(path):
    counts = {}
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) >= 3:
                counts[fields[2]] = counts.get(fields[2], 0) + 1
    return counts


def fasta_count_bp(path):
    count = 0
    bp = 0
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if raw.startswith(">"):
                count += 1
            else:
                bp += len(raw.strip())
    return count, bp


def copy_or_concat(paths, out_path, skip_gff_header=False):
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with Path(out_path).open("w", encoding="utf-8") as out:
        wrote_header = False
        for path in paths:
            with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    if skip_gff_header and line.startswith("##gff-version"):
                        if not wrote_header:
                            out.write("##gff-version 3\n")
                            wrote_header = True
                        continue
                    out.write(line)


def build_one_hap(hap, package_dir, eviann_root, out_dir):
    gff = package_input(package_dir, hap, "gff")
    cds = package_input(package_dir, hap, "cds")
    pep = package_input(package_dir, hap, "pep")
    genome = package_input(package_dir, hap, "genome")
    qc = eviann_root / hap / "postprocess" / ("pzi_%s.eviann.gene_model_qc.tsv" % hap)
    records, gene_to_mrnas, mrna_to_gene = parse_gff(gff)
    _, cds_records = read_fasta(cds)
    _, pep_records = read_fasta(pep)
    qc_rows = load_qc(qc)
    mrna_attr_by_id = {}
    for rec in records:
        if rec.get("kind") == "feature" and rec["fields"][2] == "mRNA":
            ad = rec["attr_dict"]
            if "ID" in ad:
                mrna_attr_by_id[ad["ID"]] = ad
    selected, selection_rows, skipped = choose_primary(hap, gene_to_mrnas, mrna_attr_by_id, qc_rows, cds_records, pep_records)
    selected_gene_ids = set(selected)
    selected_mrna_ids = set(selected.values())

    hap_prefix = "pzi_%s.primary" % hap
    out_gff = out_dir / ("%s.gff3" % hap_prefix)
    out_cds = out_dir / ("%s.cds.fa" % hap_prefix)
    out_pep = out_dir / ("%s.pep.fa" % hap_prefix)
    out_fa = out_dir / ("pzi_%s.fa" % hap)
    shutil.copy2(genome, out_fa)

    with out_gff.open("w", encoding="utf-8") as out:
        out.write("##gff-version 3\n")
        for rec in records:
            if rec.get("kind") != "feature":
                continue
            fields = list(rec["fields"])
            feature = fields[2]
            ad = rec["attr_dict"]
            write = False
            if feature == "gene" and ad.get("ID") in selected_gene_ids:
                write = True
            elif feature == "mRNA" and ad.get("ID") in selected_mrna_ids:
                write = True
            elif "Parent" in ad:
                parents = ad["Parent"].split(",")
                write = any(parent in selected_mrna_ids for parent in parents)
            if write:
                fields[8] = format_attrs(transform_attrs(rec["attrs"], hap, feature, selected_gene_ids, selected_mrna_ids))
                out.write("\t".join(fields) + "\n")

    cds_out_records = []
    pep_out_records = []
    for row in selection_rows:
        old = row["primary_mrna_id"]
        new = prefixed(hap, old)
        desc = "orig_id=%s gene=%s flags=%s evidence=%s" % (
            old, row["gene_id"], row["flags"] or "NA", row["evidence"] or "NA"
        )
        if old in cds_records:
            cds_out_records.append((new, desc, cds_records[old][1]))
        if old in pep_records:
            pep_out_records.append((new, desc, pep_records[old][1]))
    write_fasta(cds_out_records, out_cds)
    write_fasta(pep_out_records, out_pep)

    write_tsv(
        out_dir / ("%s.primary_isoform_selection.tsv" % hap_prefix),
        selection_rows,
        ["hap_set", "gene_id", "primary_mrna_id", "candidate_isoforms", "protein_len", "cds_len", "flags", "evidence", "selection_score"],
    )
    write_tsv(out_dir / ("%s.primary_skipped_genes.tsv" % hap_prefix), skipped, ["gene_id", "reason"])

    gff_counts = count_gff(out_gff)
    cds_count, cds_bp = fasta_count_bp(out_cds)
    pep_count, pep_aa = fasta_count_bp(out_pep)
    genome_count, genome_bp = fasta_count_bp(out_fa)
    return {
        "label": hap,
        "type": "hap",
        "members": hap,
        "genome_fasta": str(out_fa),
        "gff3": str(out_gff),
        "cds_fasta": str(out_cds),
        "pep_fasta": str(out_pep),
        "genome_sequences": str(genome_count),
        "genome_bp": str(genome_bp),
        "gff_gene": str(gff_counts.get("gene", 0)),
        "gff_mRNA": str(gff_counts.get("mRNA", 0)),
        "gff_CDS": str(gff_counts.get("CDS", 0)),
        "cds_records": str(cds_count),
        "cds_bp": str(cds_bp),
        "pep_records": str(pep_count),
        "pep_aa": str(pep_aa),
        "skipped_genes": str(len(skipped)),
        "status": "pass" if gff_counts.get("gene", 0) == cds_count == pep_count else "review",
    }


def md5sum(path):
    h = hashlib.md5()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def build_primary(args):
    package_dir = Path(args.package_dir)
    eviann_root = Path(args.eviann_root)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = []
    for hap in HAPS:
        summary.append(build_one_hap(hap, package_dir, eviann_root, out_dir))

    for label, members in COMBOS.items():
        combo_prefix = "pzi_%s.primary" % label
        member_prefixes = ["pzi_%s.primary" % hap for hap in members]
        out_fa = out_dir / ("%s.fa" % combo_prefix)
        out_gff = out_dir / ("%s.gff3" % combo_prefix)
        out_cds = out_dir / ("%s.cds.fa" % combo_prefix)
        out_pep = out_dir / ("%s.pep.fa" % combo_prefix)
        copy_or_concat([out_dir / ("pzi_%s.fa" % hap) for hap in members], out_fa)
        copy_or_concat([out_dir / ("%s.gff3" % p) for p in member_prefixes], out_gff, skip_gff_header=True)
        copy_or_concat([out_dir / ("%s.cds.fa" % p) for p in member_prefixes], out_cds)
        copy_or_concat([out_dir / ("%s.pep.fa" % p) for p in member_prefixes], out_pep)
        gff_counts = count_gff(out_gff)
        cds_count, cds_bp = fasta_count_bp(out_cds)
        pep_count, pep_aa = fasta_count_bp(out_pep)
        genome_count, genome_bp = fasta_count_bp(out_fa)
        summary.append({
            "label": label,
            "type": "combo",
            "members": ",".join(members),
            "genome_fasta": str(out_fa),
            "gff3": str(out_gff),
            "cds_fasta": str(out_cds),
            "pep_fasta": str(out_pep),
            "genome_sequences": str(genome_count),
            "genome_bp": str(genome_bp),
            "gff_gene": str(gff_counts.get("gene", 0)),
            "gff_mRNA": str(gff_counts.get("mRNA", 0)),
            "gff_CDS": str(gff_counts.get("CDS", 0)),
            "cds_records": str(cds_count),
            "cds_bp": str(cds_bp),
            "pep_records": str(pep_count),
            "pep_aa": str(pep_aa),
            "skipped_genes": "0",
            "status": "pass" if gff_counts.get("gene", 0) == cds_count == pep_count else "review",
        })

    fields = [
        "label", "type", "members", "genome_sequences", "genome_bp",
        "gff_gene", "gff_mRNA", "gff_CDS", "cds_records", "cds_bp",
        "pep_records", "pep_aa", "skipped_genes", "status",
        "genome_fasta", "gff3", "cds_fasta", "pep_fasta",
    ]
    write_tsv(out_dir / "primary_summary.tsv", summary, fields)
    with (out_dir / "CHECKSUMS.md5").open("w", encoding="utf-8") as out:
        for path in sorted(out_dir.glob("*.fa")) + sorted(out_dir.glob("*.gff3")):
            out.write("%s  %s\n" % (md5sum(path), path.name))
    (out_dir / "README_PRIMARY_EVIANN.md").write_text(
        "\n".join([
            "Primary isoform EviAnn annotation package",
            "=========================================",
            "",
            "Each gene is represented by one primary protein-coding isoform.",
            "Selection prioritizes protein presence, QC PASS status, complete ORF checks,",
            "low TE-overlap flags, EviAnn evidence type, and then protein/CDS length.",
            "",
            "Use *.primary.pep.fa for BUSCO, OrthoFinder, gene family, and RBH analyses.",
            "The full EviAnn annotation should still be retained for transcript-level review.",
            "",
        ]),
        encoding="utf-8",
    )
    print("Primary EviAnn package:", out_dir)


def parse_busco_summary(path):
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    row = {
        "busco_path": str(path),
        "complete_pct": "",
        "single_pct": "",
        "duplicated_pct": "",
        "fragmented_pct": "",
        "missing_pct": "",
        "n_buscos": "",
        "status": "parsed_without_compact_line",
    }
    match = re.search(r"C:([0-9.]+)%\[S:([0-9.]+)%,D:([0-9.]+)%\],F:([0-9.]+)%,M:([0-9.]+)%,n:(\d+)", text)
    if match:
        row.update({
            "complete_pct": match.group(1),
            "single_pct": match.group(2),
            "duplicated_pct": match.group(3),
            "fragmented_pct": match.group(4),
            "missing_pct": match.group(5),
            "n_buscos": match.group(6),
            "status": "pass",
        })
    return row


def summarize_busco(args):
    busco_dir = Path(args.busco_dir)
    rows = []
    for path in sorted(busco_dir.rglob("short_summary*.txt")):
        label = ""
        for candidate in ["A1A2B1B2", "hapA1", "hapA2", "hapB1", "hapB2", "A1A2", "B1B2", "A1B1", "A1B2", "A2B1", "A2B2"]:
            if candidate in str(path):
                label = candidate
                break
        row = parse_busco_summary(path)
        row["label"] = label or path.parent.name
        rows.append(row)
    write_tsv(
        args.out,
        rows,
        ["label", "status", "complete_pct", "single_pct", "duplicated_pct", "fragmented_pct", "missing_pct", "n_buscos", "busco_path"],
    )
    print("BUSCO summary:", args.out)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    build = subparsers.add_parser("build-primary")
    build.add_argument("--package-dir", required=True)
    build.add_argument("--eviann-root", required=True)
    build.add_argument("--out-dir", required=True)
    busco = subparsers.add_parser("summarize-busco")
    busco.add_argument("--busco-dir", required=True)
    busco.add_argument("--out", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.command == "build-primary":
        build_primary(args)
    elif args.command == "summarize-busco":
        summarize_busco(args)


if __name__ == "__main__":
    main()
