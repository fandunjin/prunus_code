#!/usr/bin/env python3
"""Prepare canonical downstream comparative-genomics inputs.

The workflow intentionally keeps one protein/CDS per gene. Downloaded or official
primary files are copied directly. NCBI-like all-isoform files are reduced to the
longest valid CDS isoform per gene using the GFF gene/transcript/CDS structure.
"""

import argparse
import csv
import gzip
import re
import shutil
from collections import defaultdict
from pathlib import Path


TRANSCRIPT_FEATURES = {"mRNA", "transcript", "lnc_RNA", "rRNA", "tRNA"}
PRIMARY_STRATEGIES = {"downloaded_primary", "official_primary"}
NUC = set("ACGTNacgtn")


def open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fieldnames):
    ensure_dir(Path(path).parent)
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def parse_attrs(text):
    attrs = {}
    for part in text.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
        elif " " in part:
            key, value = part.split(" ", 1)
            value = value.strip('"')
        else:
            key, value = part, ""
        attrs[key] = value
    return attrs


def fasta_iter(path):
    current = None
    desc = ""
    chunks = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if current is not None:
                    yield current, desc, "".join(chunks)
                header = line[1:].strip()
                current = header.split()[0]
                desc = header[len(current):].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        if current is not None:
            yield current, desc, "".join(chunks)


def read_fasta(path):
    records = {}
    for seq_id, desc, seq in fasta_iter(path):
        records[seq_id] = (desc, seq)
    return records


def write_fasta(path, records):
    ensure_dir(Path(path).parent)
    with Path(path).open("w", encoding="utf-8") as out:
        for seq_id, desc, seq in records:
            header = ">" + seq_id
            if desc:
                header += " " + desc
            out.write(header + "\n")
            for idx in range(0, len(seq), 60):
                out.write(seq[idx:idx + 60] + "\n")


def copy_plain(src, dst, force=False):
    ensure_dir(Path(dst).parent)
    src = Path(src)
    dst = Path(dst)
    if dst.is_file() and dst.stat().st_size > 0 and not force:
        return
    # Some JCVI jobs stage authoritative inputs as symlinks in canonical/.  A
    # forced refresh must replace the link itself, otherwise copy2() follows it
    # and can try to overwrite the source package.
    if force and dst.is_symlink():
        dst.unlink()
    if src.suffix == ".gz":
        with gzip.open(src, "rb") as inp, dst.open("wb") as out:
            shutil.copyfileobj(inp, out)
    else:
        shutil.copy2(src, dst)


def fasta_count(path):
    records = 0
    length = 0
    for _, _, seq in fasta_iter(path):
        records += 1
        length += len(seq)
    return records, length


def normalize_token(value):
    return value.split()[0].strip()


def candidate_ids(seq_id, desc):
    ids = {seq_id}
    text = seq_id + " " + desc
    for pattern in [
        r"(?:transcript_id|protein_id|gene|gene_id|locus_tag|Parent|ID)=([^;\s]+)",
        r"\[protein_id=([^\]]+)\]",
        r"\[gene=([^\]]+)\]",
        r"\[locus_tag=([^\]]+)\]",
    ]:
        for hit in re.findall(pattern, text):
            ids.add(normalize_token(hit))
    if "." in seq_id:
        ids.add(seq_id.rsplit(".", 1)[0])
    return ids


def build_fasta_index(records):
    index = {}
    for seq_id, (desc, seq) in records.items():
        for key in candidate_ids(seq_id, desc):
            index.setdefault(key, seq_id)
    return index


def parse_gff_models(path):
    genes = {}
    transcript_to_gene = {}
    transcript_attrs = {}
    transcript_extra_ids = defaultdict(set)
    cds_len = defaultdict(int)
    cds_lines_by_transcript = defaultdict(list)
    gene_order = []

    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid, source, feature, start, end, score, strand, phase, attr_text = fields
            attrs = parse_attrs(attr_text)
            feature_id = attrs.get("ID", "")
            parents = [p for p in attrs.get("Parent", "").split(",") if p]

            if feature == "gene":
                gene_id = feature_id or attrs.get("Name") or attrs.get("gene_id")
                if gene_id and gene_id not in genes:
                    genes[gene_id] = []
                    gene_order.append(gene_id)
                continue

            if feature in TRANSCRIPT_FEATURES:
                tx_id = feature_id or attrs.get("transcript_id")
                parent = parents[0] if parents else attrs.get("gene") or attrs.get("gene_id")
                if tx_id and parent:
                    transcript_to_gene[tx_id] = parent
                    transcript_attrs[tx_id] = attrs
                    genes.setdefault(parent, []).append(tx_id)
                    if parent not in gene_order:
                        gene_order.append(parent)
                continue

            if feature == "CDS":
                tx_parents = parents if parents else [attrs.get("transcript_id", "")]
                span = abs(int(end) - int(start)) + 1
                for tx_id in tx_parents:
                    if not tx_id:
                        continue
                    cds_len[tx_id] += span
                    cds_lines_by_transcript[tx_id].append(fields)
                    for attr_key in ["protein_id", "transcript_id", "Name", "gene", "gene_id", "locus_tag", "Parent"]:
                        value = attrs.get(attr_key)
                        if value:
                            for item in value.split(","):
                                transcript_extra_ids[tx_id].add(item)
                    if tx_id not in transcript_attrs:
                        transcript_attrs[tx_id] = attrs

    # Some GFFs attach CDS directly to genes.
    for tx_id in list(cds_len):
        if tx_id not in transcript_to_gene:
            transcript_to_gene[tx_id] = tx_id
            genes.setdefault(tx_id, []).append(tx_id)
            if tx_id not in gene_order:
                gene_order.append(tx_id)

    return gene_order, genes, transcript_attrs, transcript_extra_ids, cds_len, cds_lines_by_transcript


def attr_score(attrs):
    text = ";".join("%s=%s" % (k, v) for k, v in attrs.items()).lower()
    score = 0
    for token in ["primary", "canonical", "representative", "longest"]:
        if token in text:
            score += 100
    return score


def is_valid_cds(seq):
    clean = seq.upper().replace("-", "")
    if not clean or len(clean) % 3 != 0:
        return False
    return all(base in NUC for base in clean)


def is_valid_protein(seq):
    clean = seq.strip().rstrip("*")
    return "*" not in clean and len(clean) > 0


def select_primary_isoforms(gene_order, genes, transcript_attrs, transcript_extra_ids, cds_len, cds_records, pep_records):
    cds_index = build_fasta_index(cds_records)
    pep_index = build_fasta_index(pep_records)
    selected = []
    skipped = []

    for gene_id in gene_order:
        candidates = []
        for tx_id in genes.get(gene_id, []):
            cds_key = None
            pep_key = None
            attrs = transcript_attrs.get(tx_id, {})
            extra_ids = [tx_id] + sorted(transcript_extra_ids.get(tx_id, set()))
            for attr_key in ["ID", "protein_id", "transcript_id", "Name", "gene", "gene_id", "locus_tag"]:
                value = attrs.get(attr_key)
                if value and value not in extra_ids:
                    extra_ids.append(value)
            # NCBI CDS FASTA headers often key the CDS by protein_id. Require the
            # CDS and protein to share the same identifier when possible; this
            # avoids crossing isoforms under the same gene.
            for value in extra_ids:
                if cds_index.get(value) and pep_index.get(value):
                    cds_key = cds_index[value]
                    pep_key = pep_index[value]
                    break
            if not cds_key:
                for value in extra_ids:
                    if cds_index.get(value):
                        cds_key = cds_index[value]
                        break
            if not pep_key:
                for value in extra_ids:
                    if pep_index.get(value):
                        pep_key = pep_index[value]
                        break
            if not cds_key or not pep_key:
                continue
            cds_desc, cds_seq = cds_records[cds_key]
            pep_desc, pep_seq = pep_records[pep_key]
            if not is_valid_cds(cds_seq) or not is_valid_protein(pep_seq):
                continue
            candidates.append({
                "gene_id": gene_id,
                "transcript_id": tx_id,
                "cds_key": cds_key,
                "pep_key": pep_key,
                "cds_len": cds_len.get(tx_id, len(cds_seq)),
                "protein_len": len(pep_seq.rstrip("*")),
                "score": attr_score(attrs),
                "cds_desc": cds_desc,
                "cds_seq": cds_seq,
                "pep_desc": pep_desc,
                "pep_seq": pep_seq.rstrip("*"),
            })
        if not candidates:
            skipped.append({"gene_id": gene_id, "reason": "no_valid_cds_and_protein"})
            continue
        best = sorted(candidates, key=lambda row: (row["score"], row["cds_len"], row["protein_len"], row["transcript_id"]), reverse=True)[0]
        selected.append(best)
    return selected, skipped


def filter_gff(path, selected_tx_ids, out_path):
    selected_tx_ids = set(selected_tx_ids)
    selected_gene_ids = set()
    rows = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                rows.append(line)
                continue
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = parse_attrs(fields[8])
            feature_id = attrs.get("ID", "")
            parents = set(p for p in attrs.get("Parent", "").split(",") if p)
            if fields[2] in TRANSCRIPT_FEATURES and feature_id in selected_tx_ids:
                selected_gene_ids.update(parents)
                rows.append(line)
            elif fields[2] == "CDS" and parents.intersection(selected_tx_ids):
                rows.append(line)
            elif fields[2] == "gene":
                rows.append(("__GENE__", line, feature_id))
            elif parents.intersection(selected_tx_ids):
                rows.append(line)

    with Path(out_path).open("w", encoding="utf-8") as out:
        for row in rows:
            if isinstance(row, tuple):
                _, line, feature_id = row
                if feature_id in selected_gene_ids:
                    out.write(line)
            else:
                out.write(row)


def resolve(root, value):
    if not value:
        return None
    path = Path(value)
    if path.is_absolute():
        return path
    return root / path


def prepare_primary_species(root, row, out_root, force=False):
    label = row["label"]
    out_dir = out_root / "canonical" / label
    ensure_dir(out_dir)
    sources = {
        "genome": resolve(root, row.get("genome_source", "")),
        "gff": resolve(root, row.get("gff_source", "")),
        "cds": resolve(root, row.get("cds_source", "")),
        "pep": resolve(root, row.get("pep_source", "")),
    }
    for key in ["genome", "gff", "cds", "pep"]:
        if sources[key] and sources[key].is_file():
            suffix = {"genome": "genome.fa", "gff": "gene.gff3", "cds": "cds.fa", "pep": "pep.fa"}[key]
            copy_plain(sources[key], out_dir / suffix, force=force)

    strategy = row["isoform_strategy"]
    primary_pep = out_dir / "primary.pep.fa"
    primary_cds = out_dir / "primary.cds.fa"
    primary_gff = out_dir / "primary.gff3"

    if strategy in PRIMARY_STRATEGIES:
        if (out_dir / "pep.fa").is_file():
            copy_plain(out_dir / "pep.fa", primary_pep, force=force)
        if (out_dir / "cds.fa").is_file():
            copy_plain(out_dir / "cds.fa", primary_cds, force=force)
        if (out_dir / "gene.gff3").is_file():
            copy_plain(out_dir / "gene.gff3", primary_gff, force=force)
        selected_count = fasta_count(primary_pep)[0] if primary_pep.is_file() else 0
        skipped_count = 0
    elif strategy == "longest_CDS_from_GFF":
        required = [out_dir / "gene.gff3", out_dir / "cds.fa", out_dir / "pep.fa"]
        missing = [str(p) for p in required if not p.is_file()]
        if missing:
            raise SystemExit("ERROR: %s missing required files for isoform extraction: %s" % (label, ", ".join(missing)))
        if primary_pep.is_file() and primary_cds.is_file() and primary_gff.is_file() and not force:
            selected_count = fasta_count(primary_pep)[0]
            skipped_path = out_dir / "primary_isoform_skipped.tsv"
            skipped_count = max(0, sum(1 for _ in skipped_path.open("r", encoding="utf-8")) - 1) if skipped_path.is_file() else 0
        else:
            gene_order, genes, tx_attrs, tx_extra_ids, cds_len, _ = parse_gff_models(out_dir / "gene.gff3")
            cds_records = read_fasta(out_dir / "cds.fa")
            pep_records = read_fasta(out_dir / "pep.fa")
            selected, skipped = select_primary_isoforms(gene_order, genes, tx_attrs, tx_extra_ids, cds_len, cds_records, pep_records)
            write_fasta(primary_cds, [("%s|%s" % (x["gene_id"], x["transcript_id"]), x["cds_desc"], x["cds_seq"]) for x in selected])
            write_fasta(primary_pep, [("%s|%s" % (x["gene_id"], x["transcript_id"]), x["pep_desc"], x["pep_seq"]) for x in selected])
            filter_gff(out_dir / "gene.gff3", [x["transcript_id"] for x in selected], primary_gff)
            write_tsv(
                out_dir / "primary_isoform_selection.tsv",
                selected,
                ["gene_id", "transcript_id", "cds_key", "pep_key", "cds_len", "protein_len", "score"],
            )
            write_tsv(out_dir / "primary_isoform_skipped.tsv", skipped, ["gene_id", "reason"])
            selected_count = len(selected)
            skipped_count = len(skipped)
    elif strategy in {"needs_annotation", "transcriptome_longest_ORF"}:
        selected_count = 0
        skipped_count = 0
    else:
        raise SystemExit("ERROR: unknown isoform_strategy for %s: %s" % (label, strategy))

    genome_records, genome_len = fasta_count(out_dir / "genome.fa") if (out_dir / "genome.fa").is_file() else (0, 0)
    pep_records_count, pep_len = fasta_count(primary_pep) if primary_pep.is_file() else (0, 0)
    cds_records_count, cds_len_total = fasta_count(primary_cds) if primary_cds.is_file() else (0, 0)
    return {
        "label": label,
        "role": row["role"],
        "tier": row["tier"],
        "data_class": row["data_class"],
        "isoform_strategy": strategy,
        "include_main_orthofinder": row["include_main_orthofinder"],
        "include_sensitivity": row["include_sensitivity"],
        "include_synteny": row["include_synteny"],
        "canonical_dir": str(out_dir),
        "genome_records": genome_records,
        "genome_bp": genome_len,
        "primary_pep_records": pep_records_count,
        "primary_pep_aa": pep_len,
        "primary_cds_records": cds_records_count,
        "primary_cds_bp": cds_len_total,
        "selected_isoforms": selected_count,
        "skipped_genes": skipped_count,
    }


def link_or_copy_orthofinder(summary_rows, out_root, set_name, include_field):
    pep_dir = out_root / "orthofinder" / set_name / "pep"
    ensure_dir(pep_dir)
    rows = []
    for row in summary_rows:
        if row.get(include_field) != "yes":
            continue
        src = Path(row["canonical_dir"]) / "primary.pep.fa"
        if not src.is_file() or src.stat().st_size == 0:
            continue
        dst = pep_dir / ("%s.pep.fa" % row["label"])
        shutil.copy2(src, dst)
        n_records, n_aa = fasta_count(dst)
        rows.append({
            "label": row["label"],
            "source": str(src),
            "orthofinder_file": str(dst),
            "records": n_records,
            "aa": n_aa,
        })
    write_tsv(out_root / "orthofinder" / set_name / "pep_manifest.tsv", rows, ["label", "source", "orthofinder_file", "records", "aa"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".", help="Repository/workflow root.")
    parser.add_argument("--manifest", default="downstream_analysis/config/species_manifest.tsv")
    parser.add_argument("--out-dir", default="analysis")
    parser.add_argument("--labels", default="", help="Optional comma-separated label subset.")
    parser.add_argument("--force", action="store_true", help="Rebuild primary outputs even if they already exist.")
    parser.add_argument("--skip-orthofinder-sets", action="store_true", help="Only prepare canonical species folders.")
    args = parser.parse_args()

    root = Path(args.root).resolve()
    manifest = resolve(root, args.manifest)
    out_root = resolve(root, args.out_dir)
    ensure_dir(out_root)

    summary = []
    requested = set(label.strip() for label in args.labels.split(",") if label.strip())
    for row in read_tsv(manifest):
        if requested and row["label"] not in requested:
            continue
        summary.append(prepare_primary_species(root, row, out_root, force=args.force))

    fields = [
        "label", "role", "tier", "data_class", "isoform_strategy",
        "include_main_orthofinder", "include_sensitivity", "include_synteny",
        "canonical_dir", "genome_records", "genome_bp", "primary_pep_records",
        "primary_pep_aa", "primary_cds_records", "primary_cds_bp",
        "selected_isoforms", "skipped_genes",
    ]
    if requested:
        manifest_name = "canonical_manifest_%s.tsv" % "_".join(sorted(requested))
    else:
        manifest_name = "canonical_manifest.tsv"
    write_tsv(out_root / manifest_name, summary, fields)
    if not args.skip_orthofinder_sets:
        link_or_copy_orthofinder(summary, out_root, "main", "include_main_orthofinder")
        link_or_copy_orthofinder(summary, out_root, "sensitivity", "include_sensitivity")
    print(out_root / manifest_name)


if __name__ == "__main__":
    main()
