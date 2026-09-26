#!/usr/bin/env python3
"""Extract single-copy protein/CDS FASTA files for MAFFT and codon alignment.

This is intentionally only a glue step: OrthoFinder identifies single-copy
orthogroups; MAFFT aligns proteins; trimAl -backtrans builds frame-aware
trimmed codon alignments. The script only maps OrthoFinder gene IDs back to
canonical primary protein/CDS records and renames each single-copy sequence to
its species label.
"""

import argparse
import csv
import re
from collections import Counter
from pathlib import Path


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path):
    order = []
    records = {}
    current = None
    desc = ""
    chunks = []
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if current is not None:
                    records[current] = (desc, "".join(chunks))
                header = line[1:].strip()
                current = header.split()[0]
                desc = header[len(current):].strip()
                order.append(current)
                chunks = []
            else:
                chunks.append(line.strip())
        if current is not None:
            records[current] = (desc, "".join(chunks))
    return order, records


def write_fasta(path, order, records):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as out:
        for seq_id in order:
            out.write(">%s\n" % seq_id)
            seq = records[seq_id]
            for idx in range(0, len(seq), 80):
                out.write(seq[idx:idx + 80] + "\n")


def key_variants(value):
    variants = {value}
    variants.add(value.replace("|", "_"))
    variants.add(value.replace("|", "."))
    variants.add(re.sub(r"[^A-Za-z0-9_.:-]", "_", value))
    for suffix in (".p", ".protein", ".cds"):
        if value.endswith(suffix):
            variants.add(value[: -len(suffix)])
    if ".v" in value:
        variants.add(value.rsplit(".v", 1)[0])
    if "|" in value:
        variants.add(value.split("|", 1)[0])
        variants.add(value.split("|", 1)[1])
    return variants


def candidate_ids(seq_id, desc):
    ids = {seq_id}
    text = "%s %s" % (seq_id, desc)
    for pattern in [
        r"(?:transcript|transcript_id|protein_id|polypeptide|gene|gene_id|locus|locus_tag|Parent|ID)=([^;\s]+)",
        r"\[protein_id=([^\]]+)\]",
        r"\[gene=([^\]]+)\]",
        r"\[locus_tag=([^\]]+)\]",
    ]:
        for hit in re.findall(pattern, text):
            ids.add(hit.split()[0].strip())
    if "." in seq_id:
        ids.add(seq_id.rsplit(".", 1)[0])
    return ids


def parse_gene_cell(value):
    value = value.strip()
    if not value:
        return []
    return [part.strip() for part in value.replace(",", " ").split() if part.strip()]


def load_sequence_indexes(canonical_manifest, species):
    cds_by_species = {}
    pep_by_species = {}
    for row in read_tsv(canonical_manifest):
        label = row["label"]
        if label not in species:
            continue
        canonical_dir = Path(row["canonical_dir"])
        cds_path = canonical_dir / "primary.cds.fa"
        pep_path = canonical_dir / "primary.pep.fa"
        if not cds_path.is_file():
            raise SystemExit("ERROR: missing primary CDS for %s: %s" % (label, cds_path))
        if not pep_path.is_file():
            raise SystemExit("ERROR: missing primary protein for %s: %s" % (label, pep_path))
        _, cds_records = read_fasta(cds_path)
        _, pep_records = read_fasta(pep_path)
        cds_index = {}
        pep_index = {}
        for seq_id, (desc, seq) in cds_records.items():
            for candidate in candidate_ids(seq_id, desc):
                for key in key_variants(candidate):
                    cds_index.setdefault(key, seq.upper())
        for seq_id, (desc, seq) in pep_records.items():
            for candidate in candidate_ids(seq_id, desc):
                for key in key_variants(candidate):
                    pep_index.setdefault(key, seq.upper())
        for seq_id, (_, seq) in cds_records.items():
            for key in key_variants(seq_id):
                cds_index.setdefault(key, seq.upper())
        for seq_id, (_, seq) in pep_records.items():
            for key in key_variants(seq_id):
                pep_index.setdefault(key, seq.upper())
        cds_by_species[label] = cds_index
        pep_by_species[label] = pep_index
    missing = [sp for sp in species if sp not in cds_by_species or sp not in pep_by_species]
    if missing:
        raise SystemExit("ERROR: species absent from canonical manifest: %s" % ",".join(missing))
    return cds_by_species, pep_by_species


def find_sequence(index, gene_id):
    for key in key_variants(gene_id):
        if key in index:
            return index[key]
    for key in key_variants(gene_id.split("_", 1)[-1]):
        if key in index:
            return index[key]
    return None


def valid_cds(seq):
    clean = seq.upper().replace("-", "")
    return bool(clean) and len(clean) % 3 == 0 and re.fullmatch(r"[ACGTN]+", clean) is not None


GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate_cds(seq):
    clean = seq.upper().replace("-", "")
    amino = []
    for idx in range(0, len(clean), 3):
        codon = clean[idx:idx + 3]
        amino.append(GENETIC_CODE.get(codon, "X"))
    return "".join(amino)


def protein_cds_match(protein, cds):
    protein = protein.upper().replace("-", "").rstrip("*")
    translated = translate_cds(cds).rstrip("*")
    if "*" in protein or "*" in translated:
        return False
    if len(protein) != len(translated):
        return False
    for aa, bb in zip(protein, translated):
        if aa != bb and aa != "X" and bb != "X":
            return False
    return True


def cds_for_codon_alignment(protein, cds):
    """Return CDS sequence matched to the protein length for codon alignment."""
    protein = protein.upper().replace("-", "").rstrip("*")
    clean = cds.upper().replace("-", "")
    translated = translate_cds(clean)
    if translated.endswith("*") and len(translated) == len(protein) + 1:
        clean = clean[:-3]
    return clean


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--canonical-manifest", required=True)
    parser.add_argument("--orthogroups-tsv", required=True)
    parser.add_argument("--species", required=True, help="Comma-separated species labels in desired alignment order.")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--min-aa-aligned", type=int, default=60)
    args = parser.parse_args()

    species = [x.strip() for x in args.species.split(",") if x.strip()]
    out_dir = Path(args.out_dir)
    aa_dir = out_dir / "pal2nal_inputs" / "aa_unaligned"
    cds_dir = out_dir / "pal2nal_inputs" / "cds"
    stats_path = out_dir / "single_copy_extract_stats.tsv"
    counts_path = out_dir / "single_copy_extract_status_counts.tsv"
    single_copy_path = out_dir / "single_copy_orthogroups.txt"
    aa_dir.mkdir(parents=True, exist_ok=True)
    cds_dir.mkdir(parents=True, exist_ok=True)

    cds_indexes, pep_indexes = load_sequence_indexes(args.canonical_manifest, species)
    og_rows = read_tsv(args.orthogroups_tsv)

    used = 0
    skipped = 0
    stats = []
    selected = []

    for row in og_rows:
        og = row.get("Orthogroup") or row.get("HOG") or row.get("OG") or ""
        if not og:
            continue
        single = {}
        for sp in species:
            genes = parse_gene_cell(row.get(sp, ""))
            if len(genes) != 1:
                skipped += 1
                stats.append({"orthogroup": og, "status": "skip_not_single_copy_all_species"})
                single = {}
                break
            single[sp] = genes[0]
        if not single:
            continue

        aa_records = {}
        cds_records = {}
        try:
            for sp in species:
                gene_id = single[sp]
                prot_seq = find_sequence(pep_indexes[sp], gene_id)
                if not prot_seq:
                    raise ValueError("missing_protein_%s" % sp)
                if len(prot_seq.replace("*", "")) < args.min_aa_aligned:
                    raise ValueError("too_short_%s" % sp)
                cds_seq = find_sequence(cds_indexes[sp], gene_id)
                if not cds_seq:
                    raise ValueError("missing_cds_%s" % sp)
                if not valid_cds(cds_seq):
                    raise ValueError("invalid_cds_%s" % sp)
                if not protein_cds_match(prot_seq, cds_seq):
                    raise ValueError("translation_mismatch_%s" % sp)
                aa_records[sp] = prot_seq.upper().rstrip("*")
                cds_records[sp] = cds_for_codon_alignment(prot_seq, cds_seq)
        except ValueError as exc:
            skipped += 1
            stats.append({"orthogroup": og, "status": "skip_%s" % str(exc)})
            continue

        if set(aa_records) != set(species) or set(cds_records) != set(species):
            skipped += 1
            stats.append({"orthogroup": og, "status": "skip_taxon_mismatch"})
            continue

        write_fasta(aa_dir / ("%s.fa" % og), species, aa_records)
        write_fasta(cds_dir / ("%s.fa" % og), species, cds_records)
        selected.append(og)
        used += 1
        stats.append({"orthogroup": og, "status": "used"})

    with stats_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["orthogroup", "status"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(stats)
    status_counts = Counter(row["status"] for row in stats)
    with counts_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["status", "count"], lineterminator="\n")
        writer.writeheader()
        for status, count in status_counts.most_common():
            writer.writerow({"status": status, "count": count})
    with single_copy_path.open("w", encoding="utf-8") as handle:
        for og in selected:
            handle.write(og + "\n")

    print("used_single_copy=%d skipped=%d aa_dir=%s cds_dir=%s counts=%s" % (used, skipped, aa_dir, cds_dir, counts_path))


if __name__ == "__main__":
    main()
