#!/usr/bin/env python3
"""Build neutral JCVI inputs for copy1, copy2 and four PZI haplotypes.

The four haplotype files are treated as independent sequences.  Their source
names are retained only in ``source_manifest.tsv``; no A/B biological label is
written to BED/FASTA records.  Copy2 proteins and CDS are reconstructed from
the copy2 GFF3 and genome because no ready-made copy2 peptide/CDS package is
present in the project.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

from Bio.Seq import Seq
import pysam


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--copy1-dir", required=True)
    p.add_argument("--copy2-gff", required=True)
    p.add_argument("--copy2-genome", required=True)
    p.add_argument("--hap", action="append", nargs=2, metavar=("NEUTRAL_LABEL", "INPUT_DIR"), required=True,
                   help="Repeat four times; INPUT_DIR is an existing JCVI bed/pep/cds directory")
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def attrs(text):
    out = {}
    for item in text.strip().split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
    return out


def fasta(path):
    records = {}
    name = None
    seq = []
    with Path(path).open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(seq).upper()
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if name is not None:
        records[name] = "".join(seq).upper()
    return records


def natural(value):
    text = ""
    out = []
    for ch in value:
        if ch.isdigit():
            text += ch
        else:
            if text:
                out.append(int(text)); text = ""
            out.append(ch)
    if text:
        out.append(int(text))
    return out


def write_fasta(path, records):
    with Path(path).open("w", encoding="ascii", newline="\n") as fh:
        for name, seq in records:
            fh.write(">" + name + "\n")
            for start in range(0, len(seq), 80):
                fh.write(seq[start:start + 80] + "\n")


def read_bed_dir(input_dir, label, out_dir):
    """Neutralize an existing JCVI directory and return an audit row."""
    src = Path(input_dir)
    bed_path = src / next((x.name for x in src.glob("*.bed")), "")
    pep_path = src / next((x.name for x in src.glob("*.pep*")), "")
    cds_path = src / next((x.name for x in src.glob("*.cds*")), "")
    seqids_path = src / next((x.name for x in src.glob("*.seqids")), "")
    for path in (bed_path, pep_path, cds_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Missing JCVI input in {src}: {path}")
    pep = fasta(pep_path)
    cds = fasta(cds_path)
    seqids = []
    if seqids_path.is_file():
        seqids = [x for x in seqids_path.read_text(encoding="utf-8", errors="replace").strip().split(",") if x]
    rows = []
    with bed_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            rows.append(f[:6])
    observed = []
    for row in rows:
        if row[0] not in observed:
            observed.append(row[0])
    order = [x for x in seqids if x in observed] + [x for x in observed if x not in seqids]
    chr_map = {old: f"{label}_chr{i:02d}" for i, old in enumerate(order, 1)}
    out_dir.mkdir(parents=True, exist_ok=True)
    neutral = []
    pep_out = []
    cds_out = []
    retained = 0
    for index, row in enumerate(rows, 1):
        old = row[3]
        if old not in pep or old not in cds or not pep[old] or not cds[old]:
            continue
        new = f"{label}_g{index:06d}"
        row = [chr_map.get(row[0], row[0]), row[1], row[2], new, row[4], row[5]]
        neutral.append(row)
        pep_out.append((new, pep[old]))
        cds_out.append((new, cds[old]))
        retained += 1
    with (out_dir / f"{label}.bed").open("w", encoding="ascii", newline="\n") as fh:
        for row in neutral:
            fh.write("\t".join(row) + "\n")
    write_fasta(out_dir / f"{label}.pep", pep_out)
    write_fasta(out_dir / f"{label}.cds", cds_out)
    (out_dir / f"{label}.seqids").write_text(",".join(chr_map[x] for x in order) + "\n", encoding="ascii")
    return {"label": label, "source": str(src), "input_models": len(rows), "retained_models": retained,
            "chromosomes": len(order), "protein_records": len(pep_out), "cds_records": len(cds_out)}


def extract_copy2(gff_path, genome_path, out_dir):
    """Extract one longest coding transcript per gene from copy2."""
    gff_path = Path(gff_path)
    out_dir.mkdir(parents=True, exist_ok=True)
    genome = pysam.FastaFile(str(genome_path))
    mrna = {}
    cds_groups = defaultdict(list)
    gene_parent = {}
    seq_order = []
    with gff_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                continue
            seqid, _src, feature, start, end, _score, strand, phase, text = f
            try:
                start_i, end_i = int(start), int(end)
            except ValueError:
                continue
            if seqid not in seq_order:
                seq_order.append(seqid)
            at = attrs(text)
            if feature.lower() in {"mrna", "transcript"}:
                mid = at.get("ID")
                if not mid:
                    continue
                parent = at.get("Parent", mid)
                mrna[mid] = {"seqid": seqid, "start": min(start_i, end_i), "end": max(start_i, end_i),
                             "strand": strand, "gene": parent}
                gene_parent[mid] = parent
            elif feature.lower() == "cds":
                parents = [x for x in at.get("Parent", "").split(",") if x]
                for parent in parents:
                    cds_groups[parent].append((start_i, end_i, strand, phase))
    candidates = []
    for mid, rec in mrna.items():
        segments = cds_groups.get(mid, [])
        total = sum(abs(x[1] - x[0]) + 1 for x in segments)
        if total < 90:
            continue
        candidates.append((rec["gene"], mid, rec, segments, total))
    by_gene = defaultdict(list)
    for row in candidates:
        by_gene[row[0]].append(row)
    chosen = []
    for gene, values in by_gene.items():
        chosen.append(max(values, key=lambda x: (x[4], -x[2]["start"])))
    chosen.sort(key=lambda x: (natural(x[2]["seqid"]), x[2]["start"], x[2]["end"], x[1]))
    chr_map = {old: f"copy2_chr{i:02d}" for i, old in enumerate(seq_order, 1)}
    bed = []
    pep_out = []
    cds_out = []
    skipped = defaultdict(int)
    for index, (_gene, mid, rec, segments, _total) in enumerate(chosen, 1):
        strand = rec["strand"]
        segments = sorted(segments, key=lambda x: (x[0], x[1]), reverse=(strand == "-"))
        parts = []
        try:
            for start_i, end_i, seg_strand, _phase in segments:
                seq = genome.fetch(rec["seqid"], min(start_i, end_i) - 1, max(start_i, end_i)).upper()
                if strand == "-":
                    seq = str(Seq(seq).reverse_complement())
                parts.append(seq)
            coding = "".join(parts)
            if len(coding) < 90 or len(coding) % 3:
                skipped["nontriplet_or_short"] += 1
                continue
            protein = str(Seq(coding).translate(to_stop=False))
            if protein.endswith("*"):
                protein = protein[:-1]
            if "*" in protein or not protein:
                skipped["internal_stop"] += 1
                continue
        except Exception:
            skipped["sequence_fetch"] += 1
            continue
        new = f"copy2_g{index:06d}"
        bed.append([chr_map.get(rec["seqid"], rec["seqid"]), str(rec["start"] - 1), str(rec["end"]), new, "0", strand])
        pep_out.append((new, protein))
        cds_out.append((new, coding))
    with (out_dir / "copy2.bed").open("w", encoding="ascii", newline="\n") as fh:
        for row in bed:
            fh.write("\t".join(row) + "\n")
    write_fasta(out_dir / "copy2.pep", pep_out)
    write_fasta(out_dir / "copy2.cds", cds_out)
    (out_dir / "copy2.seqids").write_text(",".join(chr_map[x] for x in seq_order) + "\n", encoding="ascii")
    rows = [{"label": "copy2", "source": str(gff_path), "input_models": len(chosen),
             "retained_models": len(bed), "chromosomes": len(seq_order), "protein_records": len(pep_out),
             "cds_records": len(cds_out), "skipped": ";".join(f"{k}={v}" for k, v in sorted(skipped.items()))}]
    return rows[0]


def main() -> int:
    a = parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    rows.append(read_bed_dir(a.copy1_dir, "copy1", out / "copy1"))
    rows.append(extract_copy2(a.copy2_gff, a.copy2_genome, out / "copy2"))
    for label, src in a.hap:
        rows.append(read_bed_dir(src, label, out / label))
    with (out / "source_manifest.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted(set().union(*(r.keys() for r in rows))), delimiter="\t", lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow(row)
    (out / "README.txt").write_text(
        "Neutral PZI internal synteny inputs. copy1/copy2/hap01-hap04 are analysis labels only.\n"
        "The hap01-hap04 source names are recorded in source_manifest.tsv and are not biological\n"
        "A/B assignments. copy2 proteins/CDS were extracted from the Helixer GFF3 and copy2 genome\n"
        "using the longest coding transcript per gene.\n", encoding="ascii")
    print(f"[OK] Neutral internal inputs prepared: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
