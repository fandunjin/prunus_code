#!/usr/bin/env python3
"""Extract real TRF monomers and build true multiple-sequence-logo inputs.

The input ``tandem_top_arrays.tsv`` contains array coordinates and TRF periods,
not independent repeat-unit sequences.  This script fetches each array from
the unmasked genome, cuts the sequence into period-sized monomers, orients and
phase-normalizes them against the TRF consensus seed, then aligns the resulting
monomers with MAFFT.  The resulting alignments are the only inputs used for a
true information-content logo. ``--align-all`` retains every QC-passing
monomer; capped balanced sampling remains available for very large arrays.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from collections import defaultdict
from pathlib import Path


PZI_TARGETS = ("pziA1", "pziA2", "pziB1", "pziB2")


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def canonical(value):
    value = clean_sequence(value)
    if not value:
        return ""
    comp = str.maketrans("ACGT", "TGCA")
    rc = value.translate(comp)[::-1]
    return min(value, rc)


def clean_sequence(value):
    """Return only DNA bases, dropping labels such as ``long_500bp:``."""
    value = (value or "").split(":", 1)[-1].upper()
    return "".join(base for base in value if base in "ACGT")


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def fai_records(path):
    records = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 5:
                records[fields[0]] = tuple(map(int, fields[1:5]))
    if not records:
        raise ValueError(f"invalid FASTA index: {path}")
    return records


def fetch_fasta(fasta, index, seqid, start, end):
    if seqid not in index:
        raise KeyError(f"{seqid} absent from {fasta}.fai")
    length, byte_offset, line_bases, line_width = index[seqid]
    start, end = max(0, start), min(length, end)
    if end <= start:
        return ""
    pieces = []
    with fasta.open("rb") as handle:
        position = start
        while position < end:
            within = position % line_bases
            take = min(end - position, line_bases - within)
            offset = byte_offset + (position // line_bases) * line_width + within
            handle.seek(offset)
            pieces.append(handle.read(take).decode("ascii"))
            position += take
    return "".join(pieces).upper()


def resolve_genome_path(root, target):
    """Resolve either standardized target names or the original PZI names."""
    candidates = [
        root / target / f"{target}.fa",
        root / f"{target}.fa",
    ]
    if target.startswith("pzi") and len(target) == 5:
        candidates.append(root / f"pzi_hap{target[-2:]}.fa")
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return candidates[0]


def identity(a, b):
    n = min(len(a), len(b))
    if not n:
        return 0.0
    return sum(x == y for x, y in zip(a[:n], b[:n])) / n


def cyclic_variants(sequence):
    """Yield every strand/phase representation of a repeat seed."""
    sequence = clean_sequence(sequence)
    if not sequence:
        return []
    reverse = revcomp(sequence)
    return [sequence[i:] + sequence[:i] for i in range(len(sequence))] + \
           [reverse[i:] + reverse[:i] for i in range(len(reverse))]


def seed_identity(seed, consensus):
    """Best identity of a short seed against a full TRF consensus.

    TRF may report the same repeat from a different phase or strand.  Matching
    every cyclic seed and every consensus window avoids treating those labels
    as different repeat families.
    """
    seed = clean_sequence(seed)
    consensus = clean_sequence(consensus)
    if not seed or not consensus:
        return 0.0
    best = 0.0
    for variant in cyclic_variants(seed):
        if len(consensus) < len(variant):
            best = max(best, identity(variant, consensus))
            continue
        for start in range(len(consensus) - len(variant) + 1):
            best = max(best, identity(variant, consensus[start:start + len(variant)]))
    return best


def phase_orient(monomer, anchor):
    """Orient and circularly phase a monomer using a short consensus seed."""
    if not monomer:
        return ""
    anchor = clean_sequence(anchor)[: min(30, len(anchor), len(monomer))]
    if not anchor:
        return monomer
    candidates = (monomer, revcomp(monomer))
    best = (float("-inf"), monomer)
    for candidate in candidates:
        doubled = candidate + candidate
        for offset in range(len(candidate)):
            score = identity(doubled[offset: offset + len(anchor)], anchor)
            if score > best[0]:
                best = (score, doubled[offset: offset + len(candidate)])
    return best[1]


def sample_indices(n, max_n):
    if n <= max_n:
        return list(range(n))
    # Deterministic evenly spaced sampling preserves within-array variation
    # while preventing a few multi-megabase arrays from dominating the logo.
    return sorted({int(round(i * (n - 1) / (max_n - 1))) for i in range(max_n)})


def run_mafft(mafft, input_path, output_path, threads):
    command = [mafft, "--auto", "--thread", str(threads), "--quiet", str(input_path)]
    with output_path.open("w", encoding="utf-8") as handle:
        result = subprocess.run(command, stdout=handle, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"MAFFT failed for {input_path.name}: {result.stderr[-1000:]}")


def alignment_records(path):
    records = []
    name, sequence = None, []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(sequence).upper()))
                name, sequence = line[1:].split()[0], []
            else:
                sequence.append(line)
    if name is not None:
        records.append((name, "".join(sequence).upper()))
    return records


def alignment_column_rows(cluster_id, records):
    rows = []
    if not records:
        return rows
    length = max(len(sequence) for _, sequence in records)
    for position in range(length):
        column = [sequence[position] if position < len(sequence) else "-"
                  for _, sequence in records]
        counts = {letter: column.count(letter) for letter in "ACGT"}
        valid = sum(counts.values())
        gaps = column.count("-") + column.count("N")
        if valid:
            frequencies = [counts[letter] / valid for letter in "ACGT"]
            entropy = -sum(freq * math.log(freq, 2) for freq in frequencies if freq > 0)
            correction = 3.0 / (2.0 * math.log(2.0) * valid)
            information = max(0.0, 2.0 - entropy - correction)
        else:
            information = 0.0
        rows.append({
            "cluster_id": cluster_id, "alignment_position": position + 1,
            "n_sequences": len(column), "n_valid_bases": valid, "n_gaps_or_N": gaps,
            "A_count": counts["A"], "C_count": counts["C"],
            "G_count": counts["G"], "T_count": counts["T"],
            "A_frequency": f"{counts['A'] / valid:.6f}" if valid else "0",
            "C_frequency": f"{counts['C'] / valid:.6f}" if valid else "0",
            "G_frequency": f"{counts['G'] / valid:.6f}" if valid else "0",
            "T_frequency": f"{counts['T'] / valid:.6f}" if valid else "0",
            "information_bits": f"{information:.6f}",
        })
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arrays", type=Path, required=True)
    parser.add_argument("--motif-clusters", type=Path, required=True)
    parser.add_argument("--motif-members", type=Path, required=True)
    parser.add_argument("--genome-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--mafft", default="mafft")
    parser.add_argument("--targets", nargs="+", default=list(PZI_TARGETS))
    parser.add_argument("--clusters", type=int, default=4)
    parser.add_argument("--max-monomers-per-array", type=int, default=40)
    parser.add_argument("--max-monomers-per-cluster", type=int, default=400)
    parser.add_argument("--align-all", action="store_true",
                        help="Align every QC-passing extracted monomer; otherwise use balanced caps")
    parser.add_argument("--min-array-copies", type=int, default=3)
    parser.add_argument("--min-seed-identity", type=float, default=0.80,
                        help="Minimum identity for assigning a long-array consensus to a motif cluster")
    parser.add_argument("--min-anchor-identity", type=float, default=0.60,
                        help="Minimum identity after monomer strand/phase normalization")
    parser.add_argument("--coordinate-system", choices=("zero-based-half-open", "one-based-inclusive"),
                        default="zero-based-half-open",
                        help="Coordinate convention used by tandem_top_arrays.tsv")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--seed", type=int, default=20260821)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    monomer_dir = args.output_dir / "monomers"
    alignment_dir = args.output_dir / "alignments"
    monomer_dir.mkdir(exist_ok=True)
    alignment_dir.mkdir(exist_ok=True)

    cluster_rows = read_tsv(args.motif_clusters)
    member_rows = read_tsv(args.motif_members)
    members_by_cluster = defaultdict(list)
    for row in member_rows:
        members_by_cluster[row["cluster_id"]].append(row["motif"])
    cluster_period = {row["cluster_id"]: int(float(row["period_bp"])) for row in cluster_rows}
    cluster_seeds = {
        row["cluster_id"]: clean_sequence(row.get("representative_motif", ""))
        for row in cluster_rows
    }

    arrays = []
    for row in read_tsv(args.arrays):
        target = row.get("target", "")
        if target not in args.targets:
            continue
        try:
            period = int(float(row["period_bp"]))
            start, end = int(row["start"]), int(row["end"])
            if args.coordinate_system == "one-based-inclusive":
                start -= 1
            copies = int(float(row.get("copy_number", 0)))
        except (KeyError, ValueError):
            continue
        if period < 2 or end <= start or copies < args.min_array_copies:
            continue
        raw_consensus = clean_sequence(row.get("raw_consensus_motif", ""))
        motif = canonical(row.get("canonical_motif") or raw_consensus)
        arrays.append({**row, "target": target, "period": period, "start_i": start, "end_i": end,
                       "copy_number_i": copies, "motif_canonical": motif,
                       "raw_consensus_clean": raw_consensus})

    # Match array seeds to motif clusters.  A cluster is retained only when it
    # has enough coordinate-bearing arrays to estimate a logo.
    cluster_arrays = defaultdict(list)
    assignment_rows = []
    for row in arrays:
        candidates = []
        for cluster_id, motifs in members_by_cluster.items():
            cluster_period_i = cluster_period.get(cluster_id)
            if cluster_period_i is None or abs(row["period"] - cluster_period_i) > max(2, int(0.05 * cluster_period_i)):
                continue
            seed = cluster_seeds.get(cluster_id, "")
            scores = [seed_identity(seed, row["raw_consensus_clean"]),
                      seed_identity(seed, row["motif_canonical"])]
            scores.extend(seed_identity(motif, row["raw_consensus_clean"]) for motif in motifs)
            candidates.append((max(scores), cluster_id))
        if candidates:
            best_score, best_cluster = max(candidates)
            if best_score >= args.min_seed_identity:
                row["cluster_match_identity"] = best_score
                row["assigned_cluster"] = best_cluster
                cluster_arrays[best_cluster].append(row)
                assignment_rows.append({
                    "target": row["target"], "seqid": row["seqid"],
                    "start": row["start_i"], "end": row["end_i"],
                    "period_bp": row["period"], "array_length_bp": row["end_i"] - row["start_i"],
                    "assigned_cluster": best_cluster,
                    "seed_identity": f"{best_score:.4f}",
                    "canonical_motif": row.get("canonical_motif", ""),
                    "raw_consensus_length_bp": len(row.get("raw_consensus_clean", "")),
                })

    ranked = sorted(cluster_arrays, key=lambda cluster: sum(int(r["end_i"] - r["start_i"]) for r in cluster_arrays[cluster]), reverse=True)
    selected = [cluster for cluster in ranked if len(cluster_arrays[cluster]) >= args.min_array_copies][: args.clusters]
    if not selected:
        raise SystemExit("ERROR: no motif clusters have enough coordinate-bearing arrays for true logos")

    summary = []
    column_rows = []
    for cluster_id in selected:
        rows = cluster_arrays[cluster_id]
        period = cluster_period[cluster_id]
        anchor = clean_sequence(cluster_seeds.get(cluster_id) or members_by_cluster[cluster_id][0])
        raw_records = []
        target_counts = defaultdict(int)
        anchor_identities = []
        fasta_cache = {}
        for row in rows:
            fasta = resolve_genome_path(args.genome_root, row["target"])
            fai = fasta.with_suffix(fasta.suffix + ".fai")
            if not fasta.is_file() or not fai.is_file():
                raise FileNotFoundError(f"missing genome or FASTA index: {fasta} / {fai}")
            if row["target"] not in fasta_cache:
                fasta_cache[row["target"]] = (fasta, fai_records(fai))
            genome, index = fasta_cache[row["target"]]
            sequence = fetch_fasta(genome, index, row["seqid"], row["start_i"], row["end_i"])
            n_monomers = len(sequence) // period
            if n_monomers < 3:
                continue
            indices = (list(range(n_monomers)) if args.align_all else
                       sample_indices(n_monomers, args.max_monomers_per_array))
            for monomer_index in indices:
                monomer = sequence[monomer_index * period: (monomer_index + 1) * period]
                if monomer.count("N") / max(1, len(monomer)) > 0.05:
                    continue
                phased = phase_orient(monomer, anchor)
                phased_identity = identity(phased[:len(anchor)], anchor)
                if phased_identity < args.min_anchor_identity:
                    # Anchor identity is only a conservative orientation check;
                    # divergent monomers are retained when their phase is clear.
                    continue
                raw_records.append((f"{row['target']}|{row['seqid']}:{row['start_i']}-{row['end_i']}|m{monomer_index + 1}", phased))
                target_counts[row["target"]] += 1
                anchor_identities.append(phased_identity)
        if len(raw_records) < args.min_array_copies:
            continue
        sampled = raw_records[:]
        if not args.align_all and len(sampled) > args.max_monomers_per_cluster:
            sampled = [sampled[i] for i in sample_indices(len(sampled), args.max_monomers_per_cluster)]
        raw_path = monomer_dir / f"{cluster_id}.monomers.fa"
        with raw_path.open("w", encoding="utf-8") as handle:
            for name, sequence in sampled:
                handle.write(f">{name}\n{sequence}\n")
        aln_path = alignment_dir / f"{cluster_id}.aligned.fa"
        run_mafft(args.mafft, raw_path, aln_path, args.threads)
        alignment = alignment_records(aln_path)
        column_rows.extend(alignment_column_rows(cluster_id, alignment))
        summary.append({
            "cluster_id": cluster_id,
            "period_bp": period,
            "motif_member_count": len(members_by_cluster[cluster_id]),
            "source_array_count": len(rows),
            "source_arrays": ";".join(
                f"{row['target']}:{row['seqid']}:{row['start_i']}-{row['end_i']}"
                for row in rows
            ),
            "raw_monomer_count": len(raw_records),
            "aligned_monomer_count": len(sampled),
            "mean_anchor_identity": f"{sum(anchor_identities) / len(anchor_identities):.4f}" if anchor_identities else "",
            "min_array_seed_identity": f"{min(float(row.get('cluster_match_identity', 0.0)) for row in rows):.4f}" if rows else "",
            "target_counts": ",".join(f"{k}:{target_counts[k]}" for k in sorted(target_counts)),
            "monomer_fasta": raw_path.name,
            "alignment_fasta": aln_path.name,
            "alignment_method": "MAFFT --auto after TRF-period extraction and phase/orientation normalization",
        })
    if not summary:
        raise SystemExit("ERROR: no cluster produced enough valid monomers")
    write_tsv(args.output_dir / "trf_true_logo_summary.tsv", summary)
    write_tsv(args.output_dir / "trf_true_logo_array_assignments.tsv", assignment_rows)
    write_tsv(args.output_dir / "trf_true_logo_columns.tsv", column_rows)
    manifest = {
        "definition": "monomers are extracted from genome coordinates and cut by the TRF period; they are not consensus copies",
        "phase_normalization": "best circular phase and strand against the canonical TRF seed",
        "alignment": "MAFFT --auto",
        "genome_source": str(args.genome_root),
        "coordinate_system": args.coordinate_system,
        "min_seed_identity": args.min_seed_identity,
        "max_monomers_per_array": args.max_monomers_per_array,
        "max_monomers_per_cluster": args.max_monomers_per_cluster,
        "align_all": args.align_all,
        "note": "A compact logo should display a diagnostic window of the full period; full-length alignments are retained in alignments/.",
    }
    (args.output_dir / "trf_true_logo_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({"clusters": len(summary), "output_dir": str(args.output_dir)}, indent=2))


if __name__ == "__main__":
    main()
