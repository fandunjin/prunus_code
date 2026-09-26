#!/usr/bin/env python3
"""Compare neutral copy1/GTR chromosome models from screened JCVI anchors.

The comparison is deliberately label-free.  It uses copy1 ``chr01``--``chr16``
and GTR ``Chr01``--``Chr09`` only.  The reported score is a transparent
comparative metric, not a formal phylogenetic likelihood:

  score = 100 * (0.6 * anchor_retention + 0.2 * order_continuity
                 + 0.2 * structural_coherence)
          - 2 * event_burden - 0.5 * extra_assignments

All component metrics and the exact assignments are written so that model
choice is not based on the single composite score alone.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path


Q_CHRS = [f"chr{i:02d}" for i in range(1, 17)]
G_CHRS = [f"Chr{i:02d}" for i in range(1, 10)]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--anchors", required=True)
    p.add_argument("--qbed", required=True)
    p.add_argument("--gtr-bed", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--min-cell-anchors", type=int, default=8)
    p.add_argument("--min-cell-fraction", type=float, default=0.05)
    p.add_argument("--min-segment-anchors", type=int, default=8)
    return p.parse_args()


def natural_key(value: str):
    out = []
    token = ""
    for ch in value:
        if ch.isdigit():
            token += ch
        else:
            if token:
                out.append(int(token)); token = ""
            out.append(ch)
    if token:
        out.append(int(token))
    return out


def read_bed(path: Path):
    rec = {}
    order_by_chr = defaultdict(list)
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError:
                continue
            chrom, gene = f[0], f[3]
            rec[gene] = {"chrom": chrom, "start": start, "end": end}
            order_by_chr[chrom].append(gene)
    ranks = {}
    for chrom, genes in order_by_chr.items():
        for i, gene in enumerate(genes):
            ranks[gene] = i
    return rec, ranks, order_by_chr


def parse_anchors(path: Path):
    """Return unique gene pairs from JCVI anchors.new, ignoring block headers."""
    pairs = []
    seen = set()
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            f = text.split()
            if len(f) < 2:
                continue
            pair = (f[0], f[1])
            if pair in seen:
                continue
            seen.add(pair)
            pairs.append(pair)
    return pairs


def write_tsv(path: Path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def cosine(a, b):
    aa = sum(x * x for x in a)
    bb = sum(x * x for x in b)
    if aa <= 0 or bb <= 0:
        return 0.0
    return sum(x * y for x, y in zip(a, b)) / math.sqrt(aa * bb)


def supported_labels(row, min_anchors, min_fraction):
    total = sum(row)
    if total <= 0:
        return [], 0.0
    threshold = max(min_anchors, min_fraction * total)
    labels = [G_CHRS[i] for i, n in enumerate(row) if n >= threshold]
    mass = sum(row[i] for i, label in enumerate(G_CHRS) if label in labels) / total
    return labels, mass


def _average_ranks(values):
    """Return average ranks for ties without requiring scipy."""
    ranked = sorted(enumerate(values), key=lambda item: item[1])
    out = [0.0] * len(values)
    i = 0
    while i < len(ranked):
        j = i + 1
        while j < len(ranked) and ranked[j][1] == ranked[i][1]:
            j += 1
        rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            out[ranked[k][0]] = rank
        i = j
    return out


def _pearson(values_a, values_b):
    if len(values_a) < 2:
        return 0.0
    mean_a = sum(values_a) / len(values_a)
    mean_b = sum(values_b) / len(values_b)
    numerator = sum((a - mean_a) * (b - mean_b) for a, b in zip(values_a, values_b))
    denom_a = math.sqrt(sum((a - mean_a) ** 2 for a in values_a))
    denom_b = math.sqrt(sum((b - mean_b) ** 2 for b in values_b))
    return numerator / (denom_a * denom_b) if denom_a and denom_b else 0.0


def monotonicity(records):
    """Within-cell syntenic order score in either orientation.

    A previous implementation sorted reference ranks *within each query rank*,
    which made every non-empty comparison score 1.0.  Here, duplicate anchors
    are collapsed to the median reference rank for each query rank, points are
    ordered along the query chromosome, and the score combines absolute
    Spearman correlation with the dominant direction of consecutive steps.
    """
    by_q = defaultdict(list)
    for q_rank, g_rank in records:
        by_q[q_rank].append(g_rank)
    points = [(q_rank, sorted(values)[len(values) // 2]) for q_rank, values in by_q.items()]
    points.sort()
    if len(points) < 3:
        return 0.0
    query_ranks = [q for q, _ in points]
    target_ranks = [g for _, g in points]
    rho = abs(_pearson(_average_ranks(query_ranks), _average_ranks(target_ranks)))
    steps = [b - a for a, b in zip(target_ranks, target_ranks[1:]) if b != a]
    if not steps:
        return 0.0
    positive = sum(step > 0 for step in steps)
    direction = max(positive, len(steps) - positive) / len(steps)
    return (rho + direction) / 2.0


def weighted_monotonicity(records_by_target):
    """Anchor-weighted within-reference-block order continuity."""
    total = sum(len(records) for records in records_by_target.values())
    if not total:
        return 0.0
    return sum(monotonicity(records) * len(records) for records in records_by_target.values()) / total


def perfect_pairing(weights):
    """Maximum-weight perfect pairing for up to 16 items by bitmask DP."""
    n = len(weights)

    @lru_cache(maxsize=None)
    def solve(mask):
        if mask == 0:
            return 0.0, ()
        i = (mask & -mask).bit_length() - 1
        rest = mask ^ (1 << i)
        best_score, best_pairs = -1.0, ()
        jmask = rest
        while jmask:
            bit = jmask & -jmask
            j = bit.bit_length() - 1
            score, pairs = solve(rest ^ bit)
            score += weights[i][j]
            if score > best_score:
                best_score, best_pairs = score, ((i, j),) + pairs
            jmask ^= bit
        return best_score, best_pairs

    return solve((1 << n) - 1)


def assignment_score(rows, row_indices, one_to_one=True):
    """Best normalized row-to-GTR assignment score and mapping."""
    if not row_indices:
        return 0.0, {}
    weights = []
    for i in row_indices:
        total = sum(rows[i]) or 1
        weights.append([rows[i][j] / total for j in range(len(G_CHRS))])
    if one_to_one:
        # Small dynamic program for rectangular assignment.  Each row receives
        # a distinct GTR block; rows cannot be assigned to unplaced sequences.
        @lru_cache(maxsize=None)
        def dp(i, mask):
            if i == len(row_indices):
                return 0.0, ()
            best = (-1.0, ())
            for j in range(len(G_CHRS)):
                if not (mask & (1 << j)):
                    v, path = dp(i + 1, mask | (1 << j))
                    v += weights[i][j]
                    if v > best[0]:
                        best = (v, (j,) + path)
            return best
        total, path = dp(0, 0)
        return total / len(row_indices), {row_indices[i]: path[i] for i in range(len(row_indices))}
    vals = [max(row) for row in weights]
    return sum(vals) / len(vals), {row_indices[i]: max(range(len(G_CHRS)), key=lambda j: weights[i][j]) for i in range(len(row_indices))}


def main() -> int:
    a = parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    qrec, qrank, qgenes = read_bed(Path(a.qbed))
    grec, grank, ggenes = read_bed(Path(a.gtr_bed))
    q_chrs = [x for x in Q_CHRS if x in qgenes]
    g_chrs = [x for x in G_CHRS if x in ggenes]
    if q_chrs != Q_CHRS:
        raise SystemExit(f"Expected neutral copy1 chr01--chr16; found {q_chrs}")
    if g_chrs != G_CHRS:
        raise SystemExit(f"Expected GTR Chr01--Chr09; found {g_chrs}")

    matrix = {q: [0] * len(G_CHRS) for q in Q_CHRS}
    q_anchor_records = defaultdict(lambda: defaultdict(list))
    parsed = []
    for qgene, gg in parse_anchors(Path(a.anchors)):
        if qgene not in qrec or gg not in grec:
            continue
        q = qrec[qgene]["chrom"]
        g = grec[gg]["chrom"]
        if q not in Q_CHRS or g not in G_CHRS:
            continue
        qi, gi = qrank[qgene], grank[gg]
        matrix[q][G_CHRS.index(g)] += 1
        q_anchor_records[q][g].append((qi, gi))
        parsed.append({"q_gene": qgene, "gtr_gene": gg, "q_chr": q, "gtr_chr": g,
                       "q_rank": qi, "gtr_rank": gi})
    if not parsed:
        raise SystemExit("No usable copy1/GTR anchors after BED filtering")
    write_tsv(out / "anchor_pairs.tsv",
              ["q_gene", "gtr_gene", "q_chr", "gtr_chr", "q_rank", "gtr_rank"], parsed)
    write_tsv(out / "anchor_matrix.tsv", ["q_chr"] + G_CHRS,
              [{"q_chr": q, **{g: matrix[q][i] for i, g in enumerate(G_CHRS)}} for q in Q_CHRS])

    row_sets, row_mass = {}, {}
    row_cont = {}
    for q in Q_CHRS:
        row_sets[q], row_mass[q] = supported_labels(matrix[q], a.min_cell_anchors, a.min_cell_fraction)
        row_cont[q] = weighted_monotonicity(q_anchor_records[q])
    total_anchors = len(parsed)

    # K=8: profile-based neutral homolog pairing.  Pair scores use the full
    # nine-block count vector and do not assume any chromosomal names.
    weights = [[0.0] * 16 for _ in range(16)]
    norm = {}
    for i, q in enumerate(Q_CHRS):
        total = sum(matrix[q]) or 1
        norm[q] = [x / total for x in matrix[q]]
    for i in range(16):
        for j in range(i + 1, 16):
            si = set(row_sets[Q_CHRS[i]])
            sj = set(row_sets[Q_CHRS[j]])
            union = si | sj
            inter = si & sj
            jaccard = len(inter) / len(union) if union else 0.0
            weights[i][j] = weights[j][i] = 0.7 * cosine(norm[Q_CHRS[i]], norm[Q_CHRS[j]]) + 0.3 * jaccard
    pair_total, pairs = perfect_pairing(weights)
    k8_rows = []
    k8_extra = 0
    k8_events = 0
    k8_retained = 0
    for idx, (i, j) in enumerate(pairs, 1):
        q1, q2 = Q_CHRS[i], Q_CHRS[j]
        labels = sorted(set(row_sets[q1]) | set(row_sets[q2]), key=natural_key)
        retained = sum(matrix[q1][G_CHRS.index(g)] for g in labels) + sum(matrix[q2][G_CHRS.index(g)] for g in labels)
        k8_retained += retained
        extra = max(0, len(labels) - 1)
        k8_extra += extra
        k8_events += int(len(labels) > 1)
        k8_rows.append({"pair_id": f"K8_pair{idx:02d}", "copy1_chr_1": q1, "copy1_chr_2": q2,
                        "supported_gtr_blocks": ",".join(labels) if labels else "NA",
                        "pair_profile_similarity": f"{weights[i][j]:.6f}", "anchor_retained": retained,
                        "fused_block_flag": int(len(labels) > 1), "order_continuity_mean": f"{(row_cont[q1]+row_cont[q2])/2:.6f}"})
    k8_retention = k8_retained / total_anchors
    k8_order = sum(row_cont.values()) / len(row_cont)
    k8_coherence = pair_total / len(pairs)
    pair_of_q = {}
    for idx, (i, j) in enumerate(pairs, 1):
        pair_of_q[Q_CHRS[i]] = f"K8_pair{idx:02d}"
        pair_of_q[Q_CHRS[j]] = f"K8_pair{idx:02d}"

    # K=9-block: keep the nine GTR chromosomes as ancestral blocks and allow
    # each current copy1 chromosome to carry multiple supported blocks.
    k9_rows = []
    k9_extra = 0
    k9_events = 0
    k9_retained = 0
    for q in Q_CHRS:
        labels = row_sets[q]
        retained = sum(matrix[q][G_CHRS.index(g)] for g in labels)
        k9_retained += retained
        extra = max(0, len(labels) - 1)
        k9_extra += extra
        k9_events += int(len(labels) > 1)
        k9_rows.append({"copy1_chr": q, "supported_gtr_blocks": ",".join(labels) if labels else "NA",
                        "anchor_total": sum(matrix[q]), "anchor_retained": retained,
                        "retention": f"{retained/(sum(matrix[q]) or 1):.6f}",
                        "fused_block_flag": int(len(labels) > 1),
                        "order_continuity": f"{row_cont[q]:.6f}"})
    k9_retention = k9_retained / total_anchors
    k9_order = sum(row_cont.values()) / len(row_cont)
    k9_coherence = sum(row_mass.values()) / len(row_mass)

    # K=9+7: enumerate all 9/7 partitions.  The nine-member group is scored
    # against a one-to-one assignment to GTR blocks; the seven-member group
    # is allowed a block union (the expected fusion-tolerant side).
    best = None
    ranked_partitions = []
    row_vectors = [matrix[q] for q in Q_CHRS]
    q_indices = range(16)
    from itertools import combinations
    for nine_tuple in combinations(q_indices, 9):
        nine = set(nine_tuple)
        seven = [i for i in q_indices if i not in nine]
        fit9, map9 = assignment_score(row_vectors, list(nine_tuple), one_to_one=True)
        fit7, map7 = assignment_score(row_vectors, seven, one_to_one=False)
        # Assignment score is normalized per chromosome.  The reported model
        # score also includes the same order and event terms used above.
        labels9 = {i: [G_CHRS[map9[i]]] for i in map9}
        labels7 = {i: row_sets[Q_CHRS[i]] for i in seven}
        all_labels = {**labels9, **labels7}
        retained = sum(sum(matrix[Q_CHRS[i]][G_CHRS.index(g)] for g in labels) for i, labels in all_labels.items())
        retention = retained / total_anchors
        fused = sum(1 for labels in all_labels.values() if len(labels) > 1)
        extra = sum(max(0, len(labels) - 1) for labels in all_labels.values())
        order = sum(row_cont.values()) / len(row_cont)
        coherence = (fit9 * 9.0 + fit7 * 7.0) / 16.0
        score = 100 * (0.6 * retention + 0.2 * order + 0.2 * coherence) - 2 * fused - 0.5 * extra
        candidate = (score, nine_tuple, fit9, fit7, retention, fused, extra, coherence, all_labels, map9)
        ranked_partitions.append((score, nine_tuple, fit9, fit7, retention, fused, extra, coherence))
        if best is None or candidate[0] > best[0]:
            best = candidate
    assert best is not None
    k97_score, nine_tuple, fit9, fit7, k97_retention, k97_events, k97_extra, k97_coherence, labels97, map9 = best
    nine = set(nine_tuple)
    k97_rows = []
    for i, q in enumerate(Q_CHRS):
        group = "9-chromosome set" if i in nine else "7-chromosome set"
        labels = labels97[i]
        k97_rows.append({"copy1_chr": q, "set": group,
                         "supported_gtr_blocks": ",".join(labels) if labels else "NA",
                         "gtr_one_to_one_assignment": G_CHRS[map9[i]] if i in map9 else "NA",
                         "anchor_total": sum(matrix[q]), "anchor_retained": sum(matrix[q][G_CHRS.index(g)] for g in labels),
                         "order_continuity": f"{row_cont[q]:.6f}"})

    rows = [
        {"model": "K=8", "n_groups": 8, "n_complete_pairs": 8, "n_fused_blocks": k8_events,
         "n_split_blocks": "NA", "matched_copy1_chromosomes": 16, "unmatched_copy1_chromosomes": 0,
         "anchor_coverage": f"{k8_retention:.6f}", "order_continuity": f"{k8_order:.6f}",
         "structural_coherence": f"{k8_coherence:.6f}", "event_burden": k8_events, "extra_assignments": k8_extra,
         "score": f"{100*(.6*k8_retention + .2*k8_order + .2*k8_coherence) - 2*k8_events - .5*k8_extra:.6f}",
         "definition": "8 neutral profile-similarity pairs; each pair may cover >1 GTR block"},
        {"model": "K=9-block", "n_groups": 9, "n_complete_pairs": "NA", "n_fused_blocks": k9_events,
         "n_split_blocks": "reported as block multiplicity, not called fission", "matched_copy1_chromosomes": 16,
         "unmatched_copy1_chromosomes": 0, "anchor_coverage": f"{k9_retention:.6f}",
         "order_continuity": f"{k9_order:.6f}", "structural_coherence": f"{k9_coherence:.6f}",
         "event_burden": k9_events, "extra_assignments": k9_extra,
         "score": f"{100*(.6*k9_retention + .2*k9_order + .2*k9_coherence) - 2*k9_events - .5*k9_extra:.6f}",
         "definition": "9 GTR ancestral blocks; current chromosomes can carry multiple supported blocks"},
        {"model": "K=9+7", "n_groups": "9+7", "n_complete_pairs": "NA", "n_fused_blocks": k97_events,
         "n_split_blocks": "NA", "matched_copy1_chromosomes": 16, "unmatched_copy1_chromosomes": 0,
         "anchor_coverage": f"{k97_retention:.6f}", "order_continuity": f"{sum(row_cont.values())/len(row_cont):.6f}",
         "structural_coherence": f"{k97_coherence:.6f}",
         "event_burden": k97_events, "extra_assignments": k97_extra,
         "score": f"{k97_score:.6f}",
         "definition": "enumerated 9/7 partition; 9-set one-to-one to GTR, 7-set fusion-tolerant"},
    ]
    fields = list(rows[0].keys())
    write_tsv(out / "model_comparison.tsv", fields, rows)
    write_tsv(out / "k8_pairs.tsv", list(k8_rows[0].keys()) if k8_rows else [], k8_rows)
    write_tsv(out / "k9_block_assignments.tsv", list(k9_rows[0].keys()) if k9_rows else [], k9_rows)
    write_tsv(out / "k9_plus_7_partition.tsv", list(k97_rows[0].keys()) if k97_rows else [], k97_rows)
    multiplicity_rows = []
    for g in G_CHRS:
        q_supported = [q for q in Q_CHRS if g in row_sets[q]]
        multiplicity_rows.append({
            "gtr_block": g,
            "copy1_chromosomes_with_supported_anchors": ",".join(q_supported) or "NA",
            "n_copy1_chromosomes": len(q_supported),
            "k8_pairs_with_supported_anchors": ",".join(sorted({pair_of_q[q] for q in q_supported})) or "NA",
            "n_k8_pairs": len({pair_of_q[q] for q in q_supported}),
            "interpretation": "raw multiplicity; not by itself fission or gene loss",
        })
    write_tsv(out / "gtr_block_multiplicity.tsv",
              ["gtr_block", "copy1_chromosomes_with_supported_anchors", "n_copy1_chromosomes",
               "k8_pairs_with_supported_anchors", "n_k8_pairs", "interpretation"], multiplicity_rows)
    ranked_partitions.sort(key=lambda x: x[0], reverse=True)
    write_tsv(out / "k9_plus_7_top_partitions.tsv",
              ["rank", "score", "nine_chromosome_set", "seven_chromosome_set", "fit_9_set",
               "fit_7_set", "anchor_coverage", "structural_coherence", "fused_events", "extra_assignments"],
              [{"rank": i + 1, "score": f"{r[0]:.6f}",
                "nine_chromosome_set": ",".join(Q_CHRS[j] for j in r[1]),
                "seven_chromosome_set": ",".join(Q_CHRS[j] for j in q_indices if j not in set(r[1])),
                "fit_9_set": f"{r[2]:.6f}", "fit_7_set": f"{r[3]:.6f}",
                "anchor_coverage": f"{r[4]:.6f}", "structural_coherence": f"{r[7]:.6f}",
                "fused_events": r[5], "extra_assignments": r[6]} for i, r in enumerate(ranked_partitions[:20])])
    write_tsv(out / "chromosome_profiles.tsv",
              ["copy1_chr", "anchor_total", "supported_gtr_blocks", "supported_mass", "order_continuity"] + G_CHRS,
              [{"copy1_chr": q, "anchor_total": sum(matrix[q]), "supported_gtr_blocks": ",".join(row_sets[q]) or "NA",
                "supported_mass": f"{row_mass[q]:.6f}", "order_continuity": f"{row_cont[q]:.6f}",
                **{g: matrix[q][i] for i, g in enumerate(G_CHRS)}} for q in Q_CHRS])

    write_tsv(out / "model_parameters.tsv", ["parameter", "value"], [
        {"parameter": "min_cell_anchors", "value": a.min_cell_anchors},
        {"parameter": "min_cell_fraction", "value": a.min_cell_fraction},
        {"parameter": "min_segment_anchors", "value": a.min_segment_anchors},
        {"parameter": "anchor_pairs_used", "value": total_anchors},
        {"parameter": "order_continuity_definition", "value": "anchor-weighted within-cell score: mean(abs(Spearman rho), dominant step direction)"},
        {"parameter": "score_formula", "value": "100*(0.6*anchor_retention+0.2*order_continuity+0.2*structural_coherence)-2*event_burden-0.5*extra_assignments"},
        {"parameter": "interpretation", "value": "heuristic comparison; no formal likelihood or subgenome naming"},
    ])
    (out / "README.txt").write_text(
        "Neutral model comparison. copy1 labels are chr01-chr16; GTR is Chr01-Chr09 only.\n"
        "K=8 is a profile-based perfect pairing. K=9-block treats nine GTR chromosomes as\n"
        "ancestral blocks and permits a current chromosome to carry multiple supported blocks.\n"
        "K=9+7 enumerates all 9/7 partitions; the nine-member side is matched one-to-one to\n"
        "GTR and the seven-member side is fusion-tolerant. The composite score is heuristic.\n",
        encoding="utf-8")
    print(f"[OK] Neutral model comparison: {out}")
    print(f"[OK] anchors={total_anchors}; K8_pairs=8; K97_9set={','.join(Q_CHRS[i] for i in sorted(nine))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
