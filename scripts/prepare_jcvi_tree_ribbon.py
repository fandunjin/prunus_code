#!/usr/bin/env python3
"""Prepare filtered, ordered and highlighted JCVI karyotype inputs."""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import OrderedDict, defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Block:
    fields: list[str]
    qseq: str
    rseq: str
    qlo: int
    qhi: int
    rlo: int
    rhi: int
    score: int
    orientation: str


@dataclass(frozen=True)
class GenePair:
    qgene: str
    rgene: str
    score: float


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-root", required=True, type=Path)
    p.add_argument("--genome-dir", required=True, type=Path)
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--species", nargs="+", required=True)
    p.add_argument("--pzi-label", default="pzi_copy1")
    p.add_argument("--min-chr-mb", type=float, default=5.0)
    p.add_argument("--highlight-ordinals", nargs="+", type=int, default=[12, 15])
    p.add_argument("--highlight-colors", nargs="+", default=["#E69F00", "#009E73"])
    p.add_argument("--min-background-score", type=int, default=0)
    p.add_argument("--min-chr-pair-score", type=int, default=0)
    p.add_argument("--min-chr-pair-fraction", type=float, default=0.0)
    p.add_argument("--max-background-blocks-per-chr-pair", type=int, default=0)
    p.add_argument("--fixed-order-species", nargs="*", default=[])
    p.add_argument("--highlight-max-gene-gap", type=int, default=30)
    p.add_argument("--highlight-pad-genes", type=int, default=1)
    p.add_argument("--order-sweeps", type=int, default=20)
    p.add_argument("--local-swap-passes", type=int, default=10)
    p.add_argument("--no-auto-order", action="store_true")
    p.add_argument("--no-auto-orient", action="store_true")
    return p.parse_args()


def read_seqids(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8", errors="replace")
    seqids = [x for x in re.split(r"[\s,]+", text.strip()) if x]
    if not seqids:
        raise SystemExit(f"Empty seqids file: {path}")
    return [x[:-1] if x.endswith("-") else x for x in seqids]


def read_bed(
    path: Path,
) -> tuple[
    list[str],
    dict[str, str],
    list[str],
    dict[str, int],
    dict[str, int],
    dict[str, list[str]],
]:
    lines: list[str] = []
    gene_to_seq: dict[str, str] = {}
    gene_rank: dict[str, int] = {}
    seqids: list[str] = []
    gene_counts: dict[str, int] = defaultdict(int)
    genes_by_seq: dict[str, list[str]] = defaultdict(list)
    seen_seqids = set()
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise SystemExit(f"BED has fewer than six columns: {path}:{line_number}")
            seqid, gene = fields[0], fields[3]
            if gene in gene_to_seq:
                raise SystemExit(f"Duplicate BED gene ID {gene} in {path}")
            gene_to_seq[gene] = seqid
            gene_rank[gene] = gene_counts[seqid]
            gene_counts[seqid] += 1
            genes_by_seq[seqid].append(gene)
            if seqid not in seen_seqids:
                seqids.append(seqid)
                seen_seqids.add(seqid)
            lines.append(line.rstrip("\n"))
    if not lines:
        raise SystemExit(f"Empty BED file: {path}")
    return lines, gene_to_seq, seqids, gene_rank, dict(gene_counts), dict(genes_by_seq)


def fasta_lengths(path: Path, wanted: set[str]) -> dict[str, int]:
    lengths: dict[str, int] = {}
    current = None
    length = 0
    if not path.is_file():
        return lengths
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if current in wanted:
                    lengths[current] = length
                current = line[1:].strip().split()[0]
                length = 0
            elif current in wanted:
                length += len(line.strip())
    if current in wanted:
        lengths[current] = length
    return lengths


def bed_lengths(lines: list[str]) -> dict[str, int]:
    lengths: dict[str, int] = {}
    for line in lines:
        fields = line.split("\t")
        seqid = fields[0]
        end = int(fields[2])
        lengths[seqid] = max(lengths.get(seqid, 0), end)
    return lengths


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def parse_simple(path: Path, query: str, ref: str, data: dict) -> list[Block]:
    blocks: list[Block] = []
    qmap = data[query]["gene_to_seq"]
    rmap = data[ref]["gene_to_seq"]
    qrank = data[query]["gene_rank"]
    rrank = data[ref]["gene_rank"]
    qkeep = data[query]["keep_set"]
    rkeep = data[ref]["keep_set"]
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line or line.startswith("##") or line.startswith("StartGeneA"):
                continue
            if "*" in line:
                line = line.split("*", 1)[1].strip()
            elif line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 6:
                print(f"WARNING: skip non-six-column simple row: {path}:{line_number}")
                continue
            a, b, c, d, score_text, orientation = fields
            if not (a in qmap and b in qmap and c in rmap and d in rmap):
                if c in qmap and d in qmap and a in rmap and b in rmap:
                    fields = [c, d, a, b, score_text, orientation]
                    a, b, c, d = fields[:4]
                else:
                    continue
            qseq = qmap[a]
            rseq = rmap[c]
            if qmap[b] != qseq or rmap[d] != rseq:
                continue
            if qseq not in qkeep or rseq not in rkeep:
                continue
            try:
                score = int(float(score_text))
            except ValueError:
                print(f"WARNING: invalid block score: {path}:{line_number}")
                continue
            qlo, qhi = sorted((qrank[a], qrank[b]))
            rlo, rhi = sorted((rrank[c], rrank[d]))
            blocks.append(
                Block(fields, qseq, rseq, qlo, qhi, rlo, rhi, score, orientation)
            )
    return blocks


def optimize_tracks(
    species: list[str],
    data: dict,
    pair_blocks: dict,
    auto_order: bool,
    auto_orient: bool,
    fixed_order_species: set[str],
    order_sweeps: int,
    local_swap_passes: int,
):
    orders: dict[str, list[str]] = {species[0]: list(data[species[0]]["keep"])}
    reversed_seqids: dict[str, set[str]] = {species[0]: set()}
    barycenters: dict[tuple[str, str], float | None] = {}

    for i in range(len(species) - 1):
        query, ref = species[i], species[i + 1]
        qorder = orders.get(query, list(data[query]["keep"]))
        qindex = {seqid: pos for pos, seqid in enumerate(qorder)}
        qrev = reversed_seqids.get(query, set())
        blocks = pair_blocks.get((query, ref), [])
        pair_scores: dict[tuple[str, str], int] = defaultdict(int)
        for block in blocks:
            pair_scores[(block.qseq, block.rseq)] += max(1, block.score)

        original = list(data[ref]["keep"])
        assignments: dict[str, str] = {}
        if len(qorder) == len(original):
            candidates = sorted(
                ((score, qseq, rseq) for (qseq, rseq), score in pair_scores.items()),
                reverse=True,
            )
            used_q: set[str] = set()
            used_r: set[str] = set()
            for _, qseq, rseq in candidates:
                if qseq not in used_q and rseq not in used_r:
                    assignments[rseq] = qseq
                    used_q.add(qseq)
                    used_r.add(rseq)
        for rseq in original:
            if rseq in assignments:
                continue
            candidates = [
                (score, qseq)
                for (qseq, candidate_r), score in pair_scores.items()
                if candidate_r == rseq
            ]
            if candidates:
                assignments[rseq] = max(candidates)[1]

        weighted_positions: dict[str, list[tuple[float, int]]] = defaultdict(list)
        orientation_votes: dict[str, dict[bool, int]] = defaultdict(lambda: defaultdict(int))
        for block in blocks:
            if assignments.get(block.rseq) != block.qseq:
                continue
            count = max(1, data[query]["gene_counts"].get(block.qseq, 1) - 1)
            fraction = ((block.qlo + block.qhi) / 2.0) / count
            if block.qseq in qrev:
                fraction = 1.0 - fraction
            position = qindex.get(block.qseq, len(qindex)) + fraction
            weight = max(1, block.score)
            weighted_positions[block.rseq].append((position, weight))
            candidate_reverse = (block.qseq in qrev) ^ (block.orientation == "-")
            orientation_votes[block.rseq][candidate_reverse] += weight

        for seqid in original:
            values = weighted_positions.get(seqid, [])
            if values:
                barycenters[(ref, seqid)] = sum(x * w for x, w in values) / sum(w for _, w in values)
            else:
                barycenters[(ref, seqid)] = None

        if auto_order and ref not in fixed_order_species:
            original_index = {seqid: pos for pos, seqid in enumerate(original)}
            matched = [seqid for seqid in original if barycenters[(ref, seqid)] is not None]
            unmatched = [seqid for seqid in original if barycenters[(ref, seqid)] is None]
            matched.sort(key=lambda seqid: (barycenters[(ref, seqid)], original_index[seqid]))
            orders[ref] = matched + unmatched
        else:
            orders[ref] = original

        reversed_seqids[ref] = set()
        if auto_orient and ref not in fixed_order_species:
            for seqid in orders[ref]:
                votes = orientation_votes.get(seqid, {})
                if votes.get(True, 0) > votes.get(False, 0):
                    reversed_seqids[ref].add(seqid)

    for label in species:
        orders.setdefault(label, list(data[label]["keep"]))
        reversed_seqids.setdefault(label, set())

    pair_scores: dict[tuple[str, str], dict[tuple[str, str], int]] = {}
    for pair, blocks in pair_blocks.items():
        scores: dict[tuple[str, str], int] = defaultdict(int)
        for block in blocks:
            scores[(block.qseq, block.rseq)] += max(1, block.score)
        pair_scores[pair] = dict(scores)

    def pair_crossing_cost(
        query: str, ref: str, qorder: list[str], rorder: list[str]
    ) -> float:
        scores = pair_scores.get((query, ref), {})
        qindex = {seqid: pos for pos, seqid in enumerate(qorder)}
        rindex = {seqid: pos for pos, seqid in enumerate(rorder)}
        edges = [
            (qindex[qseq], rindex[rseq], math.sqrt(score))
            for (qseq, rseq), score in scores.items()
            if qseq in qindex and rseq in rindex
        ]
        total = 0.0
        for left in range(len(edges)):
            q1, r1, w1 = edges[left]
            for right in range(left + 1, len(edges)):
                q2, r2, w2 = edges[right]
                if (q1 - q2) * (r1 - r2) < 0:
                    total += w1 * w2
        return total

    def crossing_cost(candidate_orders: dict[str, list[str]]) -> float:
        return sum(
            pair_crossing_cost(
                species[i],
                species[i + 1],
                candidate_orders[species[i]],
                candidate_orders[species[i + 1]],
            )
            for i in range(len(species) - 1)
        )

    def layer_crossing_cost(layer_index: int, candidate_order: list[str]) -> float:
        label = species[layer_index]
        total = 0.0
        if layer_index > 0:
            upper = species[layer_index - 1]
            total += pair_crossing_cost(upper, label, orders[upper], candidate_order)
        if layer_index < len(species) - 1:
            lower = species[layer_index + 1]
            total += pair_crossing_cost(label, lower, candidate_order, orders[lower])
        return total

    def neighbor_barycenter(layer_index: int, seqid: str) -> float | None:
        label = species[layer_index]
        values: list[tuple[float, float]] = []
        if layer_index > 0:
            upper = species[layer_index - 1]
            upper_order = orders[upper]
            upper_index = {x: i / max(1, len(upper_order) - 1) for i, x in enumerate(upper_order)}
            for (upper_seqid, current_seqid), score in pair_scores.get((upper, label), {}).items():
                if current_seqid == seqid and upper_seqid in upper_index:
                    values.append((upper_index[upper_seqid], math.sqrt(score)))
        if layer_index < len(species) - 1:
            lower = species[layer_index + 1]
            lower_order = orders[lower]
            lower_index = {x: i / max(1, len(lower_order) - 1) for i, x in enumerate(lower_order)}
            for (current_seqid, lower_seqid), score in pair_scores.get((label, lower), {}).items():
                if current_seqid == seqid and lower_seqid in lower_index:
                    values.append((lower_index[lower_seqid], math.sqrt(score)))
        if not values:
            return None
        return sum(position * weight for position, weight in values) / sum(
            weight for _, weight in values
        )

    initial_cost = crossing_cost(orders)
    best_cost = initial_cost
    best_orders = {label: list(order) for label, order in orders.items()}
    if auto_order:
        for _ in range(max(0, order_sweeps)):
            changed = False
            for indices in (range(1, len(species)), range(len(species) - 2, -1, -1)):
                for layer_index in indices:
                    label = species[layer_index]
                    if label in fixed_order_species:
                        continue
                    original_index = {seqid: pos for pos, seqid in enumerate(orders[label])}
                    scored = []
                    unscored = []
                    for seqid in orders[label]:
                        center = neighbor_barycenter(layer_index, seqid)
                        barycenters[(label, seqid)] = center
                        if center is None:
                            unscored.append(seqid)
                        else:
                            scored.append((center, original_index[seqid], seqid))
                    proposed = [x[2] for x in sorted(scored)] + unscored
                    if proposed != orders[label]:
                        orders[label] = proposed
                        changed = True
            cost = crossing_cost(orders)
            if cost < best_cost:
                best_cost = cost
                best_orders = {label: list(order) for label, order in orders.items()}
            if not changed:
                break
        orders = {label: list(order) for label, order in best_orders.items()}
        for _ in range(max(0, local_swap_passes)):
            changed = False
            for layer_index, label in enumerate(species):
                if label in fixed_order_species or len(orders[label]) < 2:
                    continue
                baseline = layer_crossing_cost(layer_index, orders[label])
                best_layer_cost = baseline
                best_layer_order = orders[label]
                for left in range(len(orders[label]) - 1):
                    for right in range(left + 1, len(orders[label])):
                        candidate = list(orders[label])
                        candidate[left], candidate[right] = candidate[right], candidate[left]
                        cost = layer_crossing_cost(layer_index, candidate)
                        if cost + 1e-9 < best_layer_cost:
                            best_layer_cost = cost
                            best_layer_order = candidate
                if best_layer_order != orders[label]:
                    orders[label] = best_layer_order
                    changed = True
            cost = crossing_cost(orders)
            if cost < best_cost:
                best_cost = cost
                best_orders = {label: list(order) for label, order in orders.items()}
            if not changed:
                break
    orders = best_orders

    # Refine whole-chromosome orientation against both adjacent tracks. A
    # chromosome flip toggles every incident block equally, so this preserves
    # all internal inversions and rearrangements while reducing avoidable
    # direction conflicts between neighboring species.
    if auto_orient:
        for _ in range(max(1, order_sweeps)):
            changed = False
            for layer_index, label in enumerate(species):
                if label in fixed_order_species:
                    continue
                for seqid in orders[label]:
                    agreeing = 0.0
                    disagreeing = 0.0
                    incident: list[tuple[str, str, Block]] = []
                    if layer_index > 0:
                        upper = species[layer_index - 1]
                        incident.extend(
                            (upper, label, block)
                            for block in pair_blocks.get((upper, label), [])
                            if block.rseq == seqid
                        )
                    if layer_index < len(species) - 1:
                        lower = species[layer_index + 1]
                        incident.extend(
                            (label, lower, block)
                            for block in pair_blocks.get((label, lower), [])
                            if block.qseq == seqid
                        )
                    for query, ref, block in incident:
                        displayed_reverse = (
                            (block.qseq in reversed_seqids[query])
                            ^ (block.rseq in reversed_seqids[ref])
                            ^ (block.orientation == "-")
                        )
                        weight = math.sqrt(max(1, block.score))
                        if displayed_reverse:
                            disagreeing += weight
                        else:
                            agreeing += weight
                    if disagreeing > agreeing + 1e-9:
                        if seqid in reversed_seqids[label]:
                            reversed_seqids[label].remove(seqid)
                        else:
                            reversed_seqids[label].add(seqid)
                        changed = True
            if not changed:
                break
    return orders, reversed_seqids, barycenters, initial_cost, best_cost


def parse_anchor_pairs(path: Path, query: str, ref: str, data: dict) -> list[GenePair]:
    pairs: list[GenePair] = []
    qmap = data[query]["gene_to_seq"]
    rmap = data[ref]["gene_to_seq"]
    qkeep = data[query]["keep_set"]
    rkeep = data[ref]["keep_set"]
    if not path.is_file():
        return pairs
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.split()
            if len(fields) < 2:
                continue
            qgene, rgene = fields[:2]
            if qgene not in qmap or rgene not in rmap:
                if rgene in qmap and qgene in rmap:
                    qgene, rgene = rgene, qgene
                else:
                    continue
            if qmap[qgene] not in qkeep or rmap[rgene] not in rkeep:
                continue
            try:
                score = float(fields[2]) if len(fields) > 2 else 1.0
            except ValueError:
                score = 1.0
            pairs.append(GenePair(qgene, rgene, score))
    return pairs


def filter_background_blocks(
    pair_blocks: OrderedDict[tuple[str, str], list[Block]],
    min_block_score: int,
    min_pair_score: int,
    min_pair_fraction: float,
    max_blocks_per_pair: int,
) -> tuple[OrderedDict[tuple[str, str], list[Block]], list[dict[str, object]]]:
    filtered: OrderedDict[tuple[str, str], list[Block]] = OrderedDict()
    summary: list[dict[str, object]] = []
    for pair, blocks in pair_blocks.items():
        score_pass = [block for block in blocks if block.score >= min_block_score]
        pair_totals: dict[tuple[str, str], int] = defaultdict(int)
        for block in score_pass:
            pair_totals[(block.qseq, block.rseq)] += block.score
        best_by_query: dict[str, int] = defaultdict(int)
        best_by_ref: dict[str, int] = defaultdict(int)
        for (qseq, rseq), score in pair_totals.items():
            best_by_query[qseq] = max(best_by_query[qseq], score)
            best_by_ref[rseq] = max(best_by_ref[rseq], score)
        retained_pairs = {
            (qseq, rseq)
            for (qseq, rseq), score in pair_totals.items()
            if score >= min_pair_score
            and score >= min_pair_fraction * best_by_query[qseq]
            and score >= min_pair_fraction * best_by_ref[rseq]
        }
        grouped: dict[tuple[str, str], list[Block]] = defaultdict(list)
        for block in score_pass:
            chromosome_pair = (block.qseq, block.rseq)
            if chromosome_pair in retained_pairs:
                grouped[chromosome_pair].append(block)
        kept: list[Block] = []
        for chromosome_pair, chromosome_blocks in grouped.items():
            chromosome_blocks.sort(key=lambda block: block.score, reverse=True)
            if max_blocks_per_pair > 0:
                chromosome_blocks = chromosome_blocks[:max_blocks_per_pair]
            kept.extend(chromosome_blocks)
        filtered[pair] = kept
        summary.append(
            {
                "query": pair[0],
                "ref": pair[1],
                "input_blocks": len(blocks),
                "score_pass_blocks": len(score_pass),
                "retained_chr_pairs": len(retained_pairs),
                "displayed_blocks": len(kept),
            }
        )
    return filtered, summary


def group_highlight_pairs(
    pairs: list[GenePair],
    query: str,
    ref: str,
    data: dict,
    max_gap: int,
    pad: int,
) -> list[str]:
    qmap = data[query]["gene_to_seq"]
    rmap = data[ref]["gene_to_seq"]
    qrank = data[query]["gene_rank"]
    rrank = data[ref]["gene_rank"]
    qgenes_by_seq = data[query]["genes_by_seq"]
    rgenes_by_seq = data[ref]["genes_by_seq"]
    by_chromosomes: dict[tuple[str, str], list[GenePair]] = defaultdict(list)
    seen: set[tuple[str, str]] = set()
    for pair in pairs:
        key = (pair.qgene, pair.rgene)
        if key in seen:
            continue
        seen.add(key)
        by_chromosomes[(qmap[pair.qgene], rmap[pair.rgene])].append(pair)

    rows: list[str] = []
    for (qseq, rseq), chromosome_pairs in by_chromosomes.items():
        chromosome_pairs.sort(key=lambda pair: (qrank[pair.qgene], rrank[pair.rgene]))
        groups: list[list[GenePair]] = []
        directions: list[int] = []
        for pair in chromosome_pairs:
            if not groups:
                groups.append([pair])
                directions.append(0)
                continue
            previous = groups[-1][-1]
            qdelta = qrank[pair.qgene] - qrank[previous.qgene]
            rdelta = rrank[pair.rgene] - rrank[previous.rgene]
            direction = 1 if rdelta > 0 else -1 if rdelta < 0 else directions[-1]
            same_direction = directions[-1] in (0, direction)
            if qdelta <= max_gap and abs(rdelta) <= max_gap and same_direction:
                groups[-1].append(pair)
                if directions[-1] == 0 and direction:
                    directions[-1] = direction
            else:
                groups.append([pair])
                directions.append(0)

        qordered = qgenes_by_seq[qseq]
        rordered = rgenes_by_seq[rseq]
        for group, direction in zip(groups, directions):
            qranks = [qrank[pair.qgene] for pair in group]
            rranks = [rrank[pair.rgene] for pair in group]
            qlo = max(0, min(qranks) - pad)
            qhi = min(len(qordered) - 1, max(qranks) + pad)
            rlo = max(0, min(rranks) - pad)
            rhi = min(len(rordered) - 1, max(rranks) + pad)
            if direction == 0 and len(group) > 1:
                direction = 1 if rranks[-1] >= rranks[0] else -1
            orientation = "-" if direction < 0 else "+"
            a_gene, b_gene = qordered[qlo], qordered[qhi]
            if orientation == "-":
                c_gene, d_gene = rordered[rhi], rordered[rlo]
            else:
                c_gene, d_gene = rordered[rlo], rordered[rhi]
            rows.append(
                "\t".join(
                    (a_gene, b_gene, c_gene, d_gene, str(len(group)), orientation)
                )
            )
    return rows


def cascade_gene_highlights(
    species: list[str],
    data: dict,
    pair_gene_pairs: dict,
    retained_pair_blocks: dict,
    pzi_label: str,
    target_ordinals: list[int],
    max_gap: int,
    pad: int,
) -> tuple[dict[int, str], dict[tuple[str, str], list[tuple[int, str]]], list[dict[str, object]]]:
    pzi = data[pzi_label]
    target_seqids: dict[int, str] = {}
    active: dict[int, set[str]] = {}
    # Resolve targets from the original input .seqids order before any display
    # reordering or whole-chromosome reversal is applied.
    for ordinal in target_ordinals:
        if ordinal < 1 or ordinal > len(pzi["all_seqids"]):
            print(f"WARNING: pzi chromosome ordinal is out of range: {ordinal}")
            continue
        seqid = pzi["all_seqids"][ordinal - 1]
        if seqid not in pzi["keep_set"]:
            print(f"WARNING: highlighted pzi chromosome was removed by length filter: {seqid}")
            continue
        target_seqids[ordinal] = seqid
        active[ordinal] = set(pzi["genes_by_seq"].get(seqid, []))

    summary: list[dict[str, object]] = []
    highlight_rows: dict[tuple[str, str], list[tuple[int, str]]] = defaultdict(list)
    for i in range(len(species) - 1):
        query, ref = species[i], species[i + 1]
        gene_pairs = pair_gene_pairs.get((query, ref), [])
        qmap = data[query]["gene_to_seq"]
        rmap = data[ref]["gene_to_seq"]
        qrank = data[query]["gene_rank"]
        rrank = data[ref]["gene_rank"]
        retained_by_chromosomes: dict[tuple[str, str], list[Block]] = defaultdict(list)
        for block in retained_pair_blocks.get((query, ref), []):
            retained_by_chromosomes[(block.qseq, block.rseq)].append(block)

        def retained_gene_pair(pair: GenePair) -> bool:
            qseq = qmap[pair.qgene]
            rseq = rmap[pair.rgene]
            qr = qrank[pair.qgene]
            rr = rrank[pair.rgene]
            return any(
                block.qlo <= qr <= block.qhi and block.rlo <= rr <= block.rhi
                for block in retained_by_chromosomes.get((qseq, rseq), [])
            )

        next_active: dict[int, set[str]] = {}
        for ordinal in target_seqids:
            selected = [
                pair
                for pair in gene_pairs
                if pair.qgene in active.get(ordinal, set()) and retained_gene_pair(pair)
            ]
            next_active[ordinal] = {pair.rgene for pair in selected}
            rows = group_highlight_pairs(selected, query, ref, data, max_gap, pad)
            highlight_rows[(query, ref)].extend((ordinal, row) for row in rows)
            summary.append(
                {
                    "query": query,
                    "ref": ref,
                    "ordinal": ordinal,
                    "gene_pairs": len(selected),
                    "blocks": len(rows),
                    "ref_genes": len(next_active[ordinal]),
                }
            )
        active = next_active
    return target_seqids, highlight_rows, summary


def main() -> None:
    a = args()
    species = list(dict.fromkeys(a.species))
    if a.pzi_label not in species or species[0] != a.pzi_label:
        raise SystemExit("The highlighted pzi label must be the first displayed species")
    colors = {
        ordinal: a.highlight_colors[i % len(a.highlight_colors)]
        for i, ordinal in enumerate(a.highlight_ordinals)
    }
    out = a.outdir
    filtered_dir = out / "filtered"
    figures_dir = out / "figures"
    filtered_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    min_bp = int(a.min_chr_mb * 1_000_000)

    data: dict[str, dict[str, object]] = {}
    for label in species:
        input_dir = a.input_root / label
        bed_path = input_dir / f"{label}.bed"
        seqids_path = input_dir / f"{label}.seqids"
        genome_path = a.genome_dir / f"{label}.genome.fa"
        if not bed_path.is_file() or not seqids_path.is_file():
            raise SystemExit(f"Missing JCVI input for {label}: {bed_path}")
        (
            raw_lines,
            gene_to_seq,
            bed_seqids,
            gene_rank,
            gene_counts,
            genes_by_seq,
        ) = read_bed(bed_path)
        declared_seqids = read_seqids(seqids_path)
        bed_seqid_set = set(bed_seqids)
        seqids = [x for x in declared_seqids if x in bed_seqid_set]
        lengths = fasta_lengths(genome_path, set(seqids))
        source = "genome_fasta"
        if set(lengths) != set(seqids):
            lengths = bed_lengths(raw_lines)
            source = "BED_max_gene_end"
            print(f"WARNING: {label} genome FASTA is missing chromosome IDs; using BED lengths")
        keep = [seqid for seqid in seqids if lengths.get(seqid, 0) >= min_bp]
        if not keep:
            raise SystemExit(f"No chromosomes >= {a.min_chr_mb} Mb for {label}")
        keep_set = set(keep)
        kept_lines = [line for line in raw_lines if line.split("\t", 1)[0] in keep_set]
        display_alias = {
            seqid: f"{label}_chr{ordinal}"
            for ordinal, seqid in enumerate(seqids, 1)
        }
        display_lines = []
        for line in kept_lines:
            fields = line.split("\t")
            fields[0] = display_alias[fields[0]]
            display_lines.append("\t".join(fields))
        write_text(filtered_dir / f"{label}.bed", "\n".join(kept_lines) + "\n")
        write_text(filtered_dir / f"{label}.seqids", ",".join(keep) + "\n")
        write_text(figures_dir / f"{label}.filtered.bed", "\n".join(display_lines) + "\n")
        data[label] = {
            "all_seqids": seqids,
            "keep": keep,
            "keep_set": keep_set,
            "lengths": lengths,
            "gene_to_seq": gene_to_seq,
            "gene_rank": gene_rank,
            "gene_counts": gene_counts,
            "genes_by_seq": genes_by_seq,
            "display_alias": display_alias,
            "source": source,
        }

    pair_blocks: OrderedDict[tuple[str, str], list[Block]] = OrderedDict()
    pair_gene_pairs: OrderedDict[tuple[str, str], list[GenePair]] = OrderedDict()
    with a.manifest.open("r", encoding="utf-8", errors="replace") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row or row[0].startswith("#") or len(row) < 3:
                continue
            query, ref, source_path = row[:3]
            raw_anchor_path = row[3] if len(row) > 3 else "NA"
            if query not in data or ref not in data:
                continue
            source = Path(source_path)
            if not source.is_file():
                print(f"WARNING: missing anchors: {source}")
                continue
            blocks = parse_simple(source, query, ref, data)
            if blocks:
                pair_blocks[(query, ref)] = blocks
            else:
                print(f"WARNING: no retained anchors for {query} vs {ref}")
            raw_anchor = Path(raw_anchor_path)
            gene_pairs = parse_anchor_pairs(raw_anchor, query, ref, data)
            if gene_pairs:
                pair_gene_pairs[(query, ref)] = gene_pairs
            else:
                print(f"WARNING: no usable raw gene anchors for {query} vs {ref}: {raw_anchor_path}")

    display_pair_blocks, background_filter_summary = filter_background_blocks(
        pair_blocks,
        a.min_background_score,
        a.min_chr_pair_score,
        a.min_chr_pair_fraction,
        a.max_background_blocks_per_chr_pair,
    )
    orders, reversed_seqids, barycenters, initial_crossings, final_crossings = optimize_tracks(
        species,
        data,
        display_pair_blocks,
        not a.no_auto_order,
        not a.no_auto_orient,
        set(a.fixed_order_species),
        a.order_sweeps,
        a.local_swap_passes,
    )
    target_seqids, highlight_rows, highlight_summary = cascade_gene_highlights(
        species,
        data,
        pair_gene_pairs,
        display_pair_blocks,
        a.pzi_label,
        a.highlight_ordinals,
        a.highlight_max_gene_gap,
        a.highlight_pad_genes,
    )

    anchor_files: OrderedDict[tuple[str, str], str] = OrderedDict()
    for (query, ref), blocks in display_pair_blocks.items():
        destination = figures_dir / f"{query}.{ref}.filtered.anchors.simple"
        output_lines = []
        for block in blocks:
            output_lines.append("\t".join(block.fields))
        for ordinal, row in highlight_rows.get((query, ref), []):
            output_lines.append(f"{colors.get(ordinal, a.highlight_colors[0])}*{row}")
        if output_lines:
            write_text(destination, "\n".join(output_lines) + "\n")
            anchor_files[(query, ref)] = destination.name

    y_values = [0.97 - i * (0.89 / max(1, len(species) - 1)) for i in range(len(species))]
    seqids_lines = []
    layout_lines = ["# y, xstart, xend, rotation, color, label, va, bed"]
    basepair_layout_lines = ["# y, xstart, xend, rotation, color, label, va, bed"]
    displayed_bp = {
        label: sum(data[label]["lengths"].get(seqid, 0) for seqid in data[label]["keep"])
        for label in species
    }
    max_displayed_bp = max(displayed_bp.values())
    track_center = 0.50
    max_track_width = 0.84
    for i, label in enumerate(species):
        displayed = [
            data[label]["display_alias"][seqid]
            + ("-" if seqid in reversed_seqids[label] else "")
            for seqid in orders[label]
        ]
        seqids_lines.append(",".join(displayed))
        color = "#4477AA" if label == a.pzi_label else "#CC6677"
        display_label = "pmu" if label == "pmu_refseq" else label
        va = "top" if i == 0 else "center"
        layout_lines.append(
            f" {y_values[i]:.4f}, .08, .92, 0, {color}, {display_label}, {va}, {label}.filtered.bed"
        )
        relative_width = displayed_bp[label] / max_displayed_bp
        track_width = max_track_width * relative_width
        basepair_xstart = track_center - track_width / 2.0
        basepair_xend = track_center + track_width / 2.0
        basepair_layout_lines.append(
            f" {y_values[i]:.4f}, {basepair_xstart:.4f}, {basepair_xend:.4f}, 0, {color}, "
            f"{display_label}, {va}, {label}.filtered.bed"
        )

    layout_lines.append("# edges")
    basepair_layout_lines.append("# edges")
    for i in range(len(species) - 1):
        pair = (species[i], species[i + 1])
        anchor_name = anchor_files.get(pair)
        if anchor_name:
            layout_lines.append(f" e, {i}, {i + 1}, {anchor_name}")
            basepair_layout_lines.append(f" e, {i}, {i + 1}, {anchor_name}")

    write_text(figures_dir / "tree_order_pairwise.seqids", "\n".join(seqids_lines) + "\n")
    write_text(figures_dir / "tree_order_pairwise.layout", "\n".join(layout_lines) + "\n")
    write_text(
        figures_dir / "tree_order_pairwise_basepair.layout",
        "\n".join(basepair_layout_lines) + "\n",
    )
    with (out / "basepair_track_scale.tsv").open("w", encoding="utf-8") as handle:
        handle.write("species\tdisplayed_bp\tdisplayed_mb\trelative_to_max\txstart\txend\n")
        for label in species:
            relative_width = displayed_bp[label] / max_displayed_bp
            track_width = max_track_width * relative_width
            basepair_xstart = track_center - track_width / 2.0
            basepair_xend = track_center + track_width / 2.0
            handle.write(
                f"{label}\t{displayed_bp[label]}\t{displayed_bp[label] / 1_000_000:.6f}\t"
                f"{relative_width:.6f}\t{basepair_xstart:.4f}\t{basepair_xend:.4f}\n"
            )
    with (out / "highlighted_pzi_chromosomes.tsv").open("w", encoding="utf-8") as handle:
        handle.write("chromosome\tseqid\tcolor\n")
        for ordinal, seqid in target_seqids.items():
            handle.write(f"{ordinal}\t{seqid}\t{colors[ordinal]}\n")
    with (out / "highlighted_blocks.tsv").open("w", encoding="utf-8") as handle:
        handle.write("query\tref\tchromosome\thighlighted_gene_pairs\thighlighted_blocks\tref_genes\n")
        for row in highlight_summary:
            handle.write(
                f"{row['query']}\t{row['ref']}\t{row['ordinal']}\t{row['gene_pairs']}\t"
                f"{row['blocks']}\t{row['ref_genes']}\n"
            )
    with (out / "chromosome_order_orientation.tsv").open("w", encoding="utf-8") as handle:
        handle.write("species\tposition\tdisplay_seqid\toriginal_seqid\torientation\tbarycenter\n")
        for label in species:
            for position, seqid in enumerate(orders[label], 1):
                bary = barycenters.get((label, seqid))
                bary_text = "NA" if bary is None else f"{bary:.6f}"
                orientation = "-" if seqid in reversed_seqids[label] else "+"
                display_seqid = data[label]["display_alias"][seqid]
                handle.write(
                    f"{label}\t{position}\t{display_seqid}\t{seqid}\t{orientation}\t{bary_text}\n"
                )
    with (out / "chromosome_display_labels.tsv").open("w", encoding="utf-8") as handle:
        handle.write("species\toriginal_ordinal\tdisplay_seqid\toriginal_seqid\tretained\n")
        for label in species:
            for ordinal, seqid in enumerate(data[label]["all_seqids"], 1):
                handle.write(
                    f"{label}\t{ordinal}\t{data[label]['display_alias'][seqid]}\t{seqid}\t"
                    f"{int(seqid in data[label]['keep_set'])}\n"
                )
    with (out / "filtered_chromosomes.tsv").open("w", encoding="utf-8") as handle:
        handle.write("species\tkept_chromosomes\tmin_mb\tlength_source\n")
        for label in species:
            handle.write(
                f"{label}\t{len(data[label]['keep'])}\t{a.min_chr_mb}\t{data[label]['source']}\n"
            )
    with (out / "background_filtering.tsv").open("w", encoding="utf-8") as handle:
        handle.write(
            "query\tref\tinput_blocks\tscore_pass_blocks\tretained_chr_pairs\tdisplayed_blocks\n"
        )
        for row in background_filter_summary:
            handle.write(
                f"{row['query']}\t{row['ref']}\t{row['input_blocks']}\t"
                f"{row['score_pass_blocks']}\t{row['retained_chr_pairs']}\t"
                f"{row['displayed_blocks']}\n"
            )
    reduction = 0.0 if initial_crossings == 0 else 100.0 * (
        initial_crossings - final_crossings
    ) / initial_crossings
    with (out / "order_optimization.tsv").open("w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"initial_weighted_crossings\t{initial_crossings:.6f}\n")
        handle.write(f"final_weighted_crossings\t{final_crossings:.6f}\n")
        handle.write(f"reduction_percent\t{reduction:.6f}\n")
        handle.write(f"sweeps\t{a.order_sweeps}\n")
        handle.write(f"local_swap_passes\t{a.local_swap_passes}\n")

    print(f"Prepared tree-order ribbon inputs: {figures_dir}")
    print(f"Adjacent pairwise edges retained: {len(anchor_files)}")
    print(f"Highlighted pzi chromosomes: {target_seqids}")
    print(f"Auto order: {not a.no_auto_order}; auto orientation: {not a.no_auto_orient}")
    print(f"Order sweeps: {a.order_sweeps}")
    print(f"Local swap passes: {a.local_swap_passes}")
    print(f"Weighted chromosome crossings: {initial_crossings:.3f} -> {final_crossings:.3f}")
    print(f"Fixed-order species: {a.fixed_order_species}")
    print(f"Minimum background block score: {a.min_background_score}")
    print(f"Minimum chromosome-pair score: {a.min_chr_pair_score}")
    print(f"Minimum chromosome-pair fraction: {a.min_chr_pair_fraction}")
    print(
        "Maximum background blocks per chromosome pair: "
        f"{a.max_background_blocks_per_chr_pair or 'unlimited'}"
    )


if __name__ == "__main__":
    main()
