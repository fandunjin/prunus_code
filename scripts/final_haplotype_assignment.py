#!/usr/bin/env python3
import argparse
import csv
import itertools
import math
from collections import Counter, defaultdict
from pathlib import Path


HAP_ORDER = ["hapA1", "hapA2", "hapB1", "hapB2"]

COPY1_HAPS = {"hapA1", "hapB1"}
COPY2_HAPS = {"hapA2", "hapB2"}


def natural_key(text):
    parts = []
    buf = ""
    digit = False
    for char in str(text):
        is_digit = char.isdigit()
        if buf and is_digit != digit:
            parts.append(int(buf) if digit else buf)
            buf = char
        else:
            buf += char
        digit = is_digit
    if buf:
        parts.append(int(buf) if digit else buf)
    return parts


def read_tsv(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fieldnames):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def short_evidence_summary(text):
    keep = []
    for part in str(text).split("; "):
        if part.startswith("copy_balance="):
            keep.append(part)
        elif part.startswith("corrected_pair_support="):
            keep.append(part)
        elif part.startswith("ppa:"):
            keep.append(part)
        elif part.startswith("n8_refs="):
            keep.append(part)
        elif part.startswith("outgroups="):
            keep.append(part)
        elif part.startswith("hard_constraint="):
            keep.append(part)
    return "; ".join(keep)


def parse_float(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def parse_int(value, default=0):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def probability_like(score):
    return 1.0 / (1.0 + math.exp(-8.0 * (score - 0.50)))


def confidence_label(score, review_flag):
    if review_flag == "review_required":
        return "medium_review"
    if score >= 0.86:
        return "high"
    if score >= 0.72:
        return "medium"
    return "review"


def read_corrected_pairs(path):
    rows = read_tsv(path)
    pairs = []
    copy_lookup = {}
    pair_lookup = {}
    for row in rows:
        copy1_chr = row.get("copy1_chr", "") or row.get("copy1_source", "")
        copy2_chr = row.get("copy2_chr", "") or row.get("copy2_source", "")
        pair_id = row.get("pair_id", "")
        if not copy1_chr or not copy2_chr:
            continue
        pair = {
            "pair_id": pair_id,
            "copy1_chr": copy1_chr,
            "copy2_chr": copy2_chr,
            "original_copy1_chr": row.get("original_copy1_chr", ""),
            "original_copy2_chr": row.get("original_copy2_chr", ""),
            "action": row.get("action", ""),
        }
        pairs.append(pair)
        copy_lookup[copy1_chr] = "copy1"
        copy_lookup[copy2_chr] = "copy2"
        pair_lookup[copy1_chr] = pair
        pair_lookup[copy2_chr] = pair
    pairs.sort(key=lambda row: natural_key(row.get("pair_id", "")))
    return pairs, copy_lookup, pair_lookup


def read_provenance(path, copy_lookup):
    out = {}
    for row in read_tsv(path):
        chrom = row.get("chr", "")
        if not chrom:
            continue
        row = dict(row)
        row["source_copy"] = copy_lookup.get(chrom, "")
        out[chrom] = row
    return out


def normalize_member_list(values):
    members = []
    for value in values:
        if not value:
            continue
        for item in str(value).replace(";", ",").split(","):
            item = item.strip()
            if item:
                members.append(item)
    return sorted(set(members), key=natural_key)


def best_hits_by_label(evidence_dir):
    evidence_dir = Path(evidence_dir)
    labels = ["pcam", "ppa", "ppe", "par", "pdu", "pav", "pmu", "fve", "rch", "pcom"]
    by_label = {}
    for label in labels:
        path = evidence_dir / ("pzi_vs_%s_best_hits.tsv" % label)
        rows = read_tsv(path)
        if not rows:
            continue
        ref_col = "%s_chr" % label
        by_label[label] = {row.get("pzi_chr", ""): row for row in rows if row.get("pzi_chr", "")}
        for row in rows:
            row["_ref_chr"] = row.get(ref_col, "")
    return by_label


def alignment_summaries_by_label(evidence_dir):
    evidence_dir = Path(evidence_dir)
    labels = ["pcam", "ppa", "ppe", "par", "pdu", "pav", "pmu", "fve", "rch", "pcom"]
    by_label = {}
    for label in labels:
        path = evidence_dir / ("pzi_vs_%s_alignment_summary.tsv" % label)
        rows = read_tsv(path)
        if not rows:
            continue
        by_query = defaultdict(list)
        for row in rows:
            query = row.get("query_chr", "")
            ref = row.get("reference_chr", "")
            bp = parse_int(row.get("aligned_block_bp", "0"))
            covered = parse_int(row.get("aligned_covered_bp", "0"))
            if query and ref and bp > 0:
                by_query[query].append({
                    "ref_chr": ref,
                    "aligned_bp": bp,
                    "covered_bp": covered,
                    "matching_bp": parse_int(row.get("matching_bp", "0")),
                })
        for query in by_query:
            by_query[query].sort(key=lambda item: (-item["aligned_bp"], natural_key(item["ref_chr"])))
        by_label[label] = by_query
    return by_label


def read_self_synteny_matrix(evidence_dir):
    matrix = {}
    for row in read_tsv(Path(evidence_dir) / "chr32_all_vs_all_matrix.tsv"):
        query = row.get("query_chr", "")
        target = row.get("target_chr", "")
        if query and target:
            matrix[(query, target)] = parse_int(row.get("aligned_bp", "0"))
    return matrix


def chromosome_lengths_from_provenance(provenance):
    lengths = {}
    for chrom, row in provenance.items():
        if chrom.startswith("chr"):
            lengths[chrom] = parse_int(row.get("length", "0"))
    return lengths


def expected_chromosomes_from_provenance(provenance):
    chroms = sorted([chrom for chrom in provenance if chrom.startswith("chr")], key=natural_key)
    expected = ["chr%02d" % i for i in range(1, 33)]
    if chroms != expected:
        missing = sorted(set(expected) - set(chroms), key=natural_key)
        extra = sorted(set(chroms) - set(expected), key=natural_key)
        raise SystemExit("ERROR: x=8 final phasing expects chr01-32 in provenance. missing=%s extra=%s" % (
            ",".join(missing), ",".join(extra)
        ))
    return chroms


def ref_chr(row, label):
    if not row:
        return ""
    return row.get("%s_chr" % label, "") or row.get("_ref_chr", "")


def top_refs_for_members(best_hits, label, members, min_cov=0.02):
    refs = []
    for chrom in members:
        row = best_hits.get(label, {}).get(chrom, {})
        ref = ref_chr(row, label)
        cov = parse_float(row.get("coverage_ratio", "0"))
        if ref and cov >= min_cov:
            refs.append(ref)
    return refs


def ref_coherence_score(best_hits, labels, members, min_cov=0.02):
    scores = []
    details = []
    for label in labels:
        refs = top_refs_for_members(best_hits, label, members, min_cov=min_cov)
        if not refs:
            continue
        counts = Counter(refs)
        top_count = counts.most_common(1)[0][1]
        if top_count == 4:
            score = 1.0
        elif len(counts) == 2 and sorted(counts.values()) == [2, 2]:
            score = 0.72
        elif top_count == 3:
            score = 0.62
        elif top_count == 2:
            score = 0.48
        else:
            score = 0.30
        scores.append(score)
        detail = "%s:%s" % (label, ",".join("%s=%s" % (k, counts[k]) for k in sorted(counts, key=natural_key)))
        details.append(detail)
    if not scores:
        return 0.40, "no_reference_best_hit_data"
    return sum(scores) / len(scores), ";".join(details)


def ppa_pair_score(best_hits, members):
    refs = top_refs_for_members(best_hits, "ppa", members, min_cov=0.02)
    if not refs:
        return 0.40, "ppa:no_data"
    counts = Counter(refs)
    if len(counts) == 2 and sorted(counts.values()) == [2, 2]:
        return 1.00, "ppa:" + ",".join("%s=%s" % (k, counts[k]) for k in sorted(counts, key=natural_key))
    if len(counts) == 1 and list(counts.values())[0] == 4:
        return 0.70, "ppa:single_reference_group"
    if len(counts) >= 3:
        return 0.55, "ppa:fragmented_" + ",".join("%s=%s" % (k, counts[k]) for k in sorted(counts, key=natural_key))
    return 0.50, "ppa:weak"


def n8_score(best_hits, members):
    return ref_coherence_score(best_hits, ["pcam", "ppe", "par", "pdu", "pav", "pmu"], members)


def outgroup_score(best_hits, members):
    return ref_coherence_score(best_hits, ["fve", "rch", "pcom"], members, min_cov=0.004)


def reference_profile_similarity(profile_a, profile_b):
    total_a = sum(item["aligned_bp"] for item in profile_a)
    total_b = sum(item["aligned_bp"] for item in profile_b)
    if total_a <= 0 or total_b <= 0:
        return None
    by_a = defaultdict(int)
    by_b = defaultdict(int)
    for item in profile_a:
        by_a[item["ref_chr"]] += item["aligned_bp"]
    for item in profile_b:
        by_b[item["ref_chr"]] += item["aligned_bp"]
    shared = 0
    for ref in set(by_a).intersection(by_b):
        shared += min(by_a[ref], by_b[ref])
    denom = min(total_a, total_b)
    return shared / denom if denom else None


def segment_synteny_score(summaries, members):
    label_scores = []
    detail_parts = []
    for label in ["pcam", "ppe", "par", "pdu", "pav", "pmu", "ppa", "fve", "rch", "pcom"]:
        by_query = summaries.get(label, {})
        sims = []
        for a, b in itertools.combinations(members, 2):
            sim = reference_profile_similarity(by_query.get(a, []), by_query.get(b, []))
            if sim is not None:
                sims.append(sim)
        if not sims:
            continue
        score = sum(sims) / len(sims)
        if label in {"fve", "rch", "pcom"}:
            score = min(1.0, score * 1.5)
        label_scores.append(score)
        detail_parts.append("%s=%.3f" % (label, score))
    if not label_scores:
        return 0.40, "segment:no_data"
    return sum(label_scores) / len(label_scores), "segment:" + ",".join(detail_parts)


def self_synteny_score(self_matrix, lengths, members):
    if not self_matrix:
        return 0.40, "self:no_data"
    vals = []
    for a, b in itertools.combinations(members, 2):
        bp = max(self_matrix.get((a, b), 0), self_matrix.get((b, a), 0))
        denom = min(lengths.get(a, 0), lengths.get(b, 0))
        if denom > 0:
            vals.append(min(1.0, bp / denom))
    if not vals:
        return 0.40, "self:no_data"
    avg = sum(vals) / len(vals)
    return min(1.0, avg * 2.0), "self_pairwise=%.4f" % avg


def length_balance_score(lengths, members):
    vals = [lengths.get(m, 0) for m in members if lengths.get(m, 0) > 0]
    if len(vals) != 4:
        return 0.50, "length:no_data"
    mean = sum(vals) / len(vals)
    if mean <= 0:
        return 0.50, "length:no_data"
    variance = sum((v - mean) ** 2 for v in vals) / len(vals)
    cv = math.sqrt(variance) / mean
    score = max(0.25, 1.0 - min(1.0, cv))
    return score, "length_cv=%.4f" % cv


def pairs_inside_quartet(pairs, members):
    member_set = set(members)
    inside = []
    for pair in pairs:
        if pair["copy1_chr"] in member_set and pair["copy2_chr"] in member_set:
            inside.append(pair)
    inside.sort(key=lambda row: natural_key(row.get("pair_id", "")))
    return inside


def assign_hap_slots(quartet, pairs, provenance):
    members = quartet["members"]
    inside = pairs_inside_quartet(pairs, members)
    assignment = {}
    pair_for_hap = {}
    if len(inside) >= 2:
        assignment["hapA1"] = inside[0]["copy1_chr"]
        assignment["hapA2"] = inside[0]["copy2_chr"]
        assignment["hapB1"] = inside[1]["copy1_chr"]
        assignment["hapB2"] = inside[1]["copy2_chr"]
        pair_for_hap["hapA1"] = inside[0]
        pair_for_hap["hapA2"] = inside[0]
        pair_for_hap["hapB1"] = inside[1]
        pair_for_hap["hapB2"] = inside[1]
        basis = "two_corrected_copy_pairs"
    else:
        copy1_members = [m for m in members if provenance.get(m, {}).get("source_copy") == "copy1"]
        copy2_members = [m for m in members if provenance.get(m, {}).get("source_copy") == "copy2"]
        copy1_members.sort(key=natural_key)
        copy2_members.sort(key=natural_key)
        assignment["hapA1"] = copy1_members[0] if len(copy1_members) > 0 else ""
        assignment["hapB1"] = copy1_members[1] if len(copy1_members) > 1 else ""
        assignment["hapA2"] = copy2_members[0] if len(copy2_members) > 0 else ""
        assignment["hapB2"] = copy2_members[1] if len(copy2_members) > 1 else ""
        basis = "copy_balance_fallback"
    return assignment, pair_for_hap, inside, basis


def quartet_constraint_score(quartet):
    basis = quartet["basis"]
    return 0.60


def score_quartet(quartet, pairs, provenance, best_hits, summaries=None, self_matrix=None, lengths=None):
    members = quartet["members"]
    summaries = summaries or {}
    self_matrix = self_matrix or {}
    lengths = lengths or {}
    inside = pairs_inside_quartet(pairs, members)
    copies = [provenance.get(m, {}).get("source_copy", "") for m in members]
    copy_counts = Counter(copies)
    copy_balance = 1.0 if copy_counts.get("copy1", 0) == 2 and copy_counts.get("copy2", 0) == 2 else 0.25
    pair_score = min(1.0, len(inside) / 2.0)
    ppa_score, ppa_detail = ppa_pair_score(best_hits, members)
    n8, n8_detail = n8_score(best_hits, members)
    outgroup, outgroup_detail = outgroup_score(best_hits, members)
    segment, segment_detail = segment_synteny_score(summaries, members)
    self_score, self_detail = self_synteny_score(self_matrix, lengths, members)
    length_score, length_detail = length_balance_score(lengths, members)
    score = (
        0.16 * copy_balance
        + 0.14 * pair_score
        + 0.14 * ppa_score
        + 0.20 * n8
        + 0.10 * outgroup
        + 0.16 * segment
        + 0.06 * self_score
        + 0.04 * length_score
    )
    review_flag = "pass" if score >= 0.72 else "review_required"
    evidence = [
        "copy_balance=%.2f" % copy_balance,
        "corrected_pair_support=%.2f" % pair_score,
        ppa_detail,
        "n8_refs=%.2f[%s]" % (n8, n8_detail),
        "outgroups=%.2f[%s]" % (outgroup, outgroup_detail),
        "%.2f[%s]" % (segment, segment_detail),
        "%.2f[%s]" % (self_score, self_detail),
        "%.2f[%s]" % (length_score, length_detail),
        "hard_constraint=x=8_exact_cover_only",
    ]
    return {
        "copy_balance_score": copy_balance,
        "corrected_pair_score": pair_score,
        "ppa_pair_score": ppa_score,
        "n8_reference_score": n8,
        "outgroup_score": outgroup,
        "segment_synteny_score": segment,
        "self_synteny_score": self_score,
        "length_balance_score": length_score,
        "assignment_score": score,
        "assignment_probability_like": probability_like(score),
        "review_flag": review_flag,
        "confidence": confidence_label(score, review_flag),
        "evidence_summary": "; ".join(evidence),
    }


def partition_copy_pairs(pair_rows):
    seen = set()
    partitions = []
    for combo in itertools.combinations(pair_rows, 2):
        group1_pairs = list(combo)
        group1_ids = frozenset(pair["pair_id"] for pair in group1_pairs)
        group2_pairs = [pair for pair in pair_rows if pair["pair_id"] not in group1_ids]
        group2_ids = frozenset(pair["pair_id"] for pair in group2_pairs)
        if len(group2_pairs) != 2:
            continue
        key = tuple(sorted([tuple(sorted(group1_ids)), tuple(sorted(group2_ids))]))
        if key in seen:
            continue
        seen.add(key)
        partitions.append((group1_pairs, group2_pairs))
    return partitions


def members_from_pairs(pair_rows):
    members = []
    for pair in pair_rows:
        members.extend([pair["copy1_chr"], pair["copy2_chr"]])
    return sorted(members, key=natural_key)


def build_quartet_candidate_scores(chromosomes, pairs, provenance, best_hits, summaries, self_matrix, lengths):
    by_members = {}
    rows = []
    for combo in itertools.combinations(chromosomes, 4):
        members = sorted(combo, key=natural_key)
        pseudo = {
            "quartet_id": "candidate_quartet",
            "x_chr": "",
            "members": members,
            "basis": "all_chr_exact_cover_candidate",
            "review_note": "",
        }
        score_data = score_quartet(pseudo, pairs, provenance, best_hits, summaries, self_matrix, lengths)
        inside = pairs_inside_quartet(pairs, members)
        key = frozenset(members)
        by_members[key] = {
            "members": members,
            "score": score_data["assignment_score"],
            "probability_like": score_data["assignment_probability_like"],
            "evidence": score_data["evidence_summary"],
        }
        rows.append({
            "members": ",".join(members),
            "candidate_score": "%.4f" % score_data["assignment_score"],
            "candidate_probability_like": "%.4f" % score_data["assignment_probability_like"],
            "copy_balance_score": "%.4f" % score_data["copy_balance_score"],
            "corrected_pair_score": "%.4f" % score_data["corrected_pair_score"],
            "ppa_pair_score": "%.4f" % score_data["ppa_pair_score"],
            "n8_reference_score": "%.4f" % score_data["n8_reference_score"],
            "outgroup_score": "%.4f" % score_data["outgroup_score"],
            "segment_synteny_score": "%.4f" % score_data["segment_synteny_score"],
            "self_synteny_score": "%.4f" % score_data["self_synteny_score"],
            "length_balance_score": "%.4f" % score_data["length_balance_score"],
            "corrected_pair_ids_inside": ",".join(pair["pair_id"] for pair in inside),
            "evidence_summary": score_data["evidence_summary"],
        })
    rows.sort(key=lambda row: (-parse_float(row["candidate_score"]), natural_key(row["members"])))
    return by_members, rows


def enumerate_exact_cover_matchings(chromosomes, candidate_scores, keep_solutions=32, branch_per_first=220, beam_width=8000):
    chrom_count = len(chromosomes)
    if chrom_count != 32:
        raise SystemExit("ERROR: x=8 final phasing expects 32 chromosomes, observed %s" % chrom_count)
    index = {chrom: idx for idx, chrom in enumerate(chromosomes)}
    full_mask = (1 << chrom_count) - 1
    candidates_by_first = defaultdict(list)
    for key, data in candidate_scores.items():
        idxs = sorted(index[m] for m in key)
        mask = 0
        for idx in idxs:
            mask |= 1 << idx
        item = {
            "mask": mask,
            "members": data["members"],
            "score": data["score"],
        }
        for idx in idxs:
            candidates_by_first[idx].append(item)
    for first in candidates_by_first:
        candidates_by_first[first].sort(key=lambda item: (-item["score"], natural_key(",".join(item["members"]))))
        candidates_by_first[first] = candidates_by_first[first][:branch_per_first]

    states = [(0.0, full_mask, [])]
    for _depth in range(8):
        next_by_mask = {}
        for score, mask, groups in states:
            if mask == 0:
                key = (mask, " | ".join(",".join(group) for group in groups))
                next_by_mask[key] = (score, mask, groups)
                continue
            first = None
            for idx in range(chrom_count):
                if mask & (1 << idx):
                    first = idx
                    break
            for candidate in candidates_by_first.get(first, []):
                candidate_mask = candidate["mask"]
                if candidate_mask & mask != candidate_mask:
                    continue
                next_mask = mask & ~candidate_mask
                next_groups = groups + [candidate["members"]]
                next_score = score + candidate["score"]
                old = next_by_mask.get(next_mask)
                if old is None or next_score > old[0]:
                    next_by_mask[next_mask] = (next_score, next_mask, next_groups)
        states = sorted(
            next_by_mask.values(),
            key=lambda item: (-item[0], natural_key(" | ".join(",".join(group) for group in item[2]))),
        )[:beam_width]
        if not states:
            break

    covers = [(score, groups) for score, mask, groups in states if mask == 0 and len(groups) == 8]
    covers.sort(key=lambda item: (-item[0], natural_key(" | ".join(",".join(group) for group in item[1]))))
    covers = covers[:keep_solutions]
    if not covers:
        raise SystemExit("ERROR: no x=8 exact-cover haplotype solution generated by the beam search")
    return covers


def build_final_quartets(evidence_dir, pairs, provenance, best_hits, summaries, self_matrix, lengths,
                         keep_solutions=32, branch_per_first=220, beam_width=8000):
    chromosomes = expected_chromosomes_from_provenance(provenance)
    candidate_scores, candidate_rows = build_quartet_candidate_scores(
        chromosomes, pairs, provenance, best_hits, summaries, self_matrix, lengths
    )
    matchings = enumerate_exact_cover_matchings(
        chromosomes,
        candidate_scores,
        keep_solutions=keep_solutions,
        branch_per_first=branch_per_first,
        beam_width=beam_width,
    )
    top_score = matchings[0][0]
    exp_scores = [math.exp(5.0 * (score - top_score)) for score, _groups in matchings]
    exp_total = sum(exp_scores) if exp_scores else 1.0
    matching_rows = []
    for rank, (score, groups) in enumerate(matchings, start=1):
        next_score = matchings[rank][0] if rank < len(matchings) else 0.0
        quartet_desc = [",".join(group) for group in groups]
        matching_rows.append({
            "rank": str(rank),
            "total_score": "%.4f" % score,
            "mean_quartet_score": "%.4f" % (score / 8.0),
            "relative_probability_like": "%.4f" % (score / top_score if top_score else 0.0),
            "global_probability_like": "%.4f" % (exp_scores[rank - 1] / exp_total if exp_total else 0.0),
            "score_margin_to_next": "%.4f" % (score - next_score if rank < len(matchings) else score),
            "chosen": "yes" if rank == 1 else "no",
            "pair_groups": "",
            "quartets": " | ".join(quartet_desc),
        })
    selected_groups = matchings[0][1]
    quartets = []
    for idx, members in enumerate(sorted(selected_groups, key=lambda group: natural_key(",".join(group))), start=1):
        quartets.append({
            "source_id": "global_x8_exact_cover_%02d" % idx,
            "members": members,
            "basis": "global_x8_exact_cover",
            "review_note": "Selected by global x=8 exact-cover search over all chr01-32 four-chromosome candidates.",
            "consensus_key": "",
            "ppa_groups": "",
            "quartet_id": "final_quartet%02d" % idx,
            "x_chr": "chr%02d" % idx,
        })
    return quartets, candidate_rows, matching_rows


def build_assignments(final_quartets, pairs, provenance, best_hits, summaries, self_matrix, lengths):
    quartet_rows = []
    chrom_rows = []
    rename_map = {}
    review_rows = []
    seen = set()
    for quartet in final_quartets:
        members = quartet["members"]
        if len(set(members)) != 4:
            raise SystemExit("ERROR: quartet has duplicate members: %s" % quartet["quartet_id"])
        if seen.intersection(members):
            raise SystemExit("ERROR: chromosome appears in multiple final quartets: %s" % ",".join(sorted(seen.intersection(members))))
        seen.update(members)
        assignment, pair_for_hap, inside_pairs, hap_basis = assign_hap_slots(quartet, pairs, provenance)
        missing_haps = [hap for hap in HAP_ORDER if not assignment.get(hap, "")]
        if missing_haps:
            raise SystemExit(
                "ERROR: quartet %s cannot be assigned to all hapA1/hapA2/hapB1/hapB2 slots. "
                "missing=%s members=%s" % (
                    quartet["quartet_id"],
                    ",".join(missing_haps),
                    ",".join(sorted(members, key=natural_key)),
                )
            )
        scoring = score_quartet(quartet, pairs, provenance, best_hits, summaries, self_matrix, lengths)
        row = {
            "quartet_id": quartet["quartet_id"],
            "x_chr": quartet["x_chr"],
            "members": ",".join(sorted(members, key=natural_key)),
            "hapA1_chr": assignment.get("hapA1", ""),
            "hapA2_chr": assignment.get("hapA2", ""),
            "hapB1_chr": assignment.get("hapB1", ""),
            "hapB2_chr": assignment.get("hapB2", ""),
            "corrected_pair_ids": ",".join(pair["pair_id"] for pair in inside_pairs),
            "assignment_basis": quartet["basis"] + "+" + hap_basis,
            "assignment_score": "%.4f" % scoring["assignment_score"],
            "assignment_probability_like": "%.4f" % scoring["assignment_probability_like"],
            "confidence": scoring["confidence"],
            "review_flag": scoring["review_flag"],
            "evidence_summary": scoring["evidence_summary"],
            "review_note": quartet["review_note"],
        }
        quartet_rows.append(row)
        if scoring["review_flag"] != "pass":
            review_rows.append({
                "item": quartet["quartet_id"],
                "status": scoring["review_flag"],
                "score": "%.4f" % scoring["assignment_score"],
                "probability_like": "%.4f" % scoring["assignment_probability_like"],
                "evidence": scoring["evidence_summary"],
                "recommendation": quartet["review_note"],
            })
        for hap in HAP_ORDER:
            old_chr = assignment.get(hap, "")
            if not old_chr:
                continue
            suffix = hap.replace("hap", "")
            new_chr = "%s_%s" % (quartet["x_chr"], suffix)
            prov = provenance.get(old_chr, {})
            pair = pair_for_hap.get(hap, {})
            chrom_row = {
                "old_chr": old_chr,
                "new_chr": new_chr,
                "hap_set": hap,
                "quartet_id": quartet["quartet_id"],
                "x_chr": quartet["x_chr"],
                "source_copy": prov.get("source_copy", ""),
                "dominant_hap": prov.get("dominant_hap", ""),
                "dominance": prov.get("dominance", ""),
                "copy_pair_id": pair.get("pair_id", ""),
                "assignment_score": row["assignment_score"],
                "assignment_probability_like": row["assignment_probability_like"],
                "confidence": row["confidence"],
                "review_flag": row["review_flag"],
                "evidence": row["assignment_basis"],
            }
            chrom_rows.append(chrom_row)
            rename_map[old_chr] = chrom_row
    expected = {"chr%02d" % i for i in range(1, 33)}
    missing = sorted(expected - set(rename_map), key=natural_key)
    extra = sorted(set(rename_map) - expected, key=natural_key)
    if missing or extra:
        raise SystemExit("ERROR: final assignment does not cover chr01-32. missing=%s extra=%s" % (",".join(missing), ",".join(extra)))
    return quartet_rows, chrom_rows, rename_map, review_rows


def stream_renamed_fastas(pzi_fasta, out_dir, rename_map):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    hap_paths = {
        "hapA1": out_dir / "pzi_hapA1.fa",
        "hapA2": out_dir / "pzi_hapA2.fa",
        "hapB1": out_dir / "pzi_hapB1.fa",
        "hapB2": out_dir / "pzi_hapB2.fa",
    }
    all_path = out_dir / "pzi_A1A2B1B2.renamed.chr32.fa"
    hap_handles = {hap: path.open("w", encoding="utf-8") for hap, path in hap_paths.items()}
    all_handle = all_path.open("w", encoding="utf-8")
    counts = defaultdict(int)
    bases = defaultdict(int)
    found = set()
    current = None
    current_row = None

    def write_header(row):
        header = ">%s old_chr=%s hap_set=%s quartet_id=%s x_chr=%s\n" % (
            row["new_chr"], row["old_chr"], row["hap_set"], row["quartet_id"], row["x_chr"]
        )
        all_handle.write(header)
        hap_handles[row["hap_set"]].write(header)

    try:
        with Path(pzi_fasta).open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith(">"):
                    old = line[1:].strip().split()[0]
                    current = old
                    current_row = rename_map.get(old)
                    if current_row:
                        found.add(old)
                        write_header(current_row)
                    continue
                if not current_row:
                    continue
                seq = line.strip()
                if not seq:
                    continue
                all_handle.write(seq + "\n")
                hap_handles[current_row["hap_set"]].write(seq + "\n")
                counts[current_row["hap_set"]] += 0
                bases[current_row["hap_set"]] += len(seq)
                bases["all"] += len(seq)
        for row in rename_map.values():
            counts[row["hap_set"]] += 1
        missing = sorted(set(rename_map) - found, key=natural_key)
        if missing:
            raise SystemExit("ERROR: assigned chromosomes missing from FASTA: %s" % ",".join(missing))
    finally:
        all_handle.close()
        for handle in hap_handles.values():
            handle.close()
    return hap_paths, all_path, counts, bases


def write_sizes_and_checks(out_dir, hap_paths, all_path):
    rows = []
    for label, path in list(hap_paths.items()) + [("all_renamed_chr32", all_path)]:
        seq_count = 0
        total_bp = 0
        current = None
        current_bp = 0
        sizes_path = Path(str(path) + ".sizes")
        with Path(path).open("r", encoding="utf-8", errors="replace") as handle, sizes_path.open("w", encoding="utf-8") as sizes:
            for line in handle:
                if line.startswith(">"):
                    if current is not None:
                        sizes.write("%s\t%s\n" % (current, current_bp))
                        total_bp += current_bp
                    current = line[1:].strip().split()[0]
                    current_bp = 0
                    seq_count += 1
                else:
                    current_bp += len(line.strip())
            if current is not None:
                sizes.write("%s\t%s\n" % (current, current_bp))
                total_bp += current_bp
        rows.append({
            "file_label": label,
            "path": str(path),
            "seq_count": str(seq_count),
            "total_bp": str(total_bp),
            "sizes_path": str(sizes_path),
        })
    return rows


def write_annotation_next_steps(out_dir):
    out_dir = Path(out_dir)
    rows = []
    for hap in HAP_ORDER:
        rows.append({
            "hap_set": hap,
            "renamed_fasta": str(out_dir / ("pzi_%s.fa" % hap)),
            "suggested_edta_dir": "annotation/repeats/%s_edta" % hap,
            "suggested_helixer_dir": "annotation/genes/%s_helixer_edta_masked" % hap,
            "expected_gff_after_postprocess": "annotation/genes/%s_helixer_edta_masked/postprocess/pzi_%s.edta_masked.helixer.genes.gff3" % (hap, hap),
            "expected_cds_after_postprocess": "annotation/genes/%s_helixer_edta_masked/postprocess/pzi_%s.edta_masked.helixer.cds.fa" % (hap, hap),
            "expected_pep_after_postprocess": "annotation/genes/%s_helixer_edta_masked/postprocess/pzi_%s.edta_masked.helixer.protein.fa" % (hap, hap),
            "notes": "Run EDTA/Helixer/postprocess after FASTA extraction; use jobs/annotation_03 and jobs/annotation_08 as templates.",
        })
    write_tsv(out_dir / "annotation_next_steps.tsv", rows, [
        "hap_set", "renamed_fasta", "suggested_edta_dir", "suggested_helixer_dir",
        "expected_gff_after_postprocess", "expected_cds_after_postprocess",
        "expected_pep_after_postprocess", "notes",
    ])


def write_readme(out_dir):
    text = """Final A1/A2/B1/B2 haplotype assignment
========================================

This directory is the final, evidence-scored haplotype split for the
P. zippeliana 2n=32 chromosome assembly.

Core outputs:
- FINAL_HAPLOTYPE_OVERVIEW.tsv: main 8-row quartet table with hapA1/A2/B1/B2 slots, score, probability-like value, and key evidence.
- FINAL_REVIEW_NOTES.tsv: compact review table for global ambiguity, low-confidence quartets, or validation warnings.
- FINAL_FILES.tsv: manifest for final FASTA and key tables.
- final_assignment_summary.tsv: validation counts; each hap set must have 8 chromosomes.
- pzi_A1A2B1B2.rename_map.tsv: 32-row chromosome rename map.
- pzi_hapA1.fa, pzi_hapA2.fa, pzi_hapB1.fa, pzi_hapB2.fa: renamed FASTA files.
- pzi_A1A2B1B2.renamed.chr32.fa: all 32 renamed chromosomes.
- details/: detailed candidate scores, global matching alternatives, and full evidence tables.

Renaming convention:
- chr01_A1, chr01_A2, chr01_B1, chr01_B2 belong to final quartet chr01.
- Each hap FASTA contains exactly 8 chromosomes.

Important interpretation:
- copy1/copy2 are not treated as subgenomes. They are phase-copy groups
  approximately corresponding to A1+B1 and A2+B2.
- The only fixed biological constraint is x=8. The script scores all legal
  four-chromosome candidates from chr01-32 and selects a global exact-cover
  solution with eight quartets. copy1/copy2, corrected pair support,
  multi-reference synteny, outgroups, segment-level profiles, and self-synteny
  are evidence terms only, not hard quartet priors.
- A1/A2/B1/B2 are stable output slots. Without parental data, the labels should
  be interpreted as four phased haplotype sets rather than proven parental
  subgenome identities.

Annotation:
- This step only finalizes and exports FASTA. Run EDTA/Helixer/postprocess on
  the four exported FASTA files to generate final GFF/CDS/PEP.
- annotation_next_steps.tsv lists recommended output locations and templates.
"""
    (Path(out_dir) / "README_FINAL_HAPLOTYPES.md").write_text(text, encoding="utf-8")


def validate_and_summarize(quartet_rows, chrom_rows, fasta_rows, matching_rows=None):
    summary = []
    summary.append({"metric": "final_quartets", "value": str(len(quartet_rows)), "status": "pass" if len(quartet_rows) == 8 else "fail"})
    summary.append({"metric": "assigned_chromosomes", "value": str(len(chrom_rows)), "status": "pass" if len(chrom_rows) == 32 else "fail"})
    hap_counts = Counter(row["hap_set"] for row in chrom_rows)
    for hap in HAP_ORDER:
        summary.append({"metric": "%s_chromosomes" % hap, "value": str(hap_counts.get(hap, 0)), "status": "pass" if hap_counts.get(hap, 0) == 8 else "fail"})
    copy_counts = Counter(row["source_copy"] for row in chrom_rows)
    summary.append({"metric": "copy1_assigned", "value": str(copy_counts.get("copy1", 0)), "status": "pass" if copy_counts.get("copy1", 0) == 16 else "fail"})
    summary.append({"metric": "copy2_assigned", "value": str(copy_counts.get("copy2", 0)), "status": "pass" if copy_counts.get("copy2", 0) == 16 else "fail"})
    for row in fasta_rows:
        expected = "32" if row["file_label"] == "all_renamed_chr32" else "8"
        summary.append({
            "metric": "%s_fasta_seq_count" % row["file_label"],
            "value": row["seq_count"],
            "status": "pass" if row["seq_count"] == expected else "fail",
        })
    if matching_rows:
        top = matching_rows[0]
        margin = parse_float(top.get("score_margin_to_next", "0"))
        summary.append({
            "metric": "global_matching_top1_score",
            "value": top.get("total_score", ""),
            "status": "pass",
        })
        summary.append({
            "metric": "global_matching_top1_probability_like",
            "value": top.get("global_probability_like", ""),
            "status": "pass",
        })
        summary.append({
            "metric": "global_matching_top1_top2_margin",
            "value": top.get("score_margin_to_next", ""),
            "status": "review" if margin < 0.05 else "pass",
        })
    return summary


def build_overview_rows(quartet_rows):
    rows = []
    for row in sorted(quartet_rows, key=lambda r: natural_key(r["x_chr"])):
        rows.append({
            "x_chr": row["x_chr"],
            "hapA1": row["hapA1_chr"],
            "hapA2": row["hapA2_chr"],
            "hapB1": row["hapB1_chr"],
            "hapB2": row["hapB2_chr"],
            "members": row["members"],
            "score": row["assignment_score"],
            "probability_like": row["assignment_probability_like"],
            "confidence": row["confidence"],
            "review_flag": row["review_flag"],
            "key_evidence": short_evidence_summary(row["evidence_summary"]),
        })
    return rows


def build_files_manifest(out_dir, hap_paths, all_path, fasta_rows):
    rows = []
    size_by_path = {row["path"]: row for row in fasta_rows}
    for hap in HAP_ORDER:
        path = hap_paths[hap]
        stat = size_by_path.get(str(path), {})
        rows.append({
            "file_role": "%s_genome_fasta" % hap,
            "path": str(path),
            "seq_count": stat.get("seq_count", ""),
            "total_bp": stat.get("total_bp", ""),
            "description": "Renamed x=8 FASTA for %s." % hap,
        })
    stat = size_by_path.get(str(all_path), {})
    rows.append({
        "file_role": "all_renamed_chr32_fasta",
        "path": str(all_path),
        "seq_count": stat.get("seq_count", ""),
        "total_bp": stat.get("total_bp", ""),
        "description": "All 32 chromosomes renamed as chrXX_A1/A2/B1/B2.",
    })
    for name, desc in [
        ("FINAL_HAPLOTYPE_OVERVIEW.tsv", "Main human-readable quartet and haplotype assignment table."),
        ("FINAL_REVIEW_NOTES.tsv", "Global and quartet-level review notes."),
        ("FINAL_FILES.tsv", "Final FASTA and key table manifest."),
        ("details/final_quartet_assignment.tsv", "Full quartet table with detailed evidence."),
        ("details/final_chromosome_assignment.tsv", "Full chromosome rename map."),
        ("details/candidate_quartet_scores.tsv", "All 35,960 four-chromosome candidate scores."),
        ("details/global_matching_solutions.tsv", "Top global x=8 exact-cover solutions."),
    ]:
        rows.append({
            "file_role": name.replace("/", "_"),
            "path": str(Path(out_dir) / name),
            "seq_count": "",
            "total_bp": "",
            "description": desc,
        })
    return rows


def build_review_notes(quartet_rows, matching_rows, summary_rows):
    rows = []
    if matching_rows:
        top = matching_rows[0]
        margin = parse_float(top.get("score_margin_to_next", "0"))
        status = "pass" if margin >= 0.05 else "review"
        rows.append({
            "level": "global",
            "item": "chosen_x8_exact_cover",
            "status": status,
            "score_or_value": top.get("total_score", ""),
            "probability_like": top.get("global_probability_like", ""),
            "evidence": "top1_top2_margin=%s; quartets=%s" % (top.get("score_margin_to_next", ""), top.get("quartets", "")),
            "recommendation": "If margin is small, inspect FINAL_HAPLOTYPE_OVERVIEW.tsv together with dotplots/heatmaps before treating A/B labels as biological subgenomes.",
        })
        for alt in matching_rows[1:4]:
            rows.append({
                "level": "global_alternative",
                "item": "alternative_rank_%s" % alt.get("rank", ""),
                "status": "alternative",
                "score_or_value": alt.get("total_score", ""),
                "probability_like": alt.get("global_probability_like", ""),
                "evidence": "margin_to_next=%s; quartets=%s" % (alt.get("score_margin_to_next", ""), alt.get("quartets", "")),
                "recommendation": "Alternative exact-cover solution retained for review only.",
            })
    for row in quartet_rows:
        if row.get("review_flag") != "pass":
            rows.append({
                "level": "quartet",
                "item": row.get("x_chr", row.get("quartet_id", "")),
                "status": row.get("review_flag", ""),
                "score_or_value": row.get("assignment_score", ""),
                "probability_like": row.get("assignment_probability_like", ""),
                "evidence": short_evidence_summary(row.get("evidence_summary", "")),
                "recommendation": "Inspect corresponding multi-reference dotplots/heatmaps and candidate_quartet_scores.tsv.",
            })
    for row in summary_rows:
        if row.get("status") != "pass":
            rows.append({
                "level": "validation",
                "item": row.get("metric", ""),
                "status": row.get("status", ""),
                "score_or_value": row.get("value", ""),
                "probability_like": "",
                "evidence": "validation metric did not pass",
                "recommendation": "Do not use final haplotypes until this validation issue is resolved.",
            })
    if not rows:
        rows.append({
            "level": "global",
            "item": "all_checks",
            "status": "pass",
            "score_or_value": "",
            "probability_like": "",
            "evidence": "No review-required quartet or validation issue was detected.",
            "recommendation": "Use FINAL_HAPLOTYPE_OVERVIEW.tsv and FASTA files as final phased haplotype outputs.",
        })
    return rows


def parse_args():
    parser = argparse.ArgumentParser(description="Finalize PZI A1/A2/B1/B2 haplotype assignment and export renamed FASTA.")
    parser.add_argument("--root-dir", required=True)
    parser.add_argument("--evidence-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--pzi-chrom-fasta", required=True)
    parser.add_argument("--chr-hap-provenance", required=True)
    parser.add_argument("--corrected-homolog-pairs", required=True)
    parser.add_argument("--keep-solutions", type=int, default=32)
    parser.add_argument("--branch-per-first", type=int, default=220)
    parser.add_argument("--beam-width", type=int, default=8000)
    return parser.parse_args()


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pairs, copy_lookup, pair_lookup = read_corrected_pairs(args.corrected_homolog_pairs)
    provenance = read_provenance(args.chr_hap_provenance, copy_lookup)
    best_hits = best_hits_by_label(Path(args.evidence_dir))
    summaries = alignment_summaries_by_label(Path(args.evidence_dir))
    self_matrix = read_self_synteny_matrix(Path(args.evidence_dir))
    lengths = chromosome_lengths_from_provenance(provenance)

    final_quartets, candidate_rows, matching_rows = build_final_quartets(
        Path(args.evidence_dir), pairs, provenance, best_hits, summaries, self_matrix, lengths,
        keep_solutions=args.keep_solutions,
        branch_per_first=args.branch_per_first,
        beam_width=args.beam_width,
    )
    quartet_rows, chrom_rows, rename_map, review_rows = build_assignments(
        final_quartets, pairs, provenance, best_hits, summaries, self_matrix, lengths
    )
    for row in matching_rows:
        review_rows.append({
            "item": "global_matching_rank_%s" % row["rank"],
            "status": "chosen" if row["chosen"] == "yes" else "alternative",
            "score": row["total_score"],
            "probability_like": row["global_probability_like"],
            "evidence": row["quartets"] + " pair_groups=" + row["pair_groups"] + " margin_to_next=" + row["score_margin_to_next"],
            "recommendation": "The matching with chosen=yes is selected by the highest global x=8 evidence score; small score margins should be treated as ambiguous/review.",
        })

    hap_paths, all_path, counts, bases = stream_renamed_fastas(args.pzi_chrom_fasta, out_dir, rename_map)
    fasta_rows = write_sizes_and_checks(out_dir, hap_paths, all_path)
    summary_rows = validate_and_summarize(quartet_rows, chrom_rows, fasta_rows, matching_rows)
    details_dir = out_dir / "details"
    overview_rows = build_overview_rows(quartet_rows)
    review_note_rows = build_review_notes(quartet_rows, matching_rows, summary_rows)
    file_manifest_rows = build_files_manifest(out_dir, hap_paths, all_path, fasta_rows)

    write_tsv(out_dir / "FINAL_HAPLOTYPE_OVERVIEW.tsv", overview_rows, [
        "x_chr", "hapA1", "hapA2", "hapB1", "hapB2", "members", "score",
        "probability_like", "confidence", "review_flag", "key_evidence",
    ])
    write_tsv(out_dir / "FINAL_REVIEW_NOTES.tsv", review_note_rows, [
        "level", "item", "status", "score_or_value", "probability_like", "evidence", "recommendation",
    ])
    write_tsv(out_dir / "FINAL_FILES.tsv", file_manifest_rows, [
        "file_role", "path", "seq_count", "total_bp", "description",
    ])

    write_tsv(details_dir / "final_quartet_assignment.tsv", quartet_rows, [
        "quartet_id", "x_chr", "members", "hapA1_chr", "hapA2_chr", "hapB1_chr", "hapB2_chr",
        "corrected_pair_ids", "assignment_basis", "assignment_score", "assignment_probability_like",
        "confidence", "review_flag", "evidence_summary", "review_note",
    ])
    write_tsv(details_dir / "final_chromosome_assignment.tsv", sorted(chrom_rows, key=lambda r: (r["hap_set"], natural_key(r["new_chr"]))), [
        "old_chr", "new_chr", "hap_set", "quartet_id", "x_chr", "source_copy", "dominant_hap",
        "dominance", "copy_pair_id", "assignment_score", "assignment_probability_like",
        "confidence", "review_flag", "evidence",
    ])
    write_tsv(out_dir / "pzi_A1A2B1B2.rename_map.tsv", sorted(chrom_rows, key=lambda r: natural_key(r["old_chr"])), [
        "old_chr", "new_chr", "hap_set", "quartet_id", "x_chr", "source_copy", "dominant_hap",
        "dominance", "copy_pair_id", "assignment_score", "assignment_probability_like",
        "confidence", "review_flag", "evidence",
    ])
    write_tsv(details_dir / "pzi_A1A2B1B2.rename_map.tsv", sorted(chrom_rows, key=lambda r: natural_key(r["old_chr"])), [
        "old_chr", "new_chr", "hap_set", "quartet_id", "x_chr", "source_copy", "dominant_hap",
        "dominance", "copy_pair_id", "assignment_score", "assignment_probability_like",
        "confidence", "review_flag", "evidence",
    ])
    candidate_fields = [
        "members", "candidate_score", "candidate_probability_like", "copy_balance_score",
        "corrected_pair_score", "ppa_pair_score", "n8_reference_score", "outgroup_score",
        "segment_synteny_score", "self_synteny_score", "length_balance_score",
        "corrected_pair_ids_inside", "evidence_summary",
    ]
    write_tsv(details_dir / "candidate_quartet_scores.tsv", candidate_rows, candidate_fields)
    write_tsv(details_dir / "candidate_quartet_edge_scores.tsv", candidate_rows, candidate_fields)
    write_tsv(details_dir / "global_matching_solutions.tsv", matching_rows, [
        "rank", "total_score", "mean_quartet_score", "relative_probability_like", "global_probability_like",
        "score_margin_to_next", "chosen", "pair_groups", "quartets",
    ])
    write_tsv(details_dir / "final_assignment_review.tsv", review_rows, [
        "item", "status", "score", "probability_like", "evidence", "recommendation",
    ])
    write_tsv(details_dir / "final_fasta_manifest.tsv", fasta_rows, [
        "file_label", "path", "seq_count", "total_bp", "sizes_path",
    ])
    write_tsv(out_dir / "final_assignment_summary.tsv", summary_rows, ["metric", "value", "status"])
    write_annotation_next_steps(out_dir)
    write_readme(out_dir)

    failures = [row for row in summary_rows if row["status"] == "fail"]
    if failures:
        raise SystemExit("ERROR: final haplotype assignment validation failed: %s" % ",".join(row["metric"] for row in failures))
    print("Final haplotype assignment completed: %s" % out_dir)


if __name__ == "__main__":
    main()
