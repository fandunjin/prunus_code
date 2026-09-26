#!/usr/bin/env python3
"""Normalize official POInT posterior and branch-transition outputs.

The parser never infers a loss from an absent anchor.  It reports only the
posterior track assignments and the branch-state transition probabilities that
POInT itself writes.
"""

from __future__ import print_function

import argparse
import csv
import math
import re
from pathlib import Path


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def entropy(values):
    total = sum(values)
    if total <= 0:
        return 0.0
    return -sum((p / total) * math.log(p / total, 2) for p in values if p > 0)


def parse_posterior(path, outdir, pillar_map_path=None):
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 2:
        raise ValueError("POInT posterior file is empty: {}".format(path))
    header = lines[0].split("\t")
    marker = None
    for index, field in enumerate(header):
        if field.startswith("#"):
            marker = index
            break
    if marker is None:
        raise ValueError("POInT posterior header lacks # separator")
    gene_fields = ["PZI_A_gene", "PZI_B_gene", "PSP_A_gene", "PSP_B_gene", "PPA_A_gene", "PPA_B_gene"]
    pattern_names = []
    first_pattern = header[marker][1:]
    if first_pattern:
        pattern_names.append(first_pattern)
    for field in header[marker + 1:]:
        if field:
            name = field[1:] if field.startswith("#") else field
            if name:
                pattern_names.append(name)
    if not pattern_names:
        raise ValueError("POInT posterior header has no tracking patterns")
    pillar_map = {}
    if pillar_map_path:
        with Path(pillar_map_path).open("r", encoding="utf-8", newline="") as handle:
            for map_row in csv.DictReader(handle, delimiter="\t"):
                pillar_map[map_row["optimized_index"]] = map_row
    rows = []
    for index, line in enumerate(lines[1:]):
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < marker + 1:
            raise ValueError("posterior row {} is shorter than gene columns".format(index))
        genes = fields[:marker]
        probs = []
        # Data rows contain the first probability directly at the marker
        # position; only the header has a literal '#' separator token.
        for token in fields[marker:]:
            if token == "#" or token == "":
                continue
            try:
                probs.append(float(token))
            except ValueError:
                continue
        if len(genes) != 6:
            raise ValueError("posterior row {} has {} gene fields, expected six".format(index, len(genes)))
        if len(probs) != len(pattern_names):
            raise ValueError("posterior row {} has {} probabilities, expected {}".format(index, len(probs), len(pattern_names)))
        best = max(range(len(probs)), key=lambda i: probs[i])
        row = {"pillar_index": str(index), "max_posterior": "{:.12g}".format(probs[best]), "entropy_bits": "{:.12g}".format(entropy(probs))}
        if str(index) in pillar_map:
            row["pillar_id"] = pillar_map[str(index)].get("pillar_id", "")
            row["gtr_global_index"] = pillar_map[str(index)].get("gtr_global_index", "")
            row["gtr_seqid"] = pillar_map[str(index)].get("gtr_seqid", "")
            row["gtr_gene"] = pillar_map[str(index)].get("gtr_gene", "")
        for field, value in zip(gene_fields, genes):
            row[field] = value
        row["best_tracking"] = pattern_names[best]
        rows.append(row)
    map_fields = ["pillar_id", "gtr_global_index", "gtr_seqid", "gtr_gene"] if pillar_map_path else []
    fields = ["pillar_index"] + map_fields + gene_fields + ["best_tracking", "max_posterior", "entropy_bits"]
    write_tsv(Path(outdir) / "point_posterior_summary.tsv", fields, rows)
    mean_max = sum(float(row["max_posterior"]) for row in rows) / float(len(rows))
    with Path(outdir, "point_posterior_qc.tsv").open("w", encoding="ascii") as handle:
        handle.write("field\tvalue\n")
        handle.write("pillars\t{}\n".format(len(rows)))
        handle.write("mean_max_posterior\t{:.12g}\n".format(mean_max))
        handle.write("pattern_count\t{}\n".format(len(pattern_names)))
    return len(rows)


def parse_conditional(path, outdir):
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 3:
        raise ValueError("POInT conditional file is empty: {}".format(path))
    first = lines[0].split("\t")
    second = lines[1].split("\t")
    if not first or first[0] != "Site" or not second or second[0] != "Site":
        raise ValueError("unrecognized POInT conditional header")
    branches = []
    i = 1
    while i + 1 < len(first):
        branches.append((first[i], first[i + 1]))
        i += 2
    all_transition_labels = [field for field in second[1:] if field]
    if not branches or not all_transition_labels:
        raise ValueError("conditional output has no branches or transitions")
    if len(all_transition_labels) % len(branches) != 0:
        raise ValueError("conditional header has {} transition labels for {} branches".format(
            len(all_transition_labels), len(branches)
        ))
    transitions_per_branch = len(all_transition_labels) // len(branches)
    transitions = all_transition_labels[:transitions_per_branch]
    # POInT writes the same state-transition columns once for every branch.
    # Verify this contract before averaging the posterior values by branch.
    for branch_index in range(len(branches)):
        start = branch_index * transitions_per_branch
        end = start + transitions_per_branch
        if all_transition_labels[start:end] != transitions:
            raise ValueError("conditional transition labels differ across POInT branches")
    values = {branch_index: {transition_index: [] for transition_index in range(transitions_per_branch)} for branch_index in range(len(branches))}
    n_sites = 0
    for line in lines[2:]:
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 1 + len(branches) * transitions_per_branch:
            raise ValueError("conditional row has fewer probabilities than header")
        n_sites += 1
        cursor = 1
        for branch_index in range(len(branches)):
            for transition_index in range(transitions_per_branch):
                try:
                    value = float(fields[cursor])
                except ValueError:
                    value = float("nan")
                if not math.isnan(value):
                    values[branch_index][transition_index].append(value)
                cursor += 1
    rows = []
    for branch_index, (branch, branch_length) in enumerate(branches):
        for transition_index, transition in enumerate(transitions):
            vals = values[branch_index][transition_index]
            if not vals:
                continue
            rows.append({
                "branch_index": str(branch_index),
                "branch": branch,
                "branch_length": branch_length,
                "transition": transition,
                "n_pillars": str(len(vals)),
                "mean_probability": "{:.12g}".format(sum(vals) / len(vals)),
                "max_probability": "{:.12g}".format(max(vals)),
            })
    fields = ["branch_index", "branch", "branch_length", "transition", "n_pillars", "mean_probability", "max_probability"]
    write_tsv(Path(outdir) / "point_branch_transition_summary.tsv", fields, rows)
    with Path(outdir, "point_conditional_qc.tsv").open("w", encoding="ascii") as handle:
        handle.write("field\tvalue\n")
        handle.write("pillars\t{}\n".format(n_sites))
        handle.write("branches\t{}\n".format(len(branches)))
        handle.write("transitions\t{}\n".format(len(transitions)))
    return n_sites


def write_conditional_unavailable(outdir, reason):
    """Record that POInT posterior-only mode was used.

    The likelihood and posterior files are still native POInT outputs.  This
    explicit marker prevents a missing conditional matrix from being silently
    interpreted as zero transition probability.
    """
    write_tsv(
        Path(outdir) / "point_branch_transition_summary.tsv",
        ["branch_index", "branch", "branch_length", "transition", "n_pillars", "mean_probability", "max_probability"],
        [],
    )
    with Path(outdir, "point_conditional_qc.tsv").open("w", encoding="utf-8") as handle:
        handle.write("field\tvalue\n")
        handle.write("status\tunavailable\n")
        handle.write("reason\t{}\n".format(reason))
        handle.write("pillars\tNA\n")
        handle.write("branches\tNA\n")
        handle.write("transitions\tNA\n")


def parse_log_likelihood(stdout_path, tree_output_path=None):
    # POInT writes the optimized value in the .out tree description as
    # "Final ln likelihood".  Prefer that exact field over progress lines.
    if tree_output_path and Path(tree_output_path).is_file():
        tree_text = Path(tree_output_path).read_text(encoding="utf-8", errors="replace")
        exact = re.findall(r"Final\s+ln\s+likelihood:\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)", tree_text, flags=re.IGNORECASE)
        if exact:
            return exact[-1]
    text = Path(stdout_path).read_text(encoding="utf-8", errors="replace")
    values = []
    for line in text.splitlines():
        # Progress lines take the form "lnL: -1741.8 (offset: 1)".
        # Extract the number attached to the likelihood label only; a generic
        # number scan incorrectly selected the final offset value.
        matches = re.findall(
            r"(?:final\s+)?(?:ln\s*likelihood|lnl|likelihood)\s*(?:=|:)?\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)",
            line,
            flags=re.IGNORECASE,
        )
        for token in matches:
            try:
                values.append(float(token))
            except ValueError:
                pass
    if not values:
        return "NA"
    return "{:.12g}".format(values[-1])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--posterior", required=True)
    parser.add_argument("--conditional", default=None)
    parser.add_argument("--stdout", required=True)
    parser.add_argument("--tree-output", default=None)
    parser.add_argument("--model-name", required=True)
    parser.add_argument("--model-parameters", type=int, required=True)
    parser.add_argument("--pillar-map", default=None)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit("ERROR: output directory is not empty: {}".format(outdir))
    outdir.mkdir(parents=True, exist_ok=True)
    n_post = parse_posterior(args.posterior, outdir, args.pillar_map)
    conditional_status = "available"
    if args.conditional and Path(args.conditional).is_file() and Path(args.conditional).stat().st_size > 0:
        n_cond = parse_conditional(args.conditional, outdir)
    else:
        conditional_status = "unavailable"
        n_cond = "NA"
        write_conditional_unavailable(outdir, "POInT conditional branch output was not requested or was unavailable")
    log_likelihood = parse_log_likelihood(args.stdout, args.tree_output)
    write_tsv(
        outdir / "point_model_fit.tsv",
        ["model", "model_parameters", "pillars", "conditional_pillars", "log_likelihood", "conditional_status", "analysis_status"],
        [{"model": args.model_name, "model_parameters": str(args.model_parameters), "pillars": str(n_post), "conditional_pillars": str(n_cond), "log_likelihood": log_likelihood, "conditional_status": conditional_status, "analysis_status": "converged"}],
    )
    print("[OK] Normalized official POInT {}: {} posterior pillars; conditional={} ".format(args.model_name, n_post, conditional_status).rstrip())


if __name__ == "__main__":
    main()
