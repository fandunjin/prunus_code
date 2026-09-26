#!/usr/bin/env python3
"""Summarize two MCMCtree chains for the fixed Prunus species tree.

The MCMCtree time unit is 100 Myr.  This script converts ages to Ma, combines
two independent posterior chains, and reports simple split-Rhat plus a
batch-means effective sample-size approximation for each dated node.
"""

import argparse
import csv
import math
from pathlib import Path
from statistics import mean, variance


NODE_CLADE = {
    "t_n13": "Prunus crown",
    "t_n14": "racemose crown (ppa,psp,pzi)",
    "t_n15": "psp+pzi",
    "t_n16": "Prunus excluding racemose lineage",
    "t_n17": "pmal+pav+pspe+pcam",
    "t_n18": "pav+pspe+pcam",
    "t_n19": "pspe+pcam",
    "t_n20": "psal+par+pmu+pdu+ppe",
    "t_n21": "psal+par+pmu",
    "t_n22": "par+pmu",
    "t_n23": "pdu+ppe",
}

CALIBRATION = {
    "t_n13": "P. wutuensis: >55 Ma, <150 Ma",
    "t_n14": "P. hirsutipetala: >35 Ma, <55 Ma",
}


def read_log(path):
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise ValueError("No posterior records found: %s" % path)
    cols = [key for key in rows[0] if key.startswith("t_n")]
    values = {key: [float(row[key]) for row in rows] for key in cols}
    return values


def quantile(values, probability):
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower)


def hpd_interval(values, credibility=0.95):
    """Return the shortest empirical posterior interval at ``credibility``."""
    ordered = sorted(values)
    window = int(math.ceil(credibility * len(ordered)))
    if window < 2 or window > len(ordered):
        raise ValueError("Invalid HPD interval window")
    start = min(
        range(len(ordered) - window + 1),
        key=lambda index: ordered[index + window - 1] - ordered[index],
    )
    return ordered[start], ordered[start + window - 1]


def split_rhat(chains):
    halves = []
    for chain in chains:
        half = len(chain) // 2
        halves.extend([chain[:half], chain[-half:]])
    m = len(halves)
    n = len(halves[0])
    within = mean(variance(part) for part in halves)
    between = n * variance([mean(part) for part in halves])
    variance_plus = ((n - 1) / n) * within + (between / n)
    return math.sqrt(variance_plus / within) if within > 0 else float("nan")


def batch_ess(chains, batch_size=100):
    """A conservative ESS approximation based on variance of batch means."""
    all_values = [value for chain in chains for value in chain]
    total = len(all_values)
    if total < batch_size * 4:
        return float("nan")
    batches = []
    for chain in chains:
        count = len(chain) // batch_size
        batches.extend(mean(chain[i * batch_size:(i + 1) * batch_size]) for i in range(count))
    raw_variance = variance(all_values)
    batch_variance = variance(batches)
    if raw_variance <= 0 or batch_variance <= 0:
        return float("nan")
    ess = total * raw_variance / (batch_size * batch_variance)
    return min(float(total), ess)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--chain1", required=True)
    parser.add_argument("--chain2", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--unit-ma", type=float, default=100.0)
    args = parser.parse_args()

    chain1 = read_log(args.chain1)
    chain2 = read_log(args.chain2)
    if set(chain1) != set(chain2):
        raise SystemExit("Chains have different dated-node columns")

    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "node", "clade", "calibration", "posterior_samples_per_chain",
            "mean_ma", "median_ma", "hpd2.5_ma", "hpd97.5_ma",
            "chain1_mean_ma", "chain2_mean_ma", "split_rhat",
            "batch_ess_approx",
        ])
        for node in sorted(chain1, key=lambda value: int(value.split("n")[-1])):
            values = chain1[node] + chain2[node]
            hpd_low, hpd_high = hpd_interval(values)
            writer.writerow([
                node,
                NODE_CLADE.get(node, "unmapped"),
                CALIBRATION.get(node, "none"),
                len(chain1[node]),
                "%.3f" % (mean(values) * args.unit_ma),
                "%.3f" % (quantile(values, 0.5) * args.unit_ma),
                "%.3f" % (hpd_low * args.unit_ma),
                "%.3f" % (hpd_high * args.unit_ma),
                "%.3f" % (mean(chain1[node]) * args.unit_ma),
                "%.3f" % (mean(chain2[node]) * args.unit_ma),
                "%.5f" % split_rhat([chain1[node], chain2[node]]),
                "%.1f" % batch_ess([chain1[node], chain2[node]]),
            ])
    print("nodes=%d samples_per_chain=%d output=%s" % (len(chain1), len(next(iter(chain1.values()))), output))


if __name__ == "__main__":
    main()
