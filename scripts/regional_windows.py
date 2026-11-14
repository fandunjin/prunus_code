"""Window retention profiles for concordant-region enrichment."""

from pathlib import Path
import argparse
import numpy as np
import pandas as pd

SMOOTH_WINDOW = 11
ROLL_MIN_PERIODS = 6
LOW_THRESHOLD = 0.05
BLOCK_RADIUS = 5
PRIMARY = "decile"
CLASS_DEFS = {
    "decile": ("decile_class", "A_peak_B_valley_decile", "B_peak_A_valley_decile")
}
CONSENSUS_A = "A_peak_B_valley_consensus"
CONSENSUS_B = "B_peak_A_valley_consensus"


def build_windows(
    source: pd.DataFrame, density: pd.DataFrame, species: str
) -> pd.DataFrame:
    d = density[density["species"].eq(species)].copy()
    a = d[d["track"].eq("A")].sort_values("bin").reset_index(drop=True)
    b = d[d["track"].eq("B")].sort_values("bin").reset_index(drop=True)
    if not a["bin"].equals(b["bin"]):
        raise ValueError(f"A/B bins differ for {species}")

    state = source[[f"{species}_A", f"{species}_B"]].fillna(0).astype(int).to_numpy()
    starts = a["pillar_start"].to_numpy(dtype=int)
    ends = a["pillar_end"].to_numpy(dtype=int)
    n_pillars = a["n_pillars"].to_numpy(dtype=int)

    a_count = np.zeros(len(starts), dtype=int)
    b_count = np.zeros(len(starts), dtype=int)
    shared = np.zeros(len(starts), dtype=int)
    one_sided = np.zeros(len(starts), dtype=int)
    both_lost = np.zeros(len(starts), dtype=int)
    chromosome = []
    for i, (start, end) in enumerate(zip(starts, ends)):
        block = state[start:end]
        av = block[:, 0] == 1
        bv = block[:, 1] == 1
        a_count[i] = int(av.sum())
        b_count[i] = int(bv.sum())
        shared[i] = int((av & bv).sum())
        one_sided[i] = int((av ^ bv).sum())
        both_lost[i] = int((~av & ~bv).sum())
        chromosome.append(str(source.iloc[start]["gtr_seqid"]))

    windows = pd.DataFrame(
        {
            "species": species,
            "chromosome": chromosome,
            "bin": a["bin"].to_numpy(dtype=int),
            "pillar_start": starts,
            "pillar_end": ends,
            "n_pillars": n_pillars,
            "A_retained_count": a_count,
            "B_retained_count": b_count,
            "A_retention_fraction": a_count / n_pillars,
            "B_retention_fraction": b_count / n_pillars,
            "shared_AB_count": shared,
            "one_sided_count": one_sided,
            "both_lost_count": both_lost,
        }
    )
    windows["A_smoothed"] = (
        windows["A_retention_fraction"]
        .rolling(SMOOTH_WINDOW, center=True, min_periods=ROLL_MIN_PERIODS)
        .mean()
    )
    windows["B_smoothed"] = (
        windows["B_retention_fraction"]
        .rolling(SMOOTH_WINDOW, center=True, min_periods=ROLL_MIN_PERIODS)
        .mean()
    )
    windows["delta_A_minus_B"] = windows["A_smoothed"] - windows["B_smoothed"]

    eligible = (
        windows["A_retained_count"].gt(0)
        & windows["B_retained_count"].gt(0)
        & windows["A_smoothed"].gt(0)
        & windows["B_smoothed"].gt(0)
        & windows["shared_AB_count"].gt(0)
    )
    windows["eligible"] = eligible

    low = (
        windows["A_retained_count"].eq(0)
        | windows["B_retained_count"].eq(0)
        | windows["A_smoothed"].le(LOW_THRESHOLD)
        | windows["B_smoothed"].le(LOW_THRESHOLD)
    )
    blocked = np.zeros(len(windows), dtype=bool)
    for i in np.flatnonzero(low.to_numpy()):
        blocked[max(0, i - BLOCK_RADIUS) : i + BLOCK_RADIUS + 1] = True
    windows["low_retention_window"] = low
    windows["zero_adjacent_blocked"] = blocked
    clean = eligible & ~blocked
    windows["clean_eligible"] = clean

    delta = windows["delta_A_minus_B"]
    windows["trend_class"] = np.select(
        [~clean, delta.gt(0), delta.lt(0)],
        ["excluded", "A_peak_B_valley", "B_peak_A_valley"],
        default="excluded",
    )

    positive = delta[clean & delta.gt(0)]
    negative = delta[clean & delta.lt(0)]
    pos_cut = float(positive.quantile(0.75)) if len(positive) else np.nan
    neg_cut = float(negative.quantile(0.25)) if len(negative) else np.nan
    windows["extreme_class"] = np.select(
        [clean & delta.ge(pos_cut), clean & delta.le(neg_cut)],
        ["A_peak_B_valley_extreme", "B_peak_A_valley_extreme"],
        default="not_extreme",
    )
    windows["extreme_positive_cutoff"] = pos_cut
    windows["extreme_negative_cutoff"] = neg_cut

    decile_high = float(delta[clean].quantile(0.90))
    decile_low = float(delta[clean].quantile(0.10))
    windows["decile_class"] = np.select(
        [clean & delta.ge(decile_high), clean & delta.le(decile_low)],
        ["A_peak_B_valley_decile", "B_peak_A_valley_decile"],
        default="not_decile",
    )
    windows["decile_high_cutoff"] = decile_high
    windows["decile_low_cutoff"] = decile_low

    a_q75 = float(windows.loc[clean, "A_smoothed"].quantile(0.75))
    a_q25 = float(windows.loc[clean, "A_smoothed"].quantile(0.25))
    b_q75 = float(windows.loc[clean, "B_smoothed"].quantile(0.75))
    b_q25 = float(windows.loc[clean, "B_smoothed"].quantile(0.25))
    windows["strict_class"] = np.select(
        [
            clean & windows["A_smoothed"].ge(a_q75) & windows["B_smoothed"].le(b_q25),
            clean & windows["B_smoothed"].ge(b_q75) & windows["A_smoothed"].le(a_q25),
        ],
        ["A_peak_B_valley_strict", "B_peak_A_valley_strict"],
        default="not_strict",
    )
    return windows


def consensus_regions(windows_all: pd.DataFrame) -> pd.Series:
    pivot = windows_all.pivot_table(
        index="bin", columns="species", values="decile_class", aggfunc="first"
    )
    consensus_a = pivot.eq("A_peak_B_valley_decile").all(axis=1)
    consensus_b = pivot.eq("B_peak_A_valley_decile").all(axis=1)
    return pd.Series(
        np.select(
            [consensus_a, consensus_b],
            [CONSENSUS_A, CONSENSUS_B],
            default="not_consensus",
        ),
        index=pivot.index,
        name="consensus_region",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--density", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    source = pd.read_csv(args.source, sep="\t")
    density = pd.read_csv(args.density, sep="\t")
    windows = pd.concat(
        [build_windows(source, density, s) for s in ["PZI", "PSP", "PPA"]],
        ignore_index=True,
    )
    regions = consensus_regions(windows)
    windows["consensus_region"] = windows["bin"].map(regions)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    windows.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
