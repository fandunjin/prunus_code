"""Pool locus-level A/B counts for the full and diploid-control CAFE runs."""

import argparse
from pathlib import Path
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--counts", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
counts = pd.read_csv(args.counts, sep="\t").set_index("Orthogroup")
tracks = {
    "pzi": ["pzi_A", "pzi_B"],
    "psp": ["psp_A", "psp_B"],
    "ppa": ["ppa_A", "ppa_B"],
    "pcam": ["pcam"],
    "ppe": ["ppe"],
    "gtr": ["gtr"],
    "rch": ["rch"],
}
pooled = pd.DataFrame(
    {sp: counts[columns].sum(axis=1) for sp, columns in tracks.items()}
)
pooled = pooled[(pooled.gt(0).sum(axis=1) >= 2) & (pooled.max(axis=1) <= 100)]
pooled.insert(
    0, "Family ID", ["LGF" + name.removeprefix("OG") for name in pooled.index]
)
pooled.insert(0, "Desc", "(null)")
args.out.mkdir(parents=True, exist_ok=True)
pooled.to_csv(args.out / "full_counts.tsv", sep="\t", index=False)
pooled[["Desc", "Family ID", "pcam", "ppe", "gtr", "rch"]].to_csv(
    args.out / "diploid_counts.tsv", sep="\t", index=False
)
