from pathlib import Path
import csv, json
import numpy as np

import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--synteny-inputs", type=Path, required=True)
parser.add_argument("--results", type=Path, required=True)
args = parser.parse_args()
R = args.results


def rows(p):
    with Path(p).open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def save(p, rs):
    with Path(p).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rs[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rs)


ps = rows(R / "point/pillars/gtr_centered_pillars.tsv")
native = list(
    csv.reader(open(R / "point/order/point_pillars_optimized.txt"), delimiter="\t")
)
assert len(ps) == len(native) - 1 == 12739
out = {
    "pillars": len(ps),
    "unit": "one biological locus per subgenome; conditional on at least one anchored copy in every focal species",
    "tracks": {},
    "block_bootstrap": {},
}
for tr in ["PZI_A", "PZI_B", "PSP_A", "PSP_B", "PPA_A", "PPA_B"]:
    mp = {
        r["representative"]: r["gene_locus"]
        for r in rows(args.synteny_inputs / f"{tr}.locus_map.tsv")
    }
    gs = [r[tr + "_gene"] for r in ps if r[tr + "_gene"]]
    ls = [mp[g] for g in gs]
    assert len(ls) == len(set(ls)), tr
    out["tracks"][tr] = {
        "anchored_loci": len(ls),
        "duplicate_loci": 0,
        "retained_fraction": len(ls) / len(ps),
    }
rng = np.random.default_rng(20261113)
chrlist = sorted({r["gtr_seqid"] for r in ps})
draw = rng.integers(0, len(chrlist), size=(20000, len(chrlist)))
bsrows = []
for sp in ["PZI", "PSP", "PPA"]:
    vals = []
    for ch in chrlist:
        sub = [r for r in ps if r["gtr_seqid"] == ch]
        a = sum(bool(r[sp + "_A_gene"]) for r in sub)
        b = sum(bool(r[sp + "_B_gene"]) for r in sub)
        vals.append([len(sub), a, b])
        bsrows.append(
            dict(
                species=sp,
                chromosome=ch,
                pillars=len(sub),
                A_retained=a,
                B_retained=b,
                B_minus_A=b - a,
            )
        )
    arr = np.array(vals)
    sums = arr[draw].sum(axis=1)
    diff = (sums[:, 2] - sums[:, 1]) / sums[:, 0]
    out["block_bootstrap"][sp] = {
        "blocks": len(chrlist),
        "draws": len(draw),
        "B_minus_A_fraction": (arr[:, 2].sum() - arr[:, 1].sum()) / arr[:, 0].sum(),
        "percentile95": np.quantile(diff, [0.025, 0.975]).tolist(),
        "fraction_draws_positive": float(np.mean(diff > 0)),
    }
node = rows(R / "point/shared/02_tree_mapped/panel_d_tree_mapped_source_data.tsv")
v = {(r["state_retained"], r["branch"]): r for r in node if r["branch"]}
for c in ["A", "B"]:
    for b, parent in [
        ("Root", None),
        ("PPA", "Root"),
        ("NONE", "Root"),
        ("PZI", "NONE"),
        ("PSP", "NONE"),
    ]:
        r = v[(c, b)]
        initial = (
            12739 if parent is None else float(v[(c, parent)]["expected_retained"])
        )
        assert (
            abs(initial - float(r["expected_loss"]) - float(r["expected_retained"]))
            < 3e-6
        )
    out[c + "_mean_tip_retention"] = np.mean(
        [float(v[(c, s)]["expected_retained"]) / 12739 for s in ["PZI", "PSP", "PPA"]]
    )
# Kept separate from the original fitted-model outputs.
(R / "validation").mkdir(exist_ok=True)
save(R / "validation/chromosome_block_retention.tsv", bsrows)
(R / "validation/fractionation_validation.json").write_text(json.dumps(out, indent=2))
print(json.dumps(out, indent=2))
