"""Floral-family enrichment in the CAFE5 expansion/contraction results.

Inputs:
  analysis/cafe5/main/k3/{Gamma_change.tab, Gamma_family_results.txt, Gamma_clade_results.txt}
  analysis/family_turnover_20260912/og_turnover.tsv  (floral flag per OG)
  analysis/family_turnover_20260912/top_expanded_floral.tsv
Outputs:
  analysis/cafe5/main/floral_cafe5_summary.tsv
  analysis/cafe5/main/cafe5_clade_branch_summary.tsv
"""
from __future__ import annotations

import csv
import io
import sys
from pathlib import Path

import numpy as np
from scipy.stats import fisher_exact, mannwhitneyu

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

ROOT = Path(r"${PROJECT_ROOT}")
if not ROOT.exists():
    ROOT = Path("${PROJECT_ROOT}")
K3 = ROOT / "analysis" / "cafe5" / "main" / "k3"
FT = ROOT / "analysis" / "family_turnover_20260912"

# ---- read tables
change = list(csv.DictReader(open(K3 / "Gamma_change.tab", encoding="utf-8"), delimiter="\t"))
fam_res = {r["#FamilyID"]: r for r in csv.DictReader(open(K3 / "Gamma_family_results.txt", encoding="utf-8"), delimiter="\t")}
turn = list(csv.DictReader(open(FT / "og_turnover.tsv", encoding="utf-8"), delimiter="\t"))
floral_ids = {r["og"] for r in turn if r["floral"].strip()}
print("CAFE5 families:", len(change), "| floral OGs (turnover subset):", len(floral_ids))

branch_cols = [c for c in change[0] if c != "FamilyID"]
changes = np.array([[float(r[c]) for c in branch_cols] for r in change])
fam_ids = [r["FamilyID"] for r in change]
floral_mask = np.array([f in floral_ids for f in fam_ids])
print("floral families present in CAFE5 table:", int(floral_mask.sum()))

# ---- per-branch summary
rows = []
print(f"\n{'branch':8s} {'n_up':>6s} {'n_dn':>6s} | floral up/dn | background up/dn | fisher p (frac any change)")
for j, c in enumerate(branch_cols):
    up = changes[:, j] > 0
    dn = changes[:, j] < 0
    fu, fd = int((up & floral_mask).sum()), int((dn & floral_mask).sum())
    bu, bd = int((up & ~floral_mask).sum()), int((dn & ~floral_mask).sum())
    # enrichment of any change (up or down) among floral vs background
    tab = [[fu + fd, int(floral_mask.sum()) - (fu + fd)],
           [bu + bd, int((~floral_mask).sum()) - (bu + bd)]]
    try:
        _, p = fisher_exact(tab)
    except Exception:
        p = float("nan")
    rows.append({"branch": c, "families_changed_up": int(up.sum()), "families_changed_down": int(dn.sum()),
                 "floral_up": fu, "floral_down": fd, "background_up": bu, "background_down": bd,
                 "fisher_p_any_change": f"{p:.4g}"})
    print(f"{c:8s} {int(up.sum()):6d} {int(dn.sum()):6d} | {fu:4d}/{fd:<4d}      | {bu:5d}/{bd:<5d}      | p={p:.3g}")

with (ROOT / "analysis" / "cafe5" / "main" / "cafe5_clade_branch_summary.tsv").open("w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)

# ---- racemose stem branch <11>: floral vs background effect sizes
j11 = branch_cols.index("<11>")
fl = changes[floral_mask, j11]
bg = changes[~floral_mask, j11]
print(f"\nracemose stem <11>: floral n={len(fl)} median={np.median(fl):.1f} mean={fl.mean():.2f} | "
      f"background n={len(bg)} median={np.median(bg):.1f} mean={bg.mean():.2f}")
mw = mannwhitneyu(fl, bg, alternative="two-sided")
print("  Mann-Whitney p (floral vs background change on <11>):", f"{mw.pvalue:.4g}")
print("  floral > 0 on <11>:", int((fl > 0).sum()), "/", len(fl), "| background > 0:",
      int((bg > 0).sum()), "/", len(bg))

# ---- CAFE5 family-level significance
sig_set = {f for f, r in fam_res.items() if r["Significant at 0.05"] == "y"}
sig_mask = np.array([f in sig_set for f in fam_ids])
fl_present = int(floral_mask.sum())
bg_present = len(change) - fl_present
fl_sig = int((sig_mask & floral_mask).sum())
bg_sig = int((sig_mask & ~floral_mask).sum())
tab = [[fl_sig, fl_present - fl_sig], [bg_sig, bg_present - bg_sig]]
_, p_sig = fisher_exact(tab)
print(f"\nfamily-level significance (p<0.05): floral {fl_sig}/{fl_present} vs background {bg_sig}/{bg_present}; Fisher p={p_sig:.4g}")

# ---- top expanded floral families: their CAFE5 changes
top = list(csv.DictReader(open(FT / "top_expanded_floral.tsv", encoding="utf-8"), delimiter="\t"))
print("\ntop expanded floral OGs (turnover) — CAFE5 change on key branches:")
for r in top[:8]:
    og_id = r.get("og")
    if og_id in fam_ids:
        j = fam_ids.index(og_id)
        vals = {c: int(changes[j, k]) for k, c in enumerate(branch_cols)}
        print(f"  {og_id} ({r.get('floral')}, ratio={r.get('retention_ratio')}): {vals}")
    else:
        print(f"  {og_id} ({r.get('floral')}): not in CAFE5 table (filtered)")

# ---- save floral detail table
detail = []
for i, f in enumerate(fam_ids):
    if not floral_mask[i]:
        continue
    res = fam_res.get(f, {})
    detail.append({"FamilyID": f, **{c: int(changes[i, k]) for k, c in enumerate(branch_cols)},
                   "cafe5_p": res.get("pvalue", ""), "cafe5_significant": res.get("Significant at 0.05", "")})
with (ROOT / "analysis" / "cafe5" / "main" / "floral_cafe5_summary.tsv").open("w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(detail[0]), delimiter="\t")
    w.writeheader()
    w.writerows(detail)
print("\nwrote floral_cafe5_summary.tsv with", len(detail), "floral families")
