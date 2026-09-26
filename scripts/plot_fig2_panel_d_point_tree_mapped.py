#!/usr/bin/env python3
"""Map native POInT conditional-loss expectations onto two copy trees.

This renderer follows the visual logic of the reference Fig. 3d: a shared
duplication is shown once, followed by one topology for each POInT singleton
state.  It deliberately uses Copy1/Copy2, rather than modern A/B labels,
because POInT tracking states are local syntenic-copy states.  With
``--fixed-ab`` the same renderer labels the two tracks as physical A and B
and requires a POInT run made with the optional identity-permutation model.
"""

from __future__ import print_function

import argparse
import csv
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.linewidth": 0.75,
})


TOPOLOGY_BRANCHES = ("Root", "PPA", "NONE", "PZI", "PSP")
COPY_STATES = {
    "C1": {
        "loss_transition": "Dupl-->Copy2",
        "color": "#5B7FCA",
        "fill": "#E9ECF7",
        "state": "Copy1",
    },
    "C2": {
        "loss_transition": "Dupl-->Copy1",
        "color": "#E28E79",
        "fill": "#FBE9D8",
        "state": "Copy2",
    },
}
AB_STATES = {
    "C1": dict(COPY_STATES["C1"], state="A"),
    "C2": dict(COPY_STATES["C2"], state="B"),
}


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def require_file(path, label):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit("ERROR: missing {}: {}".format(label, path))


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="ascii", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def value_label(retained, loss):
    return "{:,.1f} (-{:,.1f})".format(retained, loss)


def load_transition_values(rows, n_pillars, state_defs):
    values = {}
    for row in rows:
        branch = row.get("branch", "")
        transition = row.get("transition", "")
        if branch in TOPOLOGY_BRANCHES and transition in ("Dupl-->Copy1", "Dupl-->Copy2"):
            values[(branch, transition)] = float(row["mean_probability"]) * n_pillars
    required = {(branch, state["loss_transition"])
                for branch in TOPOLOGY_BRANCHES for state in state_defs.values()}
    missing = sorted(required - set(values))
    if missing:
        raise SystemExit("ERROR: conditional POInT summary lacks transitions: {}".format(missing))
    return values


def validate_tree(path):
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    expected = {"PZI", "PSP", "PPA"}
    missing = sorted(taxon for taxon in expected if taxon not in text)
    if missing:
        raise SystemExit("ERROR: POInT tree is missing expected taxa: {}".format(", ".join(missing)))


def calculate_copy_paths(n_pillars, transition_values, state_defs):
    """Return retained-copy expectation at each descendant node.

    A Dupl-->Copy1 state retains Copy1 and therefore represents loss of Copy2;
    conversely, Dupl-->Copy2 represents loss of Copy1.  Starting from one copy
    per pillar at the shared WGD, loss expectations are propagated down the
    rooted PPA,(PSP,PZI) topology separately for C1 and C2.
    """
    results = {}
    for copy, info in state_defs.items():
        loss_transition = info["loss_transition"]
        losses = {branch: transition_values[(branch, loss_transition)] for branch in TOPOLOGY_BRANCHES}
        retained = {}
        retained["Root"] = n_pillars - losses["Root"]
        retained["PPA"] = retained["Root"] - losses["PPA"]
        retained["NONE"] = retained["Root"] - losses["NONE"]
        retained["PZI"] = retained["NONE"] - losses["PZI"]
        retained["PSP"] = retained["NONE"] - losses["PSP"]
        results[copy] = {"losses": losses, "retained": retained}
    return results


def draw_copy_tree(ax, y_root, copy, path_values, n_pillars, branch_lengths, state_defs):
    """Draw a PPA,(PSP,PZI) copy tree and branch-level POInT expectations."""
    info = state_defs[copy]
    color = info["color"]
    fill = info["fill"]
    # Fixed coordinates keep upper and lower trees structurally identical.
    x_star, x_root, x_split, x_tip = 0.13, 0.33, 0.55, 0.82
    ppa_y = y_root + 0.125
    psp_y = y_root
    pzi_y = y_root - 0.125
    pair_y = (psp_y + pzi_y) / 2.0

    ax.add_patch(plt.Rectangle((0.25, y_root - 0.18), 0.70, 0.36, color=fill, ec="none", zorder=0))
    line = {"color": "#242424", "lw": 1.05, "solid_capstyle": "butt", "zorder": 2}
    ax.plot([x_star, x_root], [y_root, y_root], **line)
    ax.plot([x_root, x_root], [pair_y, ppa_y], **line)
    ax.plot([x_root, x_tip], [ppa_y, ppa_y], **line)
    ax.plot([x_root, x_split], [pair_y, pair_y], **line)
    ax.plot([x_split, x_split], [pzi_y, psp_y], **line)
    ax.plot([x_split, x_tip], [psp_y, psp_y], **line)
    ax.plot([x_split, x_tip], [pzi_y, pzi_y], **line)

    annotations = {
        "Root": ((x_star + x_root) / 2.0, y_root + 0.027),
        "PPA": ((x_root + x_tip) / 2.0, ppa_y + 0.023),
        "NONE": ((x_root + x_split) / 2.0, pair_y + 0.023),
        "PSP": ((x_split + x_tip) / 2.0, psp_y + 0.023),
        "PZI": ((x_split + x_tip) / 2.0, pzi_y + 0.023),
    }
    for branch in TOPOLOGY_BRANCHES:
        x, y = annotations[branch]
        ax.text(x, y, value_label(path_values["retained"][branch], path_values["losses"][branch]),
                ha="center", va="bottom", fontsize=5.7, color="#242424")
    length_annotations = {
        "Root": ((x_star + x_root) / 2.0, y_root - 0.035),
        "PPA": ((x_root + x_tip) / 2.0, ppa_y - 0.035),
        "NONE": ((x_root + x_split) / 2.0, pair_y - 0.035),
        "PSP": ((x_split + x_tip) / 2.0, psp_y - 0.035),
        "PZI": ((x_split + x_tip) / 2.0, pzi_y - 0.035),
    }
    for branch in TOPOLOGY_BRANCHES:
        x, y = length_annotations[branch]
        ax.text(x, y, "{:.4f}".format(branch_lengths[branch]),
                ha="center", va="top", fontsize=4.8, color="#4D4D4D")
    for name, y in (("PPA", ppa_y), ("PSP", psp_y), ("PZI", pzi_y)):
        ax.text(x_tip + 0.017, y, name, ha="left", va="center", fontsize=6.6, fontstyle="italic")
    ax.text(0.925, y_root, info["state"], ha="center", va="center", fontsize=7.4,
            color=color, fontweight="bold")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fit", required=True, help="Native POInT model-fit summary")
    parser.add_argument("--transition", required=True, help="Native POInT conditional transition summary")
    parser.add_argument("--tree", required=True, help="Optimized native POInT tree")
    parser.add_argument("--model", required=True, help="Optimized native POInT model text")
    parser.add_argument("--partitions", required=False, help="Five shared/independent model comparison")
    parser.add_argument("--fixed-ab", action="store_true",
                        help="Render a run made with the fixed physical A/B identity model")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    state_defs = AB_STATES if args.fixed_ab else COPY_STATES
    required_files = ((args.fit, "POInT fit"), (args.transition, "POInT conditional transitions"),
                      (args.tree, "POInT tree"), (args.model, "POInT model"))
    for path, label in required_files:
        require_file(path, label)
    if not args.fixed_ab:
        if not args.partitions:
            raise SystemExit("ERROR: --partitions is required for the native Copy1/Copy2 renderer")
        require_file(args.partitions, "POInT partition comparison")
    validate_tree(args.tree)

    fit_rows = read_tsv(args.fit)
    transition_rows = read_tsv(args.transition)
    partition_rows = read_tsv(args.partitions) if args.partitions else []
    expected_model = "reference_WGD_bias_fix_fixedAB" if args.fixed_ab else "reference_WGD_bias_fix"
    selected = [row for row in fit_rows if row.get("model") == expected_model]
    if len(selected) != 1:
        raise SystemExit("ERROR: expected one converged {} POInT fit".format(expected_model))
    fit = selected[0]
    if fit.get("analysis_status", "").lower() != "converged":
        raise SystemExit("ERROR: reference_WGD_bias_fix fit is not converged")
    n_pillars = int(float(fit["pillars"]))
    if n_pillars <= 0:
        raise SystemExit("ERROR: invalid POInT pillar count")
    if (not args.fixed_ab) and (len(partition_rows) != 5 or not any(row.get("model") == "all_shared" for row in partition_rows)):
        raise SystemExit("ERROR: expected the five shared/independent POInT models")
    transition_values = load_transition_values(transition_rows, n_pillars, state_defs)
    copy_paths = calculate_copy_paths(n_pillars, transition_values, state_defs)
    model_text = Path(args.model).read_text(encoding="utf-8", errors="replace")
    model_values = {}
    for key in ("FixRate", "SwitchProb", "Copy1Bias", "Copy2Bias"):
        match = re.search(r"^{}\s+\S+\s+([^\s]+)\s*$".format(key), model_text, flags=re.MULTILINE)
        if match:
            model_values[key.lower()] = match.group(1)
    branch_lengths = {}
    for row in transition_rows:
        branch = row.get("branch", "")
        if branch in TOPOLOGY_BRANCHES:
            branch_lengths[branch] = float(row["branch_length"])
    if set(branch_lengths) != set(TOPOLOGY_BRANCHES):
        raise SystemExit("ERROR: conditional summary lacks all POInT branch lengths")

    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit("ERROR: output directory is not empty: {}".format(outdir))
    outdir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(5.2, 4.15))
    ax.set_xlim(0.03, 1.00)
    ax.set_ylim(0.03, 0.97)
    ax.axis("off")
    # One star explicitly represents the shared-WGD model.  The vertical root
    # line connects the two copy trees, as in the reference figure.
    ax.plot([0.13, 0.13], [0.25, 0.75], color="#242424", lw=1.05, zorder=1)
    ax.scatter([0.13], [0.50], s=150, marker="*", color="#C63D3D", edgecolors="none", zorder=4)
    ax.text(0.145, 0.52, "{:,}".format(n_pillars), ha="left", va="bottom", fontsize=5.7)
    draw_copy_tree(ax, 0.75, "C1", copy_paths["C1"], n_pillars, branch_lengths, state_defs)
    draw_copy_tree(ax, 0.25, "C2", copy_paths["C2"], n_pillars, branch_lengths, state_defs)
    ax.text(0.03, 0.97, "d", ha="left", va="top", fontsize=10, fontweight="bold")
    model_label = "WGD_bias_fix_genome (fixed A/B)" if args.fixed_ab else "WGD_bias_fix_genome"
    ax.text(0.03, 0.08, model_label, ha="left", va="bottom", fontsize=5.2, color="#4D4D4D")
    ax.text(0.03, 0.062, "FixRate={}".format(model_values.get("fixrate", "NA")), ha="left", va="bottom", fontsize=5.0, color="#4D4D4D")
    bias_key = "copy1bias" if args.fixed_ab and "copy1bias" in model_values else "copy2bias"
    bias_label = "Copy1Bias" if bias_key == "copy1bias" else "Copy2Bias"
    ax.text(0.03, 0.046, "{}={}".format(bias_label, model_values.get(bias_key, "NA")), ha="left", va="bottom", fontsize=5.0, color="#4D4D4D")
    if not args.fixed_ab:
        ax.text(0.03, 0.030, "SwitchProb={}".format(model_values.get("switchprob", "NA")), ha="left", va="bottom", fontsize=5.0, color="#4D4D4D")
    ax.text(0.36, 0.08, "Final ln likelihood={}".format(fit["log_likelihood"]), ha="left", va="bottom", fontsize=5.2, color="#4D4D4D")
    for copy, x in (("C1", 0.36), ("C2", 0.63)):
        terminal = [copy_paths[copy]["retained"][branch] / float(n_pillars) for branch in ("PPA", "PZI", "PSP")]
        label = state_defs[copy]["state"]
        ax.text(x, 0.046, "prop. retained {}={:.3f}".format(label, sum(terminal) / 3.0),
                ha="left", va="bottom", fontsize=5.0, color="#4D4D4D")
    ax.text(0.36, 0.030, "numbers: expected retained copies (expected losses); branch lengths below", ha="left", va="bottom", fontsize=5.0, color="#4D4D4D")
    fig.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99)
    prefix = outdir / "figure2_panel_d_point_tree_mapped"
    fig.savefig(str(prefix) + ".svg", bbox_inches="tight")
    fig.savefig(str(prefix) + ".pdf", bbox_inches="tight")
    fig.savefig(str(prefix) + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(str(prefix) + ".tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)

    source_rows = []
    for copy, info in state_defs.items():
        for branch in TOPOLOGY_BRANCHES:
            source_rows.append({
                "section": "tree_mapped_conditional_expectation",
                "copy_track": copy,
                "state_retained": info["state"],
                "loss_transition": info["loss_transition"],
                "branch": branch,
                "initial_pillars": str(n_pillars),
                "expected_loss": "{:.6f}".format(copy_paths[copy]["losses"][branch]),
                "expected_retained": "{:.6f}".format(copy_paths[copy]["retained"][branch]),
                "branch_length": "{:.8f}".format(branch_lengths[branch]),
            })
    for row in partition_rows:
        source_rows.append({
            "section": "shared_independent_model_fit",
            "model": row.get("model", ""),
            "label": row.get("label", ""),
            "log_likelihood": row.get("log_likelihood", ""),
            "AIC": row.get("AIC", ""),
            "BIC": row.get("BIC", ""),
            "delta_AIC": row.get("delta_AIC", ""),
            "delta_BIC": row.get("delta_BIC", ""),
        })
    source_rows.append({
        "section": "selected_point_model",
        "model": expected_model,
        "label": "all shared",
        "log_likelihood": fit.get("log_likelihood", ""),
        "initial_pillars": str(n_pillars),
        "fix_rate": model_values.get("fixrate", ""),
        "copy1_bias": model_values.get("copy1bias", ""),
        "copy2_bias": model_values.get("copy2bias", ""),
        "switch_prob": model_values.get("switchprob", ""),
    })
    fields = ["section", "copy_track", "state_retained", "loss_transition", "branch", "initial_pillars",
              "expected_loss", "expected_retained", "branch_length", "model", "label", "log_likelihood", "AIC", "BIC",
              "delta_AIC", "delta_BIC", "fix_rate", "copy1_bias", "copy2_bias", "switch_prob"]
    write_tsv(outdir / "panel_d_tree_mapped_source_data.tsv", fields, source_rows)
    boundary = [
        ("The two trees are fixed physical A and B duplicate tracks from the identity-permutation POInT model."
         if args.fixed_ab else
         "The two trees are POInT Copy1 and Copy2 singleton-state tracks, not direct modern A/B subgenomes."),
        "Labels are model-based expected retained copies and expected losses from native conditional transition probabilities; they are not observed gene-loss counts.",
        ("The red star marks the shared WGD root in the fixed-tree analysis; no shared-versus-independent AIC/BIC comparison is made in this fixed run."
         if args.fixed_ab else
         "The red star represents the selected all-shared event model. Its AIC/BIC comparison is stored in the accompanying source table, not plotted as a second visual claim."),
        "GTR anchors pillar homology but is not a POInT likelihood terminal and is therefore not plotted as a tip.",
    ]
    (outdir / "analysis_boundary.txt").write_text("\n".join(boundary) + "\n", encoding="ascii")
    complete_label = "POInT fixed A/B tree-mapped Figure 2d completed." if args.fixed_ab else "POInT tree-mapped Figure 2d completed."
    (outdir / "COMPLETE").write_text(complete_label + "\n", encoding="ascii")
    (outdir / "software_and_parameters.txt").write_text(
        "backend=Python matplotlib only\n"
        + "mode={}\n".format("fixed physical A/B identity permutation" if args.fixed_ab else "native Copy1/Copy2 permutations")
        + "pillars={}\n".format(n_pillars)
        + "topology=PPA,(PSP,PZI)\n"
        + "transition_source=point_branch_transition_summary.tsv\n"
        + "tree_source=optimized POInT tree output\n",
        encoding="ascii")
    print("[OK] Wrote tree-mapped POInT panel d: {}".format(outdir))


if __name__ == "__main__":
    main()
