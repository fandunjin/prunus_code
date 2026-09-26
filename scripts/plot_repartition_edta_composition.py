#!/usr/bin/env python3
"""Plot EDTA repeat composition for repartitioned Prunus genomes.

The parser accepts both EDTA summary styles present in this project: the
compact ``LTR/Copia`` style and the indented ``LTR``/``Copia`` style.  The
two panels show total genome size on the left and mutually exclusive
EDTA-annotated TE subclasses on the right. Source paths and parser QC are
written alongside the figure so that fallback inputs remain explicit.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

CATEGORIES = [
    "LTR/Gypsy", "LTR/Copia", "Other LTR", "TIR/DNA", "Helitron",
    "LINE", "SINE", "Other TE", "Non-TE",
]

# A dark, restrained palette with distinct hues and grayscale contrast.
COLORS = {
    "LTR/Gypsy": "#214E70",
    "LTR/Copia": "#743C5C",
    "Other LTR": "#9A6680",
    "TIR/DNA": "#176B6E",
    "Helitron": "#6F7F3B",
    "LINE": "#A65E2E",
    "SINE": "#C0834A",
    "Other TE": "#5D536F",
    "Non-TE": "#B8B8B8",
}

plt = None
np = None


def configure_plotting():
    global plt, np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib as mpl
    import matplotlib.pyplot as pyplot
    import numpy as numpy
    plt, np = pyplot, numpy
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7,
        "axes.linewidth": 0.75,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    })


def parse_summary(path: Path) -> dict:
    genome_bp = None
    total_interspersed_bp = None
    classes = defaultdict(int)
    parent = ""
    row_re = re.compile(r"^\s*(\S+)\s+(\d+)\s+(\d+)\s+([0-9.]+)%\s*$")

    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            m = re.search(r"Total Length:\s*([0-9,]+)\s+bp", line)
            if m:
                genome_bp = int(m.group(1).replace(",", ""))
                continue
            m = re.search(r"total interspersed\s+\d+\s+([0-9,]+)(?:\s+bp)?", line, flags=re.IGNORECASE)
            if m:
                total_interspersed_bp = int(m.group(1).replace(",", ""))
                continue

            # Parent rows in the older EDTA format: LTR -- -- --.
            pm = re.match(r"^([A-Za-z][A-Za-z0-9_-]*)\s+--\s+--\s+--", line)
            if pm:
                parent = pm.group(1)
                continue
            rm = row_re.match(line)
            if not rm:
                continue
            token, _count, bp, _pct = rm.groups()
            # A row with no indentation is a new top-level token.
            if line and not line[0].isspace():
                parent = ""
            full = f"{parent}/{token}" if parent and "/" not in token and token not in {"repeat_fragment"} else token
            category = normalize_category(full, parent)
            if category:
                classes[category] += int(bp)

    if genome_bp is None:
        raise ValueError(f"Total Length missing from {path}")
    if total_interspersed_bp is None:
        raise ValueError(f"total interspersed missing from {path}")
    mapped = sum(classes.values())
    return {
        "genome_bp": genome_bp,
        "te_bp_summary": total_interspersed_bp,
        "classes": dict(classes),
        "mapped_class_bp": mapped,
        "unmapped_te_bp": total_interspersed_bp - mapped,
    }


def parse_gff_union(path: Path) -> dict:
    """Assign overlapping EDTA features once and return class-wise unions.

    EDTA's ``bpMasked`` class totals count nested hits more than once.  For a
    stacked genome plot, each atomic interval is assigned to the active hit
    with the highest RepeatMasker score, making the class totals mutually
    exclusive and comparable with the genome length.
    """
    by_seq = defaultdict(list)
    attr_re = re.compile(r"(?:^|;)classification=([^;]+)")
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            try:
                start, end = int(fields[3]) - 1, int(fields[4])
            except ValueError:
                continue
            if end <= start:
                continue
            match = attr_re.search(fields[8])
            classification = match.group(1) if match else fields[2]
            category = normalize_category(classification, "")
            if category is None:
                continue
            try:
                score = float(fields[5]) if fields[5] not in {".", ""} else 0.0
            except ValueError:
                score = 0.0
            by_seq[fields[0]].append((start, end, score, end - start, category))

    classes = defaultdict(int)
    for intervals in by_seq.values():
        points = sorted({point for start, end, *_ in intervals for point in (start, end)})
        starts, ends = defaultdict(list), defaultdict(list)
        for idx, (start, end, score, span, category) in enumerate(intervals):
            starts[start].append(idx)
            ends[end].append(idx)
        active = set()
        for left, right in zip(points, points[1:]):
            active.difference_update(ends[left])
            active.update(starts[left])
            if left == right or not active:
                continue
            winner = max(active, key=lambda idx: (intervals[idx][2], intervals[idx][3]))
            classes[intervals[winner][4]] += right - left
    return dict(classes)


def normalize_category(token: str, parent: str) -> str | None:
    t = token.lower().replace("_", "-")
    p = parent.lower().replace("_", "-")
    if "gypsy" in t:
        return "LTR/Gypsy"
    if "copia" in t:
        return "LTR/Copia"
    if t.startswith("ltr/") or p == "ltr":
        return "Other LTR"
    if "helitron" in t:
        return "Helitron"
    if t.startswith("tir/") or t.startswith("dna/") or p in {"tir", "dna"}:
        return "TIR/DNA"
    if t.startswith("line") or p == "line":
        return "LINE"
    if t.startswith("sine") or p == "sine" or "trna" in t:
        return "SINE"
    if t in {"snrna/unknown", "rrna/unknown", "repeat-fragment", "repeat_fragment"}:
        return "Other TE"
    if t.startswith("nonltr/") or t.startswith("non-tir/") or t.startswith("nonltr"):
        return "Other TE"
    if t.startswith("unknown/") or t == "unknown":
        return "Other TE"
    # Explicitly exclude rDNA, 45S and generic RNA from interspersed TE bars.
    if "rdna" in t or t.startswith("45s") or t.startswith("rrna"):
        return None
    return "Other TE" if "/" in t else None


def choose_sources(input_root: Path, discovery_root: Path) -> list[tuple[str, str, Path, str]]:
    rows: list[tuple[str, str, Path, str]] = []
    # Repartitioned haplotypes: use their own EDTA summaries when available.
    for species, pattern, group in (
        ("pzi", "pzi_hapA1.fa.mod.EDTA.TEanno.sum", "PZI"),
        ("pzi", "pzi_hapB1.fa.mod.EDTA.TEanno.sum", "PZI"),
    ):
        found = sorted(input_root.glob(f"pzi/edta/*/{pattern}"))
        if not found:
            raise FileNotFoundError(f"missing PZI EDTA summary for {pattern}")
        path = found[0]
        label = re.search(r"hap([AB]\d)", path.name).group(1)
        rows.append((f"PZI-{label}", group, path, "repartition_input"))

    # PSP/PPA summaries may be stored as oldrepeat files in this local bundle.
    # Historical EDTA results are stored either beside the repartition inputs
    # or at the project root; accept both layouts for local and HPC reruns.
    old_dirs = [
        input_root.parent / "oldrepeat",
        input_root.parent.parent / "oldrepeat",
        Path("oldrepeat"),
    ]
    old_dirs = list(dict.fromkeys(p.resolve() for p in old_dirs if p.is_dir()))
    for species, labels, group in (
        ("psp", ["A1", "B1"], "PSP"),
        ("ppa", ["A", "B"], "PPA"),
    ):
        for label in labels:
            if species == "psp":
                exact = f"{species}_hap{label}.fa.mod.EDTA.TEanno.sum"
            else:
                exact = f"{species}{label}.fa.mod.EDTA.TEanno.sum"
            candidates = list(input_root.glob(f"{species}/**/{exact}"))
            for old in old_dirs:
                candidates.extend(old.glob(exact))
            if not candidates:
                raise FileNotFoundError(f"missing {species}-{label} EDTA summary")
            rows.append((f"{group}-{label}", group, sorted(candidates)[0], "repartition_input_or_oldrepeat"))

    for species, group in (("pcam", "PCAM"), ("ppe", "PPE")):
        candidates = sorted(discovery_root.glob(f"{species}/**/*TEanno.sum"))
        source_kind = "discovery"
        if not candidates:
            candidates = sorted((input_root.parent / "TE").glob(f"{species}*TEanno.sum"))
            for old in old_dirs:
                candidates += sorted(old.glob(f"{species}*TEanno.sum"))
            source_kind = "fallback_0903_TE_or_oldrepeat"
        if not candidates:
            raise FileNotFoundError(f"missing {species} EDTA summary in discovery and fallback locations")
        rows.append((species.upper(), group, candidates[0], source_kind))
    return rows


def save_figure(fig, outdir: Path, stem: str):
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.tiff", dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(outdir / f"{stem}.png", dpi=300, bbox_inches="tight")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-root", type=Path, default=Path("0903/haploid_repartition/inputs"))
    ap.add_argument("--discovery-root", type=Path, default=Path("analysis/prunus_pan_te/discovery"))
    ap.add_argument("--output-dir", type=Path, default=Path("0903/figures/figure1c_te_composition"))
    ap.add_argument("--stats-only", action="store_true", help="write statistics and description without importing plotting packages")
    args = ap.parse_args()

    sources = choose_sources(args.input_root, args.discovery_root)
    records = []
    for label, group, path, source_kind in sources:
        parsed = parse_summary(path)
        genome = parsed["genome_bp"]
        gff_path = path.with_name(path.name.replace(".sum", ".gff3"))
        if not gff_path.is_file():
            raise FileNotFoundError(f"missing EDTA GFF3 paired with {path}: {gff_path}")
        classes = parse_gff_union(gff_path)
        te = sum(classes.values())
        if te > genome:
            raise ValueError(f"TE span exceeds genome for {label}: {te}>{genome}")
        records.append({
            "label": label, "group": group, "source_kind": source_kind,
            "source": str(path), "gff_source": str(gff_path), "genome_bp": genome, "te_bp": te,
            "te_bp_summary": parsed["te_bp_summary"], "non_te_bp": genome - te,
            "mapped_class_bp": te, "gff_minus_summary_bp": te - parsed["te_bp_summary"],
            "classes": classes,
        })

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fields = ["label", "group", "source_kind", "source", "gff_source", "genome_bp", "te_bp", "te_bp_summary", "non_te_bp", "mapped_class_bp", "gff_minus_summary_bp"]
    with (args.output_dir / "edta_composition.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields + CATEGORIES[:-1], delimiter="\t")
        writer.writeheader()
        for rec in records:
            row = {key: rec[key] for key in fields}
            row.update({cat: rec["classes"].get(cat, 0) for cat in CATEGORIES[:-1]})
            writer.writerow(row)

    manifest = {
        "definition": "EDTA GFF3 intervals are de-overlapped by highest RepeatMasker score and split by normalized classes; Non-TE = Total Length - de-overlapped TE union",
        "warning": "This is an independent-EDTA composition comparison, not the standardized pan-TE library result used for the final cross-species main analysis. PZI-A1/B1 and PSP-A1/B1 are representative haplotypes; PPA-A/B are the two repartitioned packages.",
        "targets": [{"label": r["label"], "group": r["group"], "source_kind": r["source_kind"], "source": r["source"], "gff_source": r["gff_source"]} for r in records],
    }
    with (args.output_dir / "edta_composition_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)

    # Plain-text result summary for the manuscript drafting stage.
    pzi = [r for r in records if r["group"] == "PZI"]
    psp = [r for r in records if r["group"] == "PSP"]
    mean = lambda rows, key: sum(r[key] for r in rows) / len(rows) / 1e6
    pct = lambda rows: 100 * sum(r["te_bp"] for r in rows) / sum(r["genome_bp"] for r in rows)
    desc = (
        f"EDTA independently annotated {mean(pzi, 'te_bp') / mean(pzi, 'genome_bp') * 100:.1f}% of the PZI haplotype packages on average, "
        f"compared with {mean(psp, 'te_bp') / mean(psp, 'genome_bp') * 100:.1f}% for PSP. "
        f"The corresponding mean genome sizes were {mean(pzi, 'genome_bp'):.1f} Mb and {mean(psp, 'genome_bp'):.1f} Mb. "
        "The right-hand stacked bars show how EDTA-annotated sequence is distributed among LTR, TIR/DNA, Helitron, LINE, SINE and other repeat classes; the left-hand bars show total genome size. "
        "PCAM and PPE are included as diploid Prunus references. These values are independent-EDTA estimates and should be interpreted separately from the standardized pan-TE library comparison.\n"
    )
    desc += (
        f"\n中文结果说明：以 PZI-A1、PZI-B1 和 PSP-A1、PSP-B1 作为 A/B 亚基因组代表，独立 EDTA 注释显示，PZI 代表包平均约 {pct(pzi):.1f}% 的基因组序列被识别为转座相关重复，"
        f"高于 PSP（{pct(psp):.1f}%）、PPA（{pct([r for r in records if r['group'] == 'PPA']):.1f}%）、"
        f"PCAM（{pct([r for r in records if r['group'] == 'PCAM']):.1f}%）和 PPE（{pct([r for r in records if r['group'] == 'PPE']):.1f}%）。PZI 的基因组长度也明显更大，"
        f"PZI 代表包平均约 {mean(pzi, 'genome_bp'):.1f} Mb，而 PSP 平均约 {mean(psp, 'genome_bp'):.1f} Mb。PZI 的扩张主要表现为 LTR、TIR/DNA、LINE 及其他重复类别的增加。"
        "由于 EDTA GFF3 中存在嵌套重复，本图对重叠区按最高 RepeatMasker score 进行一次性归属，故右侧 TE 堆叠部分严格加和为 TE 注释总长度。"
        "PCAM 和 PPE 在当前本地目录中使用独立 EDTA 结果作为 discovery 结果缺失时的回退输入；该图适合作为独立 EDTA 的探索性比较，"
        "最终跨物种主结论仍应以统一 pan-TE library 的 RepeatMasker 结果为准。\n"
    )
    (args.output_dir / "figure1c_edta_te_description.txt").write_text(desc, encoding="utf-8")
    if args.stats_only:
        print(json.dumps({"targets": len(records), "output_dir": str(args.output_dir), "mode": "stats-only"}, ensure_ascii=False))
        return

    configure_plotting()
    target_order = [rec["label"] for rec in records]
    fig, axes = plt.subplots(
        1, 2, figsize=(7.2, 4.5),
        gridspec_kw={"width_ratios": (1.0, 1.35), "wspace": 0.34},
        sharex=True,
    )
    ax_genome, ax_te = axes
    x = np.arange(len(records))
    genome_mb = np.array([rec["genome_bp"] / 1e6 for rec in records])
    te_mb = np.array([rec["te_bp"] / 1e6 for rec in records])

    # Left: total genome size for each representative package.
    ax_genome.bar(x, genome_mb, width=0.72, color="#4A5568", linewidth=0)
    for i, value in enumerate(genome_mb):
        ax_genome.text(i, value + max(genome_mb) * 0.018, f"{value:.0f}",
                       ha="center", va="bottom", fontsize=5.8)
    ax_genome.set_ylabel("Genome size (Mb)")
    ax_genome.set_ylim(0, max(genome_mb) * 1.14)
    ax_genome.grid(axis="y", color="#E4E4E4", linewidth=0.5)
    ax_genome.set_axisbelow(True)
    ax_genome.text(-0.14, 1.04, "c", transform=ax_genome.transAxes,
                   fontsize=10, fontweight="bold", va="bottom")

    # Right: TE sequence only; Non-TE is intentionally omitted.
    bottom = np.zeros(len(records))
    for category in CATEGORIES[:-1]:
        vals = np.array([rec["classes"].get(category, 0) / 1e6 for rec in records])
        ax_te.bar(x, vals, bottom=bottom, width=0.72, color=COLORS[category],
                  label=category, linewidth=0)
        bottom += vals
    for i, rec in enumerate(records):
        ax_te.text(i, bottom[i] + max(te_mb) * 0.018,
                   f"{100 * rec['te_bp'] / rec['genome_bp']:.1f}%",
                   ha="center", va="bottom", fontsize=5.8)
    ax_te.set_ylabel("TE-annotated sequence (Mb)")
    ax_te.set_ylim(0, max(te_mb) * 1.18)
    ax_te.grid(axis="y", color="#E4E4E4", linewidth=0.5)
    ax_te.set_axisbelow(True)
    ax_te.legend(ncol=2, loc="lower center", bbox_to_anchor=(0.5, 1.01),
                 fontsize=5.2, columnspacing=0.9, handlelength=1.2)
    ax_te.text(-0.10, 1.04, "d", transform=ax_te.transAxes,
               fontsize=10, fontweight="bold", va="bottom")

    for ax in axes:
        ax.set_xticks(x, target_order, rotation=45, ha="right")
        ax.tick_params(axis="both", length=2.5, width=0.6, pad=2)
        ax.spines["left"].set_linewidth(0.75)
        ax.spines["bottom"].set_linewidth(0.75)
    fig.subplots_adjust(left=0.09, right=0.99, bottom=0.25, top=0.80)
    save_figure(fig, args.output_dir, "figure1c_edta_te_composition")

    print(json.dumps({"targets": len(records), "output_dir": str(args.output_dir)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
