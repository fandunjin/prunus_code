#!/usr/bin/env python3
"""Prepare and draw Fig. 5c-f for inflorescence-related gene families."""
from __future__ import annotations
import argparse
import csv
import math
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CANON = ROOT / "analysis" / "canonical"
OF = ROOT / "analysisphylogenetics" / "02_orthofinder" / "main" / "Results_Aug25"
IDMAP = ROOT / "analysisphylogenetics" / "01_frozen_inputs" / "main" / "id_maps"
SYNTENY = ROOT / "analysis" / "wgd_kaks" / "block_representative_YN" / "summary" / "syntenic_pair_plan.tsv"
REPART = ROOT / "haploid_repartition" / "maps" / "pzi_repartition_map.tsv"
JCVI_ROOT = ROOT / "downstream_analysis" / "circos" / "work" / "jcvi"
OUT = ROOT / "article" / "figure5_inflorescence"
PHASES = ["A1", "A2", "B1", "B2"]
TRACK = {"A1": "A", "A2": "A", "B1": "B", "B2": "B"}
TREE_PRIORITY = ["WUS", "KNOX/BELL", "AG", "AP2", "PIN", "SPL", "TCP"]
F5_DISPLAY = ["LFY", "AP1", "FUL", "TFL1/CEN", "WUS", "CLV", "LAX", "KNOX/BELL", "PIN", "AG"]

# Explicitly documented, function-driven candidates. These are not GO terms.
# Boundaries deliberately exclude substrings such as "splicing" from SPL.
CANDIDATES = [
    ("LFY", r"\bLEAFY\b|\bLFY\b"), ("AP1", r"\bAP1\b|\bAPETALA[- ]?1\b"),
    ("FUL", r"\bFUL\b|\bFRUITFULL\b"), ("TFL1/CEN", r"\bTERMINAL FLOWER\b|\bTFL1\b|\bCEN\b|CENTRORADIALIS"),
    ("FT", r"\bFLOWERING LOCUS T\b|\bFT\b"), ("SOC1", r"\bSUPPRESSOR OF OVEREXPRESSION OF CONSTANS\b|\bSOC1\b"),
    ("SVP", r"\bSHORT VEGETATIVE PHASE\b|\bSVP\b"), ("AP2", r"\bAPETALA[- ]?2\b|\bAP2(?:/ERF)?\b"),
    ("SEP", r"\bSEPALLATA\b|\bSEP[1-4]?\b"), ("AGL6", r"\bAGAMOUS[- ]?LIKE[- ]?6\b|\bAGL6\b"),
    ("AG", r"\bFLORAL HOMEOTIC\b|\bAGAMOUS\b|\bAG[1-4]\b"),
    ("WUS", r"\bWUSCHEL\b|\bWUS\b"), ("CLV", r"\bCLAVATA\b|\bCLV[1-3]\b"),
    ("BRC1", r"\bBRANCHED[- ]?1\b|\bBRC1\b"), ("LAX", r"\bLAX PANICLE\b|\bLAX\b"),
    ("KNOX/BELL", r"\bKNOX\b|\bKNAT\d*\b|\bKNOTTED\b|\bBELL(?:1)?\b|BELL HOMEODOMAIN"),
    ("PIN", r"\bPIN[- ]FORMED\b|\bPIN[1-9]?\b|\bPINOID\b|\bPIN[- ]LIKE[S]?\b|\bPILS\d*\b"),
    ("TCP", r"\bTEOSINTE BRANCHED\b|\bTCP\d+[A-Z-]*\b"), ("SPL", r"\bSQUAMOSA PROMOTER BINDING\b|\bSPL\d+[A-Z-]*\b"),
]

@dataclass
class Gene:
    phase: str; gene_id: str; chrom: str; start: int; end: int; strand: str; note: str
    original_id: str = ""; protein_seq: str = ""; stable_id: str = ""; orthogroup: str = ""
    transcript_id: str = ""; old_gene_id: str = ""; mapping_score: float = 0.0
    families: list[str] = field(default_factory=list)

def attrs(text: str) -> dict[str, str]:
    out = {}
    for item in text.rstrip(";\n").split(";"):
        if "=" in item:
            k, v = item.split("=", 1); out[k] = v
    return out

def fasta(path: Path) -> dict[str, str]:
    seqs, key, chunks = {}, None, []
    with path.open(encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                if key is not None: seqs[key] = "".join(chunks)
                key, chunks = line[1:].strip().split()[0], []
            else: chunks.append(line.strip())
    if key is not None: seqs[key] = "".join(chunks)
    return seqs

def parse_gff(phase: str) -> list[Gene]:
    genes = []; by_id = {}
    with (CANON / f"pzi_hap{phase}" / "primary.gff3").open(encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("#"): continue
            p = line.rstrip("\n").split("\t")
            if len(p) != 9: continue
            a = attrs(p[8])
            if p[2] == "gene":
                note = a.get("Note", "")
                fam = [name for name, pat in CANDIDATES if re.search(pat, note, re.I)]
                if fam:
                    g = Gene(phase, a.get("ID", ""), p[0], int(p[3]), int(p[4]), p[6], note, families=fam)
                    genes.append(g); by_id[g.gene_id] = g
            elif p[2] == "mRNA":
                parent = a.get("Parent", "")
                if parent in by_id and not by_id[parent].transcript_id:
                    by_id[parent].transcript_id = a.get("ID", "")
    return genes

def load_phase_seqmap() -> dict[str, dict[str, str]]:
    """Map canonical chrXX_phase names to the old CnA/CnB chromosomes."""
    out = {p: {} for p in PHASES}
    with REPART.open(encoding="utf-8", errors="ignore") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            phase = row.get("haplotype", "").replace("hap", "", 1)
            if phase not in out: continue
            # Canonical GFF uses chr01_A1..chr08_A1 (ordinal phase tracks),
            # whereas frozen GFF uses C1A..C8A. Preserve that ordinal mapping.
            new_seqid = row.get("new_seqid", "")
            m = re.match(r"C(\d+)[AB]$", new_seqid)
            if m:
                out[phase][f"chr{int(m.group(1)):02d}_{phase}"] = new_seqid
    return out

def old_models(phase: str) -> dict[str, list[dict[str, object]]]:
    """Read frozen gene intervals and primary transcript IDs used by OrthoFinder."""
    path = ROOT / "haploid_repartition" / "inputs" / "pzi" / "gff" / f"pzi_hap{phase}.gff3"
    models = defaultdict(list); current = None
    with path.open(encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("#"): continue
            p = line.rstrip("\n").split("\t")
            if len(p) != 9: continue
            a = attrs(p[8])
            if p[2] == "gene":
                current = {"id": a.get("ID", ""), "seqid": p[0], "start": int(p[3]), "end": int(p[4]), "transcript": ""}
                models[p[0]].append(current)
            elif p[2] == "mRNA" and current is not None and a.get("Parent") == current["id"] and not current["transcript"]:
                current["transcript"] = a.get("ID", "")
    return models

def load_id_map(track: str) -> dict[str, str]:
    out = {}
    with (IDMAP / f"pzi_{track}.id_map.tsv").open(encoding="utf-8", errors="ignore") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["original_id"]] = row["stable_id"]; out[row["original_id"].split(".")[0]] = row["stable_id"]
    return out

def map_canonical_to_stable(genes_by_phase: dict[str, list[Gene]]) -> None:
    id_maps = {"A": load_id_map("A"), "B": load_id_map("B")}
    seqmap = load_phase_seqmap()
    for phase, genes in genes_by_phase.items():
        models = old_models(phase)
        # Exact protein sequence bridge is robust to copy-2 coordinate systems.
        seq_to_old: dict[str, list[str]] = defaultdict(list)
        source_phase = phase
        for old, seq in fasta(ROOT / "haploid_repartition" / "inputs" / "pzi" / "pep" / f"pzi_hap{source_phase}.pep.fa").items():
            seq_to_old[seq.rstrip("*")].append(old)
        loc_to_seq: dict[str, str] = {}
        for header, seq in fasta(CANON / f"pzi_hap{phase}" / "primary.pep.fa").items():
            m = re.search(r"(LOC_[0-9]+)", header)
            if m and len(seq) > len(loc_to_seq.get(m.group(1), "")):
                loc_to_seq[m.group(1)] = seq.rstrip("*")
        for gene in genes:
            loc = re.search(r"(LOC_[0-9]+)", gene.gene_id)
            seq = loc_to_seq.get(loc.group(1), "") if loc else ""
            exact = seq_to_old.get(seq, []) if seq else []
            if exact:
                gene.original_id = exact[0]; gene.mapping_score = 1.0
                gene.stable_id = id_maps[TRACK[phase]].get(gene.original_id, id_maps[TRACK[phase]].get(gene.original_id.split(".")[0], ""))
                gene.old_gene_id = gene.original_id.rsplit(".", 1)[0]
                continue
            old_chr = seqmap.get(phase, {}).get(gene.chrom, "")
            best, best_score = None, 0.0
            for model in models.get(old_chr, []):
                overlap = max(0, min(gene.end, model["end"]) - max(gene.start, model["start"]) + 1)
                denom = min(gene.end - gene.start + 1, model["end"] - model["start"] + 1)
                score = overlap / denom if denom else 0.0
                if score > best_score: best, best_score = model, score
            if best is not None and best_score >= 0.50:
                gene.old_gene_id = str(best["id"]); gene.original_id = str(best["transcript"])
                gene.mapping_score = best_score
                gene.stable_id = id_maps[TRACK[phase]].get(gene.original_id, id_maps[TRACK[phase]].get(gene.original_id.split(".")[0], ""))

def load_orthogroups() -> tuple[dict[str, str], dict[str, set[str]]]:
    gene_to_og, og_members = {}, {}
    with (OF / "Orthogroups" / "Orthogroups.txt").open(encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if ":" not in line: continue
            og, rest = line.rstrip("\n").split(":", 1); members = set(rest.split()); og_members[og] = members
            for m in members:
                if m.startswith(("pzi_A__", "pzi_B__")): gene_to_og[m] = og
    return gene_to_og, og_members

def candidate_tables(genes_by_phase, gene_to_og, og_members):
    records = []
    for phase, genes in genes_by_phase.items():
        for g in genes:
            g.orthogroup = gene_to_og.get(g.stable_id, "")
            for family in g.families:
                records.append({"family": family, "phase": phase, "track": TRACK[phase], "gene_id": g.gene_id, "transcript_id": g.transcript_id, "old_gene_id": g.old_gene_id, "original_id": g.original_id, "stable_id": g.stable_id, "orthogroup": g.orthogroup, "mapping_score": f"{g.mapping_score:.3f}", "chrom": g.chrom, "start": g.start, "end": g.end, "strand": g.strand, "note": g.note, "mapping_status": "mapped_to_orthofinder" if g.orthogroup else "annotation_only"})
    summaries = []
    for family in sorted({r["family"] for r in records}):
        fr = [r for r in records if r["family"] == family]; pres = {p: int(any(r["phase"] == p for r in fr)) for p in PHASES}; ac = sum(r["track"] == "A" for r in fr); bc = sum(r["track"] == "B" for r in fr)
        state = "expanded" if ac and bc and (ac > 2 or bc > 2) else "both_retained" if ac and bc else "A_only" if ac else "B_only"
        summaries.append({"family": family, **pres, "A_count": ac, "B_count": bc, "state": state, "orthogroups": ";".join(sorted({r["orthogroup"] for r in fr if r["orthogroup"]})), "evidence": "annotation_note;exact_protein_bridge" if any(r["stable_id"] for r in fr) else "annotation_note_only"})
    bg = []
    for og, members in og_members.items():
        a = sum(m.startswith("pzi_A__") for m in members); b = sum(m.startswith("pzi_B__") for m in members)
        if a or b: bg.append({"orthogroup": og, "A_present": int(a > 0), "B_present": int(b > 0), "state": "both_retained" if a and b else "A_only" if a else "B_only", "A_count": a, "B_count": b})
    return records, summaries, bg

def fisher_two_sided(a, b, c, d):
    n = a + b + c + d
    def prob(x): return math.comb(a + b, x) * math.comb(c + d, a + c - x) / math.comb(n, a + c)
    p0 = prob(a); lo, hi = max(0, a - d), min(a + b, a + c)
    return min(1.0, sum(prob(x) for x in range(lo, hi + 1) if prob(x) <= p0 + 1e-15))

def effect_sizes(summaries, bg):
    out, nbg = [], len(bg)
    bg_map = {x["orthogroup"]: x for x in bg}
    all_ogs = set(bg_map)
    for row in summaries:
        candidates = {x for x in row.get("orthogroups", "").split(";") if x}
        if not candidates: continue
        background = all_ogs - candidates
        for comp, key in [("A_vs_background", "A_present"), ("B_vs_background", "B_present")]:
            aa = sum(bg_map[og][key] for og in candidates if og in bg_map)
            bb = len(candidates) - aa
            cc = sum(bg_map[og][key] for og in background)
            dd = len(background) - cc
            OR = (aa + .5) * (dd + .5) / ((bb + .5) * (cc + .5)); se = math.sqrt(sum(1 / (x + .5) for x in (aa, bb, cc, dd)))
            out.append({"family": row["family"], "state": row["state"], "comparison": comp, "odds_ratio": OR, "log2_effect": math.log2(OR), "ci_low": math.log2(OR) - 1.96 * se, "ci_high": math.log2(OR) + 1.96 * se, "p_value": fisher_two_sided(aa, bb, cc, dd), "candidate_orthogroups": len(candidates), "candidate_present": aa, "candidate_absent": bb, "background_orthogroups": len(background), "background_present": cc, "background_absent": dd})
    order = sorted(range(len(out)), key=lambda i: out[i]["p_value"]); q = [1.0] * len(out)
    for rank, idx in enumerate(order, 1): q[idx] = min(1.0, out[idx]["p_value"] * len(out) / rank)
    for i in range(len(order) - 2, -1, -1): q[order[i]] = min(q[order[i]], q[order[i + 1]])
    for row, qv in zip(out, q): row["q_value"] = qv
    return out

def sequence_id_map():
    out = {}
    with (OF / "WorkingDirectory" / "SequenceIDs.txt").open(encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if ":" in line: k, v = line.rstrip().split(":", 1); out[k.strip()] = v.strip()
    return out

def annotated_tree(og, seqmap):
    path = OF / "WorkingDirectory" / "Trees_ids" / f"{og}.txt"
    if not path.exists() or path.stat().st_size == 0: return ""
    s = path.read_text(encoding="utf-8", errors="ignore").strip()
    return re.sub(r"(?<![A-Za-z0-9_])(?:\d+_\d+)(?=[:),])", lambda m: seqmap.get(m.group(0), m.group(0)), s)

def synteny_rows(candidate_rows):
    """Extract exact candidate-containing JCVI anchors across the four phase pairs."""
    wanted = {r.get("transcript_id", ""): r for r in candidate_rows if r.get("transcript_id")}
    out = []
    pairs = [("A1", "A2"), ("A1", "B1"), ("A2", "B2"), ("B1", "B2")]
    for left, right in pairs:
        folder = JCVI_ROOT / f"{left}_{right}"
        path = folder / f"hap{left}.hap{right}.anchors"
        if not path.exists(): continue
        with path.open(encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"): continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2: continue
                a, b = parts[0], parts[1]
                if a not in wanted and b not in wanted: continue
                for candidate, partner in ((a, b), (b, a)):
                    if candidate not in wanted: continue
                    hit = wanted[candidate]
                    out.append({"family": hit.get("family", ""), "candidate_phase": hit.get("phase", ""), "candidate_gene": candidate, "partner_gene": partner, "phase_pair": f"{left}_{right}", "anchor_score": parts[2] if len(parts) > 2 else "", "anchor_file": str(path.relative_to(ROOT))})
    return out

def write_tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows: path.write_text("", encoding="utf-8"); return
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)

def write_qc_manifest(records, summaries, effects, syn, outdir):
    high = re.compile(r"\b(LEAFY|LFY|AP1|APETALA[- ]?1|FRUITFULL|FUL|TFL1|CENTRORADIALIS|CEN|FLOWERING LOCUS T|SOC1|SHORT VEGETATIVE PHASE|SVP|CLAVATA[123]?|BRANCHED[- ]?1|LAX PANICLE)\b", re.I)
    manifest = []
    for row in records:
        tier = "high" if high.search(row["note"]) else "moderate"
        if row["mapping_status"] == "annotation_only": tier = "annotation_only_review"
        manifest.append({"family": row["family"], "phase": row["phase"], "gene_id": row["gene_id"], "transcript_id": row["transcript_id"], "old_gene_id": row["old_gene_id"], "stable_id": row["stable_id"], "orthogroup": row["orthogroup"], "mapping_score": row["mapping_score"], "review_tier": tier, "note": row["note"]})
    write_tsv(outdir / "candidate_gene_manifest.tsv", manifest)
    by_status = defaultdict(int)
    for row in records: by_status[(row["phase"], row["mapping_status"])] += 1
    summary = [{"metric": "candidate_records", "value": len(records)}, {"metric": "mapped_records", "value": sum(1 for r in records if r["mapping_status"] == "mapped_to_orthofinder")}, {"metric": "orthogroup_effect_rows", "value": len(effects)}, {"metric": "jcvi_anchor_endpoint_hits", "value": len(syn)}, {"metric": "unique_candidate_families", "value": len(summaries)}]
    for (phase, status), n in sorted(by_status.items()): summary.append({"metric": f"{phase}_{status}", "value": n})
    write_tsv(outdir / "analysis_qc_summary.tsv", summary)

def draw_figures(records, summaries, effects, trees, syn, outdir):
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    import numpy as np
    import seaborn as sns
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 7, "axes.spines.right": False, "axes.spines.top": False, "axes.linewidth": 0.8})
    sns.set_theme(style="white", context="paper")
    fig = plt.figure(figsize=(11.5, 8.3), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.08, 0.92], height_ratios=[1.0, 1.15])
    ax = fig.add_subplot(gs[0, 0]); fams = [r["family"] for r in summaries]
    mat = np.array([[r[p] for p in PHASES] for r in summaries], dtype=float) if summaries else np.zeros((1, 4))
    sns.heatmap(mat, cmap=sns.color_palette(["#f2f2f2", "#3c78a8"], as_cmap=True), vmin=0, vmax=1, cbar=False, linewidths=0.4, linecolor="white", annot=mat.astype(int), fmt="d", xticklabels=PHASES, yticklabels=fams or ["none"], ax=ax)
    ax.set_title("5c  Candidate-family retention", loc="left", fontweight="bold"); ax.set_xlabel("Operational PZI phase track"); ax.set_ylabel("Function-driven family")
    ax = fig.add_subplot(gs[0, 1]); plot_rows = [e for e in effects if e["comparison"] == "A_vs_background"] + [e for e in effects if e["comparison"] == "B_vs_background"]; ys = np.arange(len(plot_rows))
    for i, e in enumerate(plot_rows):
        color = "#b45f06" if e["comparison"].startswith("A") else "#3c78a8"; ax.plot([e["ci_low"], e["ci_high"]], [i, i], color=color, lw=1.2); ax.scatter(e["log2_effect"], i, color=color, s=18, zorder=3)
    ax.axvline(0, color="#777777", lw=0.8, ls="--"); ax.set_yticks(ys); ax.set_yticklabels([f"{e['family']} ({e['comparison'][0]})" for e in plot_rows], fontsize=6); ax.set_xlabel("log2 odds ratio vs. all PZI-containing OrthoFinder families"); ax.set_title("5d  Candidate-family effect sizes", loc="left", fontweight="bold")

    class Node:
        def __init__(self, name="", length=0.0): self.name, self.length, self.children = name, length, []
    def parse_newick(s):
        i = 0
        def sub():
            nonlocal i
            n = Node()
            if i < len(s) and s[i] == "(":
                i += 1
                while True:
                    n.children.append(sub())
                    if s[i] == ",": i += 1; continue
                    if s[i] == ")": i += 1; break
            st = i
            while i < len(s) and s[i] not in ":,();": i += 1
            n.name = s[st:i]
            if i < len(s) and s[i] == ":":
                i += 1; st = i
                while i < len(s) and s[i] not in ",();": i += 1
                try: n.length = float(s[st:i])
                except ValueError: n.length = 0.0
            return n
        return sub()
    def draw_tree(ax, newick, title):
        if not newick: ax.text(.5, .5, "No usable OrthoFinder tree", ha="center", va="center"); ax.axis("off"); return
        root = parse_newick(newick); leaves = []
        def collect(n, depth=0):
            if not n.children: leaves.append((n, depth)); return
            for c in n.children: collect(c, depth + c.length)
        collect(root); leaves = leaves[-80:]; ymap = {id(n): i for i, (n, _) in enumerate(leaves)}
        def xy(n, x=0):
            if not n.children: return x + n.length, ymap.get(id(n), 0)
            coords = [xy(c, x + n.length) for c in n.children]; yy = [q[1] for q in coords]; xx = [q[0] for q in coords]; mid = sum(yy) / len(yy)
            ax.plot([min(xx), max(xx)], [mid, mid], color="#444444", lw=.45)
            for cx, cy in coords: ax.plot([cx, cx], [cy, mid], color="#444444", lw=.45)
            return x + n.length, mid
        xy(root)
        for n, _ in leaves:
            yy = ymap[id(n)]; color = "#b45f06" if "pzi_A__" in n.name else "#3c78a8" if "pzi_B__" in n.name else "#666666"; ax.text(0, yy, n.name.replace("__", "|"), fontsize=4.5, color=color, va="center", ha="right", transform=ax.get_yaxis_transform())
        ax.set_title(title, loc="left", fontweight="bold", fontsize=7); ax.set_ylim(-1, max(1, len(leaves))); ax.set_yticks([]); ax.set_xticks([])
    ax = fig.add_subplot(gs[1, 0]); ax.axis("off"); inner = gs[1, 0].subgridspec(max(1, len(trees)), 1)
    for i, t in enumerate(trees): draw_tree(fig.add_subplot(inner[i, 0]), t["newick"], f"5e  {t['family']} | {t['orthogroup']} | {t['state']}")
    ax = fig.add_subplot(gs[1, 1]); ax.set_title("5f  Local synteny and copy fate", loc="left", fontweight="bold"); ax.set_xlim(0, 1.15); ax.set_ylim(0, max(1, len(summaries)) * 1.2); ax.axis("off")
    anchor_counts = defaultdict(int)
    for hit in syn:
        anchor_counts[(hit["family"], hit["candidate_phase"])] += 1
    preferred = [s for fam in F5_DISPLAY for s in summaries if s["family"] == fam and any(anchor_counts[(fam, p)] for p in PHASES)]
    display = preferred or [s for s in summaries if any(anchor_counts[(s["family"], p)] for p in PHASES)] or summaries
    for i, s in enumerate(display):
        y = len(display) * 1.2 - i * 1.2 - .8; ax.text(0, y + .25, s["family"], fontsize=6, va="center")
        for j, p in enumerate(PHASES):
            x = .34 + j * .13; present = s[p]; count = anchor_counts[(s["family"], p)]
            color = ("#b45f06" if p.startswith("A") else "#3c78a8") if present else "#eeeeee"
            ax.add_patch(Rectangle((x, y), .09, .34, facecolor=color, edgecolor="#444444", lw=.4))
            ax.text(x + .045, y + .17, str(count) if count else ("-" if not present else "0"), ha="center", va="center", fontsize=5, color="white" if present else "#555555")
        ax.text(1.04, y + .17, s["state"], fontsize=5.5, va="center", ha="left")
    ax.text(.385, len(display) * 1.2 + .05, "A1", fontsize=5, ha="center"); ax.text(.515, len(display) * 1.2 + .05, "A2", fontsize=5, ha="center"); ax.text(.645, len(display) * 1.2 + .05, "B1", fontsize=5, ha="center"); ax.text(.775, len(display) * 1.2 + .05, "B2", fontsize=5, ha="center")
    ax.text(0, .12, "Numbers: exact JCVI anchor hits for candidate transcripts; 0 = annotated candidate without a retained local anchor.", fontsize=5.2)
    fig.suptitle("Inflorescence-related gene retention and copy fate in operational PZI A/B tracks", x=.01, ha="left", fontsize=10, fontweight="bold")
    for ext, dpi in [("svg", None), ("pdf", None), ("png", 600), ("tiff", 600)]:
        kw = {"bbox_inches": "tight"};
        if dpi: kw["dpi"] = dpi
        fig.savefig(outdir / f"figure5c_5f_inflorescence.{ext}", **kw)
    plt.close(fig)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--outdir", type=Path, default=OUT); ap.add_argument("--skip-plot", action="store_true"); args = ap.parse_args()
    genes = {p: parse_gff(p) for p in PHASES}; map_canonical_to_stable(genes); gene_to_og, og_members = load_orthogroups(); records, summaries, bg = candidate_tables(genes, gene_to_og, og_members); effects = effect_sizes(summaries, bg); syn = synteny_rows(records); seqmap = sequence_id_map()
    trees = []
    ordered_summaries = [s for fam in TREE_PRIORITY for s in summaries if s["family"] == fam] + [s for s in summaries if s["family"] not in TREE_PRIORITY]
    for s in ordered_summaries:
        for og in s["orthogroups"].split(";") if s["orthogroups"] else []:
            nwk = annotated_tree(og, seqmap)
            if nwk: trees.append({"family": s["family"], "orthogroup": og, "state": s["state"], "newick": nwk}); break
        if len(trees) >= 3: break
    args.outdir.mkdir(parents=True, exist_ok=True); write_tsv(args.outdir / "candidate_gene_records.tsv", records); write_tsv(args.outdir / "candidate_family_summary.tsv", summaries); write_tsv(args.outdir / "candidate_effect_sizes.tsv", effects); write_tsv(args.outdir / "local_synteny_candidates.tsv", syn); write_tsv(args.outdir / "representative_gene_trees.tsv", trees); write_qc_manifest(records, summaries, effects, syn, args.outdir)
    (args.outdir / "README.md").write_text("# Fig. 5c-f data package\n\nPython-generated. A/B are operational PZI tracks; A1/A2 and B1/B2 are phase copies. Candidates come from explicit function descriptions in canonical PZI GFF3 Notes and should be manually reviewed before publication. Exact protein matching or coordinate overlap bridges canonical LOC annotations to the frozen OrthoFinder IDs; copy-2 records without a stable ID remain annotation-only. Effect sizes compare candidate orthogroup sets with all PZI-containing OrthoFinder families; this is not a GO enrichment test. 5f contains exact JCVI anchor hits across A1-A2, A1-B1, A2-B2 and B1-B2; counts are candidate endpoint hits and can include both endpoints of one anchor. Empty heatmap cells are not called losses. Direct parental species are not inferred.\n", encoding="utf-8")
    (args.outdir / "figure5_contract.md").write_text("# Figure contract\n\nCore conclusion: inflorescence/meristem-associated families show phase-resolved retention patterns linked to copy fate without assigning parental species.\n\n5c: A1/A2/B1/B2 presence heatmap. 5d: orthogroup-level log2 odds ratios versus the non-candidate PZI-containing OrthoFinder background with Fisher tests and BH q-values. 5e: representative OrthoFinder trees with PZI tracks highlighted. 5f: exact JCVI anchor counts for candidate transcripts across four phase pairs, paired with family copy-fate states.\n\nReview risks: ascertainment, annotation-note ambiguity, incomplete phase mapping, coordinate-bridge uncertainty, and lack of direct parental evidence.\n", encoding="utf-8")
    if not args.skip_plot:
        try: draw_figures(records, summaries, effects, trees, syn, args.outdir)
        except ImportError as exc: print(f"Data tables written; plotting dependency unavailable: {exc}")
    print(f"Wrote {args.outdir}")

if __name__ == "__main__": main()
