#!/usr/bin/env python3
"""Cross-species inflorescence-family pipeline (PZI + PSP + PPA).

Evidence-first screen across phased haplotype protein sets:
  PZI  A1/A2/B1/B2 (0903/haploid_repartition/inputs/pzi/pep/pzi_hap*.pep.fa)
  PSP  A1/A2/B1/B2 (0903/haploid_repartition/inputs/psp/pep/psp_hap*.pep.fa)
  PPA  A/B         (0903/haploid_repartition/inputs/ppa/pep/ppa_hap*.pep.fa)

Steps: HMMER domain screen -> reference-anchored ML trees -> placement ->
tree-supported candidate records -> homeolog anchor components (A1-B1/A2-B2
for PZI/PSP, A-B for PPA) -> per-species copy counts and exact tests.

PZI uses the revised 0903 haplotype packages throughout (the June canonical
tracks are superseded).  PZI anchor files carry the same LOC-style identifiers
as the frozen packages, so the sequence bridge reduces to an identity map.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
import subprocess
import sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import run_inflorescence_family_pipeline as rifp  # noqa: E402

DEFAULT_OUT = ROOT / "analysis" / "inflorescence_family_pipeline_crossspecies_20260912"

SPECIES_PREFIX = ("pzi_hap", "psp_hap", "ppa_hap", "pcam_hap", "ppe_hap", "gtr_hap")
PHASE_RE = re.compile(r"^(pzi|psp|ppa|pcam|ppe|gtr)_hap(A1|A2|B1|B2|A|B|S)_")
LOC_RE = re.compile(r"(LOC_\d+)")
TRACK = {"A1": "A", "A2": "A", "B1": "B", "B2": "B", "A": "A", "B": "B", "S": "S"}
ALL_SPECIES = ["pzi", "psp", "ppa", "pcam", "ppe", "gtr"]
SPECIES_LABEL = {"pzi": "PZI", "psp": "PSP", "ppa": "PPA"}

TARGET_FILES = {
    "pzi": {
        "A1": ROOT / "0903" / "haploid_repartition" / "inputs" / "pzi" / "pep" / "pzi_hapA1.pep.fa",
        "A2": ROOT / "0903" / "haploid_repartition" / "inputs" / "pzi" / "pep" / "pzi_hapA2.pep.fa",
        "B1": ROOT / "0903" / "haploid_repartition" / "inputs" / "pzi" / "pep" / "pzi_hapB1.pep.fa",
        "B2": ROOT / "0903" / "haploid_repartition" / "inputs" / "pzi" / "pep" / "pzi_hapB2.pep.fa",
    },
    "psp": {
        "A1": ROOT / "0903" / "haploid_repartition" / "inputs" / "psp" / "pep" / "psp_hapA1.pep.fa",
        "A2": ROOT / "0903" / "haploid_repartition" / "inputs" / "psp" / "pep" / "psp_hapA2.pep.fa",
        "B1": ROOT / "0903" / "haploid_repartition" / "inputs" / "psp" / "pep" / "psp_hapB1.pep.fa",
        "B2": ROOT / "0903" / "haploid_repartition" / "inputs" / "psp" / "pep" / "psp_hapB2.pep.fa",
    },
    "ppa": {
        "A": ROOT / "0903" / "haploid_repartition" / "inputs" / "ppa" / "pep" / "ppa_hapA.pep.fa",
        "B": ROOT / "0903" / "haploid_repartition" / "inputs" / "ppa" / "pep" / "ppa_hapB.pep.fa",
    },
    "pcam": {
        "S": ROOT / "analysis" / "canonical" / "pcam" / "primary.pep.fa",
    },
    "ppe": {
        "S": ROOT / "analysis" / "canonical" / "ppe" / "primary.pep.fa",
    },
    "gtr": {
        "S": ROOT / "0903" / "karyotype" / "prepared" / "gtr" / "gtr.pep.fa",
    },
}

LEGACY_PZI_PEP = {
    "A1": ROOT / "0903" / "haploid_repartition" / "inputs" / "pzi" / "pep" / "pzi_hapA1.pep.fa",
    "A2": ROOT / "0903" / "haploid_repartition" / "inputs" / "pzi" / "pep" / "pzi_hapA2.pep.fa",
    "B1": ROOT / "0903" / "haploid_repartition" / "inputs" / "pzi" / "pep" / "pzi_hapB1.pep.fa",
    "B2": ROOT / "0903" / "haploid_repartition" / "inputs" / "pzi" / "pep" / "pzi_hapB2.pep.fa",
}

ANCHOR_FILES = {
    "pzi": {
        ("A1", "B1"): ROOT / "0903/downstream_analysis/figure2_repartition_ks/01_jcvi_pairs/work/pzi_hapA1_vs_pzi_hapB1/pzi_hapA1.pzi_hapB1.anchors",
        ("A2", "B2"): ROOT / "0903/downstream_analysis/figure2_repartition_ks/01_jcvi_pairs/work/pzi_hapA2_vs_pzi_hapB2/pzi_hapA2.pzi_hapB2.anchors",
    },
    "psp": {
        ("A1", "B1"): ROOT / "0903/downstream_analysis/figure2_repartition_ks/01_jcvi_pairs/work/psp_hapA1_vs_psp_hapB1/psp_hapA1.psp_hapB1.anchors",
        ("A2", "B2"): ROOT / "0903/downstream_analysis/figure2_repartition_ks/01_jcvi_pairs/work/psp_hapA2_vs_psp_hapB2/psp_hapA2.psp_hapB2.anchors",
    },
    "ppa": {
        ("A", "B"): ROOT / "0903/downstream_analysis/figure2_repartition_ks/01_jcvi_pairs/work/ppa_hapA_vs_ppa_hapB/ppa_hapA.ppa_hapB.anchors",
    },
}

CONFIRM_BOOT = 80.0
CONFIRM_ALRT = 80.0
FH_WORKERS = 4


def log(msg: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as fh:
        fh.write(msg + "\n")
    print(msg, flush=True)


def is_target_id(key: str) -> bool:
    return key.startswith(SPECIES_PREFIX)


def gene_id_from_key(key: str) -> str:
    m = LOC_RE.search(key)
    if m:
        return m.group(1)
    return re.sub(r"[^A-Za-z0-9_.-]", "_", key)


def build_targets(out: Path, logfile: Path, species_list: list[str]):
    """Concatenate phased/single-track proteins (all isoforms) + references."""
    combined: dict[str, tuple[str, str]] = {}
    phase_by_id: dict[str, str] = {}
    species_by_id: dict[str, str] = {}
    loc_by_id: dict[str, str] = {}
    for species in species_list:
        phases = TARGET_FILES[species]
        for phase, path in phases.items():
            if not path.exists():
                raise FileNotFoundError(path)
            for key, (header, seq) in rifp.read_fasta(path).items():
                gid = gene_id_from_key(key)
                new_key = key if key.startswith(f"{species}_hap") else f"{species}_hap{phase}_{re.sub(r'[^A-Za-z0-9_.-]', '_', key)}"
                if new_key in combined:
                    raise RuntimeError(f"duplicate target id: {new_key}")
                combined[new_key] = (new_key, seq)
                phase_by_id[new_key] = phase
                species_by_id[new_key] = species
                loc_by_id[new_key] = gid
    at_fa = rifp.REF_DIR / "Arabidopsis_thaliana.TAIR10.pep.all.fa.gz"
    at_records: dict[str, tuple[str, str]] = {}
    with gzip.open(at_fa, "rt", encoding="utf-8", errors="ignore") as fh:
        key = None
        header = ""
        chunks: list[str] = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if key is not None:
                    at_records[key] = (header, "".join(chunks).replace("*", ""))
                header = line[1:].strip()
                key = header.split()[0]
                chunks = []
            elif key is not None:
                chunks.append(line.strip())
        if key is not None:
            at_records[key] = (header, "".join(chunks).replace("*", ""))
    combined.update(at_records)
    rosaceae_meta, rosaceae_records = rifp.load_rosaceae_references()
    if not rosaceae_meta or not rosaceae_records:
        raise FileNotFoundError("Rosaceae references missing; run prepare_rosaceae_references.py")
    overlap = set(combined).intersection(rosaceae_records)
    if overlap:
        raise RuntimeError(f"reference/target ID overlap: {sorted(overlap)[:3]}")
    combined.update(rosaceae_records)
    target_fa = out / "targets" / "all_species.primary.pep.fa"
    rifp.write_fasta(target_fa, combined)
    n_targets = sum(1 for k in combined if is_target_id(k))
    log(f"[targets] {n_targets} phased proteins + {len(combined) - n_targets} references -> {target_fa}", logfile)
    return phase_by_id, species_by_id, loc_by_id, at_records, rosaceae_meta


def complete_targets(family: str, family_hits) -> list[str]:
    requirements = rifp.FAMILY_REQUIREMENTS[family]
    if len(requirements) == 1:
        profile = requirements[0]
        return sorted(
            x for x, hit in family_hits[profile].items()
            if is_target_id(x) and float(hit["profile_coverage"]) >= 0.5
        )
    common = set(family_hits[requirements[0]])
    for profile in requirements[1:]:
        common &= set(family_hits[profile])
    return sorted(
        x for x in common
        if is_target_id(x)
        and all(float(family_hits[p][x]["profile_coverage"]) >= 0.5 for p in requirements)
    )


def run_hmmsearch(out: Path, combined_fa: Path, logfile: Path, threads: int):
    rifp.extract_profiles(rifp.HMM_DB, out / "reference" / "selected_profiles.hmm")

    def one(profile: str):
        hmm = out / "reference" / f"{profile}.hmm"
        domtbl = out / "hmmsearch" / f"{profile}.domtblout"
        domtbl.parent.mkdir(parents=True, exist_ok=True)
        cpu = max(2, threads // FH_WORKERS)
        cmd = [rifp.HMMSEARCH, "--cut_ga", "--cpu", str(cpu), "--domtblout", str(domtbl), "--noali", str(hmm), str(combined_fa)]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return profile, domtbl

    family_hits: dict[str, dict[str, dict[str, object]]] = {}
    with ThreadPoolExecutor(max_workers=FH_WORKERS) as pool:
        for profile, domtbl in pool.map(one, list(rifp.PROFILES)):
            rows = rifp.parse_domtbl(domtbl, profile, {})
            family_hits[profile] = rifp.best_hits(rows)
            log(f"[hmmsearch] {profile}: {len(family_hits[profile])} best hits", logfile)
    return family_hits


def select_tree_records(records: list[str], combined) -> list[str]:
    """Longest isoform per (species, phase, LOC); references pass through."""
    best: dict[str, str] = {}
    for key in records:
        if not is_target_id(key):
            best[key] = key
            continue
        m = PHASE_RE.match(key)
        gene = f"{m.group(1)}_hap{m.group(2)}_{gene_id_from_key(key)}"
        old = best.get(gene)
        if old is None or len(combined[key][1]) > len(combined[old][1]):
            best[gene] = key
    return sorted(best.values())


def make_tree_inputs(out: Path, family_hits, combined, reference_metadata, logfile: Path):
    tree_targets: dict[str, list[str]] = {}
    summary = []
    for family in rifp.FAMILY_REQUIREMENTS:
        complete = complete_targets(family, family_hits)
        records = select_tree_records(complete, combined)
        refs = rifp.reference_ids_for_family(family, reference_metadata)
        ids = sorted(set(records) | set(refs))
        tree_targets[family] = ids
        rifp.write_fasta(out / "sequences" / f"{family}.fa", {x: combined[x] for x in ids if x in combined})
        n_target = sum(1 for x in records)
        summary.append({"family": family, "complete_records": len(complete), "tree_targets": n_target, "total_with_refs": len(ids)})
        log(f"[tree-input] {family}: complete={len(complete)} gene-level={n_target} total={len(ids)}", logfile)
    rifp.write_rows(out / "sequences" / "tree_input_summary.tsv", summary)
    return tree_targets


def run_align_trees(out: Path, tree_targets, logfile: Path, threads: int, only: list[str] | None = None) -> None:
    def one(family: str):
        if only and family not in only:
            return
        ids = tree_targets[family]
        if len(ids) < 4:
            return
        inp = out / "sequences" / f"{family}.fa"
        aln = out / "alignments" / f"{family}.mafft.fa"
        trimmed = out / "alignments" / f"{family}.trimmed.fa"
        prefix = out / "trees" / family
        aln.parent.mkdir(parents=True, exist_ok=True)
        prefix.parent.mkdir(parents=True, exist_ok=True)
        treefile = out / "trees" / f"{family}.treefile"
        if not treefile.exists() or treefile.stat().st_size == 0:
            with aln.open("w", encoding="utf-8") as fh:
                subprocess.run([rifp.MAFFT, "--auto", "--thread", "4", str(inp)], check=True, stdout=fh, stderr=subprocess.PIPE, text=True)
            try:
                rifp.run([rifp.TRIMAL, "-in", str(aln), "-out", str(trimmed), "-automated1"], logfile)
            except subprocess.CalledProcessError:
                import shutil
                shutil.copyfile(aln, trimmed)
            model = rifp.IQTREE_FIXED_MODELS.get(family, "MFP")
            cpu = max(2, threads // FH_WORKERS)
            rifp.run([rifp.IQTREE, "-s", str(trimmed), "-m", model, "-B", "1000", "--alrt", "1000", "--bnni",
                      "-T", str(cpu), "-seed", rifp.IQTREE_SEED, "--prefix", str(prefix)], logfile)
        log(f"[tree-done] {family}", logfile)

    with ThreadPoolExecutor(max_workers=FH_WORKERS) as pool:
        list(pool.map(one, list(tree_targets)))


def tree_placements(out: Path, tree_targets, reference_metadata, loc_by_id, logfile: Path):
    rows: list[dict[str, object]] = []
    for family, ids in tree_targets.items():
        treefile = out / "trees" / f"{family}.treefile"
        if not treefile.exists() or treefile.stat().st_size == 0:
            log(f"[placement] {family}: no tree", logfile)
            continue
        root = rifp.parse_newick(treefile.read_text(encoding="utf-8", errors="ignore").strip().rstrip(";"))
        leaves: dict[str, object] = {}

        def walk(node, distance):
            node.root_distance = distance
            if not node.children:
                leaves[node.name] = node
                return
            for child in node.children:
                walk(child, distance + child.length)

        walk(root, 0.0)
        refs = [x for x in leaves if x in reference_metadata and reference_metadata[x].get("tree_family") == family]
        targets = [x for x in leaves if is_target_id(x)]
        for target in targets:
            m = PHASE_RE.match(target)
            species = m.group(1)
            phase = m.group(2)
            gene = f"{species}_hap{phase}_{loc_by_id.get(target, '')}"
            if not refs:
                rows.append({"family": family, "target_id": target, "gene_key": gene, "species": species, "phase": phase,
                             "nearest_reference": "", "nearest_reference_symbol": "", "nearest_reference_group": "",
                             "patristic_distance": "", "mrca_bootstrap": "", "mrca_sh_alrt": "",
                             "reference_group_count": 0, "group_support_count": 0, "subfamily_group": "",
                             "placement_status": "no_reference_in_tree"})
                continue
            a = leaves[target]
            ancestors = {}
            n = a
            while n is not None:
                ancestors[id(n)] = n
                n = n.parent
            scored = []
            for ref in refs:
                b = leaves[ref]
                n = b
                while n is not None and id(n) not in ancestors:
                    n = n.parent
                if n is None:
                    continue
                dist = a.root_distance + b.root_distance - 2 * n.root_distance
                scored.append((dist, ref, n))
            if not scored:
                continue
            dist, ref, mrca = min(scored, key=lambda x: x[0])
            boot, alrt = rifp.support_values(mrca)
            ref_meta = reference_metadata[ref]
            group = ref_meta.get("subfamily_group", "")
            group_refs = [r for r in refs if reference_metadata[r].get("subfamily_group") == group]
            group_support_count = 0
            for group_ref in group_refs:
                b = leaves[group_ref]
                n = b
                while n is not None and id(n) not in ancestors:
                    n = n.parent
                if n is None:
                    continue
                gb, ga = rifp.support_values(n)
                if gb is not None and ga is not None and gb >= CONFIRM_BOOT and ga >= CONFIRM_ALRT:
                    group_support_count += 1
            supported = boot is not None and alrt is not None and boot >= CONFIRM_BOOT and alrt >= CONFIRM_ALRT
            rows.append({"family": family, "target_id": target, "gene_key": gene, "species": species, "phase": phase,
                         "nearest_reference": ref, "nearest_reference_symbol": ref_meta.get("gene_symbol", ""),
                         "nearest_reference_group": group, "patristic_distance": f"{dist:.6g}",
                         "mrca_bootstrap": "" if boot is None else f"{boot:.1f}",
                         "mrca_sh_alrt": "" if alrt is None else f"{alrt:.1f}",
                         "reference_group_count": len(group_refs), "group_support_count": group_support_count,
                         "subfamily_group": group,
                         "placement_status": "supported" if supported else "weak"})
    rifp.write_rows(out / "trees" / "tree_placements_crossspecies.tsv", rows)
    return rows


def collapse_gene_evidence(out: Path, placements, loc_by_id):
    """Best placement per (family, gene_key); write gene-level evidence table."""
    best: dict[tuple[str, str], dict[str, object]] = {}
    for p in placements:
        key = (str(p["family"]), str(p["gene_key"]))
        old = best.get(key)
        rank = (p.get("placement_status") == "supported", float(p["mrca_bootstrap"] or 0), float(p["mrca_sh_alrt"] or 0))
        if old is None or rank > old["_rank"]:
            row = dict(p)
            row["_rank"] = rank
            best[key] = row
    rows = []
    for (family, gene), p in sorted(best.items()):
        boot = float(p["mrca_bootstrap"]) if p.get("mrca_bootstrap") else None
        alrt = float(p["mrca_sh_alrt"]) if p.get("mrca_sh_alrt") else None
        gc = int(p.get("reference_group_count", 0) or 0)
        gs = int(p.get("group_support_count", 0) or 0)
        strict = boot is not None and alrt is not None and boot >= CONFIRM_BOOT and alrt >= CONFIRM_ALRT
        supported = bool(p.get("subfamily_group")) and strict and gc >= 2 and gs >= 2
        rows.append({
            "family": family, "gene_key": gene, "record_id": p["target_id"],
            "species": p["species"], "phase": p["phase"], "track": TRACK.get(str(p["phase"]), ""),
            "loc": loc_by_id.get(str(p["target_id"]), ""),
            "nearest_reference_symbol": p.get("nearest_reference_symbol", ""),
            "subfamily_group": p.get("subfamily_group", ""),
            "mrca_bootstrap": p.get("mrca_bootstrap", ""), "mrca_sh_alrt": p.get("mrca_sh_alrt", ""),
            "reference_group_count": gc, "group_support_count": gs,
            "placement_status": p.get("placement_status", ""),
            "evidence_tier": "tree_supported" if supported else ("tree_weak" if p.get("subfamily_group") else "domain_only"),
        })
    rifp.write_rows(out / "candidates_crossspecies.tsv", rows)
    return rows


def pzi_legacy_bridge(combined, logfile: Path):
    seq_to_legacy: dict[str, str] = {}
    for phase, path in LEGACY_PZI_PEP.items():
        for key, (_h, seq) in rifp.read_fasta(path).items():
            seq_to_legacy.setdefault(seq, key)
    forward: dict[str, str] = {}
    for key, (_h, seq) in combined.items():
        if key.startswith("pzi_hap"):
            legacy = seq_to_legacy.get(seq)
            if legacy:
                forward[key] = legacy
    reverse = {v: k for k, v in forward.items()}
    log(f"[bridge] canonical PZI -> legacy ids: {len(forward)}", logfile)
    return forward, reverse


def build_anchor_edges(out: Path, combined, legacy_reverse, logfile: Path):
    rows: list[dict[str, object]] = []
    for species, pairs in ANCHOR_FILES.items():
        for (left, right), path in pairs.items():
            if not path.exists():
                log(f"[anchor] missing {path}", logfile)
                continue
            kept = 0
            with path.open(encoding="utf-8", errors="ignore") as fh:
                for line in fh:
                    if not line.strip() or line.startswith("#"):
                        continue
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 2:
                        continue
                    a_raw, b_raw = f[0], f[1]
                    loc_a = LOC_RE.search(a_raw)
                    loc_b = LOC_RE.search(b_raw)
                    if not loc_a or not loc_b:
                        continue
                    if species == "pzi":
                        a = legacy_reverse.get(a_raw)
                        b = legacy_reverse.get(b_raw)
                        a_label = f"LOC={loc_a.group(1)}"
                        b_label = f"LOC={loc_b.group(1)}"
                        a_present = a is not None
                        b_present = b is not None
                        if not a_present and not b_present:
                            continue
                        rows.append({
                            "species": species, "phase_pair": f"{left}_{right}", "left_phase": left, "right_phase": right,
                            "left_id": a or a_label, "right_id": b or b_label,
                            "left_loc": loc_a.group(1), "right_loc": loc_b.group(1),
                            "left_in_target": int(a_present), "right_in_target": int(b_present),
                            "anchor_score": f[2] if len(f) > 2 else "",
                        })
                    else:
                        a = f"{species}_hap{left}_{a_raw}"
                        b = f"{species}_hap{right}_{b_raw}"
                        a_present = a in combined
                        b_present = b in combined
                        if not a_present and not b_present:
                            continue
                        rows.append({
                            "species": species, "phase_pair": f"{left}_{right}", "left_phase": left, "right_phase": right,
                            "left_id": a, "right_id": b, "left_loc": loc_a.group(1), "right_loc": loc_b.group(1),
                            "left_in_target": int(a_present), "right_in_target": int(b_present),
                            "anchor_score": f[2] if len(f) > 2 else "",
                        })
                    kept += 1
            log(f"[anchor] {species} {left}-{right}: {kept} edges kept", logfile)
    rifp.write_rows(out / "anchors" / "homeolog_anchor_edges.tsv", rows)
    return rows


def gene_key_of(record_id: str) -> str:
    return re.sub(r"-mRNA-\d+$", "", record_id)


def build_components(out: Path, edges, evidence, logfile: Path):
    supported = {r["gene_key"]: r for r in evidence if r["evidence_tier"] == "tree_supported"}
    by_species: dict[str, list[dict[str, object]]] = defaultdict(list)
    norm_edges = []
    for e in edges:
        e = dict(e)
        e["left_key"] = gene_key_of(str(e["left_id"])) if e["left_in_target"] else str(e["left_id"])
        e["right_key"] = gene_key_of(str(e["right_id"])) if e["right_in_target"] else str(e["right_id"])
        norm_edges.append(e)
        by_species[e["species"]].append(e)
    edges = norm_edges
    component_rows: list[dict[str, object]] = []
    for species, sedges in by_species.items():
        parent: dict[str, str] = {}

        def find(x: str) -> str:
            parent.setdefault(x, x)
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(a: str, b: str) -> None:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[rb] = ra

        for e in sedges:
            if e["left_in_target"] and e["right_in_target"]:
                union(str(e["left_key"]), str(e["right_key"]))
            elif e["left_in_target"]:
                parent.setdefault(str(e["left_key"]), str(e["left_key"]))
            elif e["right_in_target"]:
                parent.setdefault(str(e["right_key"]), str(e["right_key"]))
        groups: dict[str, list[str]] = defaultdict(list)
        for node in list(parent):
            groups[find(node)].append(node)
        for idx, (root, members) in enumerate(sorted(groups.items(), key=lambda kv: kv[0]), 1):
            sup = [m for m in members if m in supported]
            if not sup:
                continue
            fams = sorted({str(supported[m]["family"]) for m in sup})
            groups_sub = sorted({str(supported[m]["subfamily_group"]) for m in sup})
            phases = sorted({str(supported[m]["phase"]) for m in sup})
            has_a = any(TRACK.get(p) == "A" for p in phases)
            has_b = any(TRACK.get(p) == "B" for p in phases)
            fate = "A_and_B" if (has_a and has_b) else "A_only" if has_a else "B_only"
            pairs = sorted({str(e["phase_pair"]) for e in sedges if str(e["left_key"]) in members or str(e["right_key"]) in members})
            component_rows.append({
                "component_id": f"{species.upper()}_COMP_{idx:05d}", "species": species,
                "n_members": len(members), "n_tree_supported": len(sup),
                "families": ";".join(fams), "subfamily_groups": ";".join(groups_sub),
                "phase_pattern": "+".join(phases), "copy_fate": fate,
                "anchor_pairs": ";".join(pairs), "members": ";".join(members),
            })
    rifp.write_rows(out / "copy_fate" / "components_crossspecies.tsv", component_rows)
    log(f"[components] {len(component_rows)} components with tree-supported members", logfile)
    return component_rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--threads", type=int, default=24)
    ap.add_argument("--stage", default="all", choices=["all", "screen", "trees", "post"])
    ap.add_argument("--only", default="", help="comma-separated families for tree stage")
    ap.add_argument("--species", default="pzi,psp,ppa", help="comma-separated species subset")
    args = ap.parse_args()
    out = args.outdir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    logfile = out / "pipeline.log"
    species_list = [s for s in args.species.split(",") if s]
    log(f"=== crossspecies pipeline start stage={args.stage} species={','.join(species_list)} ===", logfile)

    phase_by_id, species_by_id, loc_by_id, at_records, rosaceae_meta = build_targets(out, logfile, species_list)
    combined = rifp.read_fasta(out / "targets" / "all_species.primary.pep.fa")

    family_hits = run_hmmsearch(out, out / "targets" / "all_species.primary.pep.fa", logfile, args.threads)

    at_symbols = {}
    with gzip.open(rifp.REF_DIR / "Arabidopsis_thaliana.TAIR10.pep.all.fa.gz", "rt", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                key = line[1:].split()[0]
                at_symbols[key] = rifp.parse_gene_symbol(line.strip())
    reference_metadata = rifp.build_reference_metadata(at_symbols, at_records, rosaceae_meta)

    tree_targets = make_tree_inputs(out, family_hits, combined, reference_metadata, logfile)
    if args.stage in {"all", "trees"}:
        only = [x for x in args.only.split(",") if x] or None
        run_align_trees(out, tree_targets, logfile, args.threads, only)
    if args.stage == "screen":
        log("=== stage screen done ===", logfile)
        return 0

    placements = tree_placements(out, tree_targets, reference_metadata, loc_by_id, logfile)
    evidence = collapse_gene_evidence(out, placements, loc_by_id)
    forward, reverse = pzi_legacy_bridge(combined, logfile)
    edges = build_anchor_edges(out, combined, reverse, logfile)
    build_components(out, edges, evidence, logfile)

    from crossspecies_stats import run_stats
    run_stats(out, evidence, edges, logfile)
    log("=== pipeline done ===", logfile)
    return 0


if __name__ == "__main__":
    sys.exit(main())
