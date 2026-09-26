#!/usr/bin/env python3
"""Evidence-first inflorescence-family workflow.

The workflow uses Pfam profiles for established inflorescence regulators.
K-box and SRF-TF are treated as complementary MADS-box domain evidence;
KNOX1 and KNOX2 are treated as complementary KNOX evidence.  Broad-domain
families are only accepted after reference-anchored phylogenetic placement.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "analysis" / "inflorescence_family_pipeline"
HMM_DB = ROOT / "software" / "Pfam-A.hmm.gz"
REF_DIR = DEFAULT_OUT / "reference"
ROS_REF_FA = REF_DIR / "rosaceae_reference.fa"
ROS_REF_META = REF_DIR / "rosaceae_reference_metadata.tsv"
HMMSEARCH = shutil.which("hmmsearch") or "/usr/bin/hmmsearch"
MAFFT = shutil.which("mafft") or "${CONDA_ROOT}/bin/mafft"
TRIMAL = shutil.which("trimal") or "${CONDA_ROOT}/bin/trimal"
IQTREE = shutil.which("iqtree") or "${CONDA_ROOT}/bin/iqtree"
IQTREE_SEED = "20260829"
# AP2 is a broad family; its model was already selected by BIC on the same
# trimmed alignment (Q.PFAM+G4).  Reusing that model avoids an expensive and
# unnecessary second all-model search while retaining full support resampling.
IQTREE_FIXED_MODELS = {
    "PEBP": "Q.PLANT+G4",
    "MADS": "Q.PLANT+I+G4",
    "C_LFY_FLO": "JTTDCMUT+G4",
    "AP2": "Q.PFAM+G4",
    "SPL": "Q.MAMMAL+I+R3",
    "TCP": "JTT+G4",
    "KNOX": "Q.PLANT+G4",
    "YABBY": "Q.BIRD+G4",
}

PHASES = ["A1", "A2", "B1", "B2"]
PROFILES = {
    "PEBP": "PF01161",
    "K_box": "PF01486",
    "SRF_TF": "PF00319",
    "C_LFY_FLO": "PF17538",
    "AP2": "PF00847",
    "SPL": "PF03110",
    "TCP": "PF03634",
    "KNOX1": "PF03790",
    "KNOX2": "PF03791",
    "Homeobox": "PF00046",
    "YABBY": "PF04690",
}

# A family enters its tree only when all required domains are complete.  The
# profile names are deliberately kept separate from tree-family labels so that
# composite-domain families can be audited without conflating HMM hits.
FAMILY_REQUIREMENTS = {
    "PEBP": ("PEBP",),
    "MADS": ("K_box", "SRF_TF"),
    "C_LFY_FLO": ("C_LFY_FLO",),
    "AP2": ("AP2",),
    "SPL": ("SPL",),
    "TCP": ("TCP",),
    "KNOX": ("KNOX1", "KNOX2"),
    # Homeobox (PF00046) is retained in PROFILES for a complete audit of
    # possible WUS/WOX hits, but is intentionally not a formal tree family:
    # the domain is too broad without a dedicated WOX reference panel and
    # would otherwise dominate the analysis with unrelated homeobox proteins.
    "YABBY": ("YABBY",),
}
TREE_FAMILIES = tuple(FAMILY_REQUIREMENTS)
PROFILE_TREE_FAMILIES = {
    profile: {family for family, reqs in FAMILY_REQUIREMENTS.items() if profile in reqs}
    for profile in PROFILES
}

# The two LFY/FLO records that had passed the automatic evidence filter were
# removed from the formal confirmed set at the user's request.  They remain
# in strict_candidates.tsv and final_candidate_evidence.tsv for auditability.
USER_EXCLUDED_TARGETS = {
    "pzi_hapA1_LOC_00021795-mRNA-1",
    "pzi_hapA2_LOC_00024365-mRNA-1",
}

# Confirmed Arabidopsis gene symbols used as anchors for subfamily placement.
REFERENCE_SYMBOLS = {
    "PEBP": {"FT", "TSF", "TFL1", "BFT", "ATC", "MFT"},
    "MADS": {
        "AP1", "CAL", "FUL", "AG", "AGL6", "AP3", "PI", "SEP1", "SEP2",
        "SEP3", "SEP4", "SVP", "SOC1", "AGL24", "AGL15", "AGL18", "AGL20",
        "AGL42", "AGL62", "AGL80", "SHP1", "SHP2", "STK", "ANR1", "FLC",
    },
    "LFY": {"LFY"},
    "AP2": {"AP2", "TOE1", "TOE2", "TOE3"},
    "SPL": {"SPL2", "SPL3", "SPL4", "SPL5", "SPL9", "SPL15"},
    "TCP": {"TCP2", "TCP3", "TCP4", "TCP10", "TCP24", "BRC1", "BRC2"},
    "KNOX": {"STM", "KNAT1", "KNAT2", "KNAT6"},
    "WOX": {"WUS", "WOX5", "WOX9"},
    "YABBY": {"YAB1", "YAB2", "YAB3", "YAB4", "YAB5"},
}


def load_rosaceae_references() -> tuple[dict[str, dict[str, str]], dict[str, tuple[str, str]]]:
    """Load curated Rosaceae reference sequences and metadata."""
    if not ROS_REF_FA.exists() or not ROS_REF_META.exists():
        return {}, {}
    metadata: dict[str, dict[str, str]] = {}
    with ROS_REF_META.open(encoding="utf-8", errors="ignore", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("protein_id"):
                metadata[row["protein_id"]] = dict(row)
    records = read_fasta(ROS_REF_FA)
    records = {key: value for key, value in records.items() if key in metadata}
    return metadata, records


def run(cmd: list[str], log: Path) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("a", encoding="utf-8") as fh:
        fh.write("$ " + " ".join(str(x) for x in cmd) + "\n")
        subprocess.run(cmd, check=True, stdout=fh, stderr=subprocess.STDOUT)


def read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    """Return id -> (header, sequence), keeping the first token as id."""
    out: dict[str, tuple[str, str]] = {}
    key = None
    header = ""
    seq: list[str] = []
    with path.open(encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if key is not None:
                    out[key] = (header, "".join(seq).replace("*", ""))
                header = line[1:].strip()
                key = header.split()[0]
                seq = []
            elif key is not None:
                seq.append(line.strip())
    if key is not None:
        out[key] = (header, "".join(seq).replace("*", ""))
    return out


def write_fasta(path: Path, records: dict[str, tuple[str, str]] | list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    items = records.items() if isinstance(records, dict) else records
    with path.open("w", encoding="utf-8") as fh:
        for key, value in items:
            if isinstance(value, tuple):
                header, seq = value
            else:
                header, seq = key, value
            fh.write(f">{header or key}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def parse_gene_symbol(header: str) -> str:
    m = re.search(r"(?:gene_symbol|gene_name):([^\s;]+)", header)
    return m.group(1) if m else ""


def extract_profiles(db_gz: Path, out_path: Path) -> dict[str, dict[str, str]]:
    """Extract requested Pfam records from a gzipped Pfam-A database."""
    wanted = set(PROFILES.values())
    found: dict[str, dict[str, str]] = {}
    records: dict[str, str] = {}
    current: list[str] = []
    with gzip.open(db_gz, "rt", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            current.append(line)
            if line.strip() != "//":
                continue
            record = "".join(current)
            current = []
            acc_m = re.search(r"^ACC\s+(\S+)", record, re.M)
            name_m = re.search(r"^NAME\s+(\S+)", record, re.M)
            if not acc_m or not name_m:
                continue
            acc = acc_m.group(1).split(".")[0]
            if acc in wanted:
                records[acc] = record
                found[acc] = {"name": name_m.group(1), "accession": acc_m.group(1)}
    if set(found) != wanted:
        missing = ",".join(sorted(wanted - set(found)))
        raise RuntimeError(f"Missing Pfam profiles: {missing}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("".join(records[a] for a in sorted(records)), encoding="utf-8")
    for family, acc in PROFILES.items():
        (out_path.parent / f"{family}.hmm").write_text(records[acc], encoding="utf-8")
    return found


def make_pzi_fasta(out: Path) -> dict[str, str]:
    records: dict[str, tuple[str, str]] = {}
    phase_by_id: dict[str, str] = {}
    for phase in PHASES:
        path = ROOT / "analysis" / "canonical" / f"pzi_hap{phase}" / "primary.pep.fa"
        if not path.exists():
            raise FileNotFoundError(path)
        for key, (header, seq) in read_fasta(path).items():
            if key in records:
                raise RuntimeError(f"Duplicate PZI protein ID: {key}")
            records[key] = (header, seq)
            phase_by_id[key] = phase
    write_fasta(out, records)
    return phase_by_id


def make_combined_target(
    pzi_fa: Path,
    at_fa_gz: Path,
    out: Path,
    external_records: dict[str, tuple[str, str]] | None = None,
) -> tuple[dict[str, str], dict[str, tuple[str, str]]]:
    pzi = read_fasta(pzi_fa)
    at_records: dict[str, tuple[str, str]] = {}
    with gzip.open(at_fa_gz, "rt", encoding="utf-8", errors="ignore") as fh:
        key = None
        header = ""
        seq: list[str] = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if key is not None:
                    at_records[key] = (header, "".join(seq).replace("*", ""))
                header = line[1:].strip()
                key = header.split()[0]
                seq = []
            elif key is not None:
                seq.append(line.strip())
        if key is not None:
            at_records[key] = (header, "".join(seq).replace("*", ""))
    combined = dict(pzi)
    combined.update(at_records)
    if external_records:
        overlap = set(combined).intersection(external_records)
        if overlap:
            raise RuntimeError(f"Reference IDs overlap target IDs: {sorted(overlap)[:3]}")
        combined.update(external_records)
    write_fasta(out, combined)
    symbols = {k: parse_gene_symbol(v[0]) for k, v in at_records.items()}
    return symbols, combined


def parse_domtbl(path: Path, profile: str, symbols: dict[str, str]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 23:
                continue
            # HMMER3 domtblout columns are fixed through accession/description.
            target, target_acc, tlen = f[0], f[1], int(f[2])
            query, query_acc, qlen = f[3], f[4], int(f[5])
            row = {
                "target_id": target,
                "target_accession": target_acc,
                "target_length": tlen,
                "profile": profile,
                "profile_name": query,
                "profile_accession": query_acc,
                "profile_length": qlen,
                "full_evalue": float(f[6]),
                "full_score": float(f[7]),
                "full_bias": float(f[8]),
                "domain_number": f"{f[9]}/{f[10]}",
                "domain_c_evalue": float(f[11]),
                "domain_i_evalue": float(f[12]),
                "domain_score": float(f[13]),
                "domain_bias": float(f[14]),
                "hmm_from": int(f[15]),
                "hmm_to": int(f[16]),
                "ali_from": int(f[17]),
                "ali_to": int(f[18]),
                "env_from": int(f[19]),
                "env_to": int(f[20]),
                "domain_accuracy": float(f[21]),
                "description": " ".join(f[22:]),
                "arabidopsis_symbol": symbols.get(target, ""),
            }
            row["profile_coverage"] = (row["hmm_to"] - row["hmm_from"] + 1) / qlen
            rows.append(row)
    return rows


def best_hits(rows: list[dict[str, object]]) -> dict[str, dict[str, object]]:
    best: dict[str, dict[str, object]] = {}
    for row in rows:
        key = str(row["target_id"])
        old = best.get(key)
        rank = (float(row["domain_i_evalue"]), -float(row["domain_score"]))
        if old is None or rank < (float(old["domain_i_evalue"]), -float(old["domain_score"])):
            best[key] = row
    return best


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def build_reference_metadata(
    at_symbols: dict[str, str],
    at_records: dict[str, tuple[str, str]],
    rosaceae_meta: dict[str, dict[str, str]],
) -> dict[str, dict[str, str]]:
    """Create a single metadata map for Arabidopsis and Rosaceae anchors."""
    metadata: dict[str, dict[str, str]] = {}
    for family, symbols in REFERENCE_SYMBOLS.items():
        tree_family = "C_LFY_FLO" if family == "LFY" else family
        # Keep the longest TAIR10 isoform for each curated gene symbol.
        chosen: dict[str, tuple[str, str]] = {}
        for key, symbol in at_symbols.items():
            if symbol not in symbols or key not in at_records:
                continue
            old = chosen.get(symbol)
            if old is None or len(at_records[key][1]) > len(old[1]):
                chosen[symbol] = (key, at_records[key][1])
        for symbol, (key, _seq) in chosen.items():
            if tree_family == "PEBP":
                group = "FT-like" if symbol in {"FT", "TSF", "BFT", "ATC"} else "TFL1-like" if symbol == "TFL1" else "MFT-like" if symbol == "MFT" else "PEBP-other"
            elif tree_family == "MADS":
                if symbol in {"AP1", "CAL"}: group = "AP1/CAL"
                elif symbol == "FUL": group = "FUL"
                elif symbol in {"AG", "SHP1", "SHP2", "STK"}: group = "AG/SHP/STK"
                elif symbol in {"SEP1", "SEP2", "SEP3", "SEP4"}: group = "SEP"
                elif symbol == "SVP": group = "SVP/DAM"
                elif symbol == "SOC1": group = "SOC1"
                elif symbol == "AGL6": group = "AGL6"
                else: group = "other-MADS"
            elif tree_family == "AP2":
                group = "AP2/TOE"
            elif tree_family == "SPL":
                group = "SPL3/4/5" if symbol in {"SPL3", "SPL4", "SPL5"} else "SPL9/15" if symbol in {"SPL9", "SPL15"} else "SPL-other"
            elif tree_family == "TCP":
                group = "BRC1/2" if symbol in {"BRC1", "BRC2"} else "TCP2/3/4" if symbol in {"TCP2", "TCP3", "TCP4"} else "TCP-other"
            elif tree_family == "KNOX":
                group = "STM" if symbol == "STM" else "KNAT1/2/6"
            elif tree_family == "WOX":
                group = "WUS" if symbol == "WUS" else "WOX5/9"
            elif tree_family == "YABBY":
                group = "YABBY1/2/3" if symbol in {"YAB1", "YAB2", "YAB3"} else "YABBY4/5"
            else:
                group = "LFY/FLO"
            metadata[key] = {
                "protein_id": key,
                "species": "Arabidopsis thaliana",
                "accession": key,
                "gene_symbol": symbol,
                "tree_family": tree_family,
                "subfamily_group": group,
                "curation_source": "TAIR10_gene_symbol_longest_isoform",
            }
    metadata.update(rosaceae_meta)
    return metadata


def reference_ids_for_family(family: str, reference_metadata: dict[str, dict[str, str]]) -> list[str]:
    return sorted(
        key for key, row in reference_metadata.items()
        if row.get("tree_family") == family
    )


def complete_pzi_targets(
    family: str,
    family_hits: dict[str, dict[str, dict[str, object]]],
) -> list[str]:
    """Return PZI proteins with complete required HMM domains for a family."""
    requirements = FAMILY_REQUIREMENTS[family]
    if len(requirements) == 1:
        profile = requirements[0]
        return sorted(
            x for x, hit in family_hits[profile].items()
            if x.startswith("pzi_hap") and float(hit["profile_coverage"]) >= 0.5
        )
    common = set(family_hits[requirements[0]])
    for profile in requirements[1:]:
        common &= set(family_hits[profile])
    return sorted(
        x for x in common
        if x.startswith("pzi_hap")
        and all(float(family_hits[p][x]["profile_coverage"]) >= 0.5 for p in requirements)
    )


def make_alignment_inputs(
    out: Path,
    family_hits: dict[str, dict[str, dict[str, object]]],
    combined: dict[str, tuple[str, str]],
    reference_metadata: dict[str, dict[str, str]],
) -> dict[str, list[str]]:
    tree_targets: dict[str, list[str]] = {}
    # Build each tree from complete required-domain hits plus curated
    # Arabidopsis/Rosaceae references.  This keeps broad HMM families out of
    # the formal candidate set until the tree can resolve their subfamily.
    for family, requirements in FAMILY_REQUIREMENTS.items():
        ids = complete_pzi_targets(family, family_hits)
        ids.extend(reference_ids_for_family(family, reference_metadata))
        ids = sorted(set(ids))
        tree_targets[family] = ids
        records = {x: (combined[x][0], combined[x][1]) for x in ids if x in combined}
        write_fasta(out / "sequences" / f"{family}.fa", records)

    ref_rows: list[dict[str, object]] = []
    for target, row in sorted(reference_metadata.items()):
        if row.get("tree_family") in set(TREE_FAMILIES):
            ref_rows.append({
                "tree_family": row.get("tree_family", ""),
                "species": row.get("species", ""),
                "gene_symbol": row.get("gene_symbol", ""),
                "subfamily_group": row.get("subfamily_group", ""),
                "accession": row.get("accession", ""),
                "protein_id": target,
                "selection": row.get("curation_source", ""),
            })
    write_rows(out / "sequences" / "reference_proteins.tsv", ref_rows)
    return tree_targets


def run_align_trees(out: Path, tree_targets: dict[str, list[str]], log: Path, force: bool = False) -> None:
    for family, ids in tree_targets.items():
        inp = out / "sequences" / f"{family}.fa"
        aln = out / "alignments" / f"{family}.mafft.fa"
        trimmed = out / "alignments" / f"{family}.trimmed.fa"
        prefix = out / "trees" / family
        if len(ids) < 4:
            (out / "trees" / f"{family}.status.txt").write_text(f"SKIPPED: only {len(ids)} sequences\n", encoding="utf-8")
            continue
        if not force and (out / "trees" / f"{family}.treefile").exists() and (out / "trees" / f"{family}.treefile").stat().st_size > 0:
            continue
        aln.parent.mkdir(parents=True, exist_ok=True)
        with aln.open("w", encoding="utf-8") as fh:
            subprocess.run([MAFFT, "--auto", "--thread", "4", str(inp)], check=True, stdout=fh, stderr=subprocess.PIPE, text=True)
        try:
            run([TRIMAL, "-in", str(aln), "-out", str(trimmed), "-automated1"], log)
        except subprocess.CalledProcessError:
            shutil.copyfile(aln, trimmed)
        prefix.parent.mkdir(parents=True, exist_ok=True)
        redo = ["-redo"] if force else []
        model = IQTREE_FIXED_MODELS.get(family, "MFP")
        run([IQTREE, "-s", str(trimmed), "-m", model, "-B", "1000", "--alrt", "1000", "--bnni", "-T", "AUTO", "-seed", IQTREE_SEED, *redo, "--prefix", str(prefix)], log)


def parse_newick_tips(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8", errors="ignore")
    # Tree leaves have names immediately before ':' or ',' or ')'.
    return sorted(set(re.findall(r"(?:(?<=\()|(?<=,))([^():;,]+)(?=[:),])", text)))


class TreeNode:
    def __init__(self, name: str = "", length: float = 0.0):
        self.name = name
        self.length = length
        self.children: list[TreeNode] = []
        self.parent: TreeNode | None = None
        self.root_distance = 0.0


def parse_newick(text: str) -> TreeNode:
    """Small Newick parser supporting IQ-TREE support labels such as 98/100."""
    i = 0

    def parse_sub() -> TreeNode:
        nonlocal i
        node = TreeNode()
        if i < len(text) and text[i] == "(":
            i += 1
            while True:
                child = parse_sub()
                child.parent = node
                node.children.append(child)
                if i < len(text) and text[i] == ",":
                    i += 1
                    continue
                if i < len(text) and text[i] == ")":
                    i += 1
                    break
        start = i
        while i < len(text) and text[i] not in ":,();":
            i += 1
        node.name = text[start:i].strip()
        if i < len(text) and text[i] == ":":
            i += 1
            start = i
            while i < len(text) and text[i] not in ",();":
                i += 1
            try:
                node.length = float(text[start:i])
            except ValueError:
                node.length = 0.0
        return node

    return parse_sub()


def support_values(node: TreeNode) -> tuple[float | None, float | None]:
    if not node.name:
        return None, None
    bits = node.name.split("/")
    try:
        first = float(bits[0])
    except ValueError:
        return None, None
    second = None
    if len(bits) > 1:
        try:
            second = float(bits[1])
        except ValueError:
            pass
    return first, second


def tree_placements(
    out: Path,
    tree_targets: dict[str, list[str]],
    reference_metadata: dict[str, dict[str, str]],
) -> None:
    rows: list[dict[str, object]] = []

    def phase_of(target: str) -> str:
        m = re.match(r"pzi_hap(A1|A2|B1|B2)_", target)
        return m.group(1) if m else ""

    for family, ids in tree_targets.items():
        treefile = out / "trees" / f"{family}.treefile"
        if not treefile.exists():
            continue
        root = parse_newick(treefile.read_text(encoding="utf-8", errors="ignore").strip().rstrip(";"))
        leaves: dict[str, TreeNode] = {}

        def walk(node: TreeNode, distance: float) -> None:
            node.root_distance = distance
            if not node.children:
                leaves[node.name] = node
                return
            for child in node.children:
                walk(child, distance + child.length)

        walk(root, 0.0)
        refs = [
            x for x in leaves
            if x in reference_metadata and reference_metadata[x].get("tree_family") == family
        ]
        pzi = [x for x in leaves if x.startswith("pzi_hap")]
        for target in pzi:
            if not refs:
                rows.append({
                    "family": family, "target_id": target, "phase": phase_of(target),
                    "nearest_reference": "", "nearest_reference_symbol": "",
                    "nearest_reference_species": "", "nearest_reference_accession": "",
                    "patristic_distance": "", "mrca_bootstrap": "", "mrca_sh_alrt": "",
                    "reference_group_count": 0, "group_support_count": 0,
                    "placement_status": "no_reference_in_tree", "subfamily_group": "",
                })
                continue
            a = leaves[target]
            ancestors: dict[int, TreeNode] = {}
            n: TreeNode | None = a
            while n is not None:
                ancestors[id(n)] = n
                n = n.parent
            scored: list[tuple[float, str, TreeNode]] = []
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
            boot, alrt = support_values(mrca)
            ref_meta = reference_metadata[ref]
            symbol = ref_meta.get("gene_symbol", "")
            group = ref_meta.get("subfamily_group", "")
            group_refs = [
                r for r in refs
                if reference_metadata[r].get("subfamily_group") == group
            ]
            group_support_count = 0
            for group_ref in group_refs:
                b = leaves[group_ref]
                n = b
                while n is not None and id(n) not in ancestors:
                    n = n.parent
                if n is None:
                    continue
                gb, ga = support_values(n)
                if gb is not None and ga is not None and gb >= 80 and ga >= 80:
                    group_support_count += 1
            supported = boot is not None and alrt is not None and boot >= 70 and alrt >= 70
            rows.append({
                "family": family, "target_id": target, "phase": phase_of(target),
                "nearest_reference": ref, "nearest_reference_symbol": symbol,
                "nearest_reference_species": ref_meta.get("species", ""),
                "nearest_reference_accession": ref_meta.get("accession", ""),
                "patristic_distance": f"{dist:.6g}",
                "mrca_bootstrap": "" if boot is None else f"{boot:.1f}",
                "mrca_sh_alrt": "" if alrt is None else f"{alrt:.1f}",
                "reference_group_count": len(group_refs),
                "group_support_count": group_support_count,
                "placement_status": "supported_reference_placement" if supported else "weak_or_unresolved_placement",
                "subfamily_group": group,
            })
    write_rows(out / "trees" / "tree_placements.tsv", rows)


def build_copy_fate(out: Path, candidates: list[dict[str, object]], phase_by_id: dict[str, str]) -> None:
    wanted_rows = [r for r in candidates if str(r["target_id"]).startswith("pzi_hap")]
    wanted = {str(r["target_id"]): r for r in wanted_rows}
    families_by_gene: dict[str, set[str]] = defaultdict(set)
    for row in wanted_rows:
        families_by_gene[str(row["target_id"])].add(str(row["family"]))
    pairs = [("A1", "A2"), ("A1", "B1"), ("A2", "B2"), ("B1", "B2")]
    hits: list[dict[str, object]] = []
    for left, right in pairs:
        folder = ROOT / "downstream_analysis" / "circos" / "work" / "jcvi" / f"{left}_{right}"
        path = folder / f"hap{left}.hap{right}.anchors.screened"
        if not path.exists():
            path = folder / f"hap{left}.hap{right}.anchors"
        if not path.exists():
            continue
        with path.open(encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 2:
                    continue
                a, b = f[0], f[1]
                for cand, partner in ((a, b), (b, a)):
                    if cand in wanted:
                        hits.append({"target_id": cand, "candidate_phase": phase_by_id.get(cand, ""), "partner_id": partner, "partner_phase_pair": f"{left}_{right}", "anchor_score": f[2] if len(f) > 2 else "", "source": str(path.relative_to(ROOT))})
    write_rows(out / "copy_fate" / "candidate_anchor_hits.tsv", hits)
    by_gene: dict[str, dict[str, object]] = {}
    for gene in wanted:
        fams = families_by_gene[gene]
        row: dict[str, object] = {
            "target_id": gene,
            "phase": phase_by_id.get(gene, ""),
            "PEBP": int("PEBP" in fams),
            "K_box": int("K_box" in fams),
            "SRF_TF": int("SRF_TF" in fams),
            "MADS_both_domains": int("MADS" in fams or {"K_box", "SRF_TF"}.issubset(fams)),
            "C_LFY_FLO": int("C_LFY_FLO" in fams),
            "anchor_hit_count": 0,
            "anchor_pairs": "",
        }
        for family in TREE_FAMILIES:
            if family not in {"PEBP", "MADS", "C_LFY_FLO"}:
                row[family] = int(family in fams)
        by_gene[gene] = row
    for h in hits:
        row = by_gene[h["target_id"]]
        row["anchor_hit_count"] = int(row["anchor_hit_count"]) + 1
        existing = set(filter(None, str(row["anchor_pairs"]).split(";")))
        existing.add(str(h["partner_phase_pair"]))
        row["anchor_pairs"] = ";".join(sorted(existing))
    write_rows(out / "copy_fate" / "candidate_copy_fate.tsv", list(by_gene.values()))
    family_rows: list[dict[str, object]] = []
    family_members: dict[str, set[str]] = {
        family: {g for g, fs in families_by_gene.items() if family in fs}
        for family in TREE_FAMILIES if family != "MADS"
    }
    family_members["MADS_both_domains"] = family_members.get("MADS", set())
    for family, members in family_members.items():
        for phase in PHASES:
            phase_members = {g for g in members if phase_by_id.get(g) == phase}
            expected_pairs = {"A1": {"A1_A2", "A1_B1"}, "A2": {"A1_A2", "A2_B2"}, "B1": {"A1_B1", "B1_B2"}, "B2": {"A2_B2", "B1_B2"}}[phase]
            anchored = {h["target_id"] for h in hits if h["target_id"] in phase_members and h["partner_phase_pair"] in expected_pairs}
            family_rows.append({"family": family, "phase": phase, "candidate_count": len(phase_members), "anchored_candidate_count": len(anchored), "unanchored_candidate_count": len(phase_members - anchored), "anchor_rate": f"{len(anchored) / len(phase_members):.4f}" if phase_members else "", "interpretation": "local_synteny_supported" if anchored else "no_local_anchor_or_unresolved"})
    write_rows(out / "copy_fate" / "copy_fate_family_summary.tsv", family_rows)


def merge_final_evidence(out: Path) -> None:
    fate_path = out / "copy_fate" / "candidate_copy_fate.tsv"
    placement_path = out / "trees" / "tree_placements.tsv"
    if not fate_path.exists():
        return
    fate = list(csv.DictReader(fate_path.open(encoding="utf-8"), delimiter="\t"))
    placements = defaultdict(list)
    if placement_path.exists():
        for row in csv.DictReader(placement_path.open(encoding="utf-8"), delimiter="\t"):
            placements[row["target_id"]].append(row)
    rows: list[dict[str, object]] = []
    for row in fate:
        ps = placements.get(row["target_id"], [])
        p = next((item for item in ps if item.get("family") == row.get("family")), ps[0] if ps else {})
        mads = row.get("MADS_both_domains", "0") == "1"
        boot = float(p["mrca_bootstrap"]) if p.get("mrca_bootstrap") else None
        alrt = float(p["mrca_sh_alrt"]) if p.get("mrca_sh_alrt") else None
        group_count = int(p.get("reference_group_count", "0") or 0)
        group_support = int(p.get("group_support_count", "0") or 0)
        anchor_count = int(row.get("anchor_hit_count", "0") or 0)
        strict_support = boot is not None and alrt is not None and boot >= 80 and alrt >= 80
        confirmed_by_evidence = bool(
            p
            and p.get("placement_status") == "supported_reference_placement"
            and strict_support
            and group_count >= 2
            and group_support >= 2
            and anchor_count > 0
            and p.get("subfamily_group")
        )
        user_excluded = row["target_id"] in USER_EXCLUDED_TARGETS
        confirmed = confirmed_by_evidence and not user_excluded
        reasons: list[str] = []
        if user_excluded:
            reasons.append("user_excluded_lfy_flo")
        if not p:
            reasons.append("no_tree_placement")
        if p and not strict_support:
            reasons.append("support_below_80_80")
        if p and group_count < 2:
            reasons.append("fewer_than_two_same_group_references")
        if p and group_support < 2:
            reasons.append("fewer_than_two_supported_same_group_references")
        if anchor_count == 0:
            reasons.append("no_screened_local_anchor")
        family = row.get("family", "")
        output_row: dict[str, object] = {
            "target_id": row["target_id"],
            "phase": row["phase"],
            "domain_class": "MADS_both_domains" if mads else family,
            "PEBP_domain": row.get("PEBP", "0"),
            "K_box_domain": row.get("K_box", "0"),
            "SRF_TF_domain": row.get("SRF_TF", "0"),
            "C_LFY_FLO_domain": row.get("C_LFY_FLO", "0"),
            "tree_family": p.get("family", family),
            "nearest_reference_symbol": p.get("nearest_reference_symbol", ""),
            "nearest_reference_species": p.get("nearest_reference_species", ""),
            "nearest_reference_accession": p.get("nearest_reference_accession", ""),
            "subfamily_group": p.get("subfamily_group", ""),
            "mrca_bootstrap": p.get("mrca_bootstrap", ""),
            "mrca_sh_alrt": p.get("mrca_sh_alrt", ""),
            "reference_group_count": p.get("reference_group_count", "0"),
            "group_support_count": p.get("group_support_count", "0"),
            "placement_status": p.get("placement_status", "no_full_length_tree"),
            "anchor_hit_count": row.get("anchor_hit_count", "0"),
            "anchor_pairs": row.get("anchor_pairs", ""),
            "copy_fate_interpretation": "local_synteny_supported" if int(row.get("anchor_hit_count", "0") or 0) > 0 else "no_local_anchor_or_unresolved",
            "final_status": "confirmed" if confirmed else "excluded",
            "exclusion_reason": ";".join(reasons),
        }
        for tree_family in TREE_FAMILIES:
            if tree_family not in {"PEBP", "MADS", "C_LFY_FLO"}:
                output_row[f"{tree_family}_domain"] = row.get(tree_family, "0")
        rows.append(output_row)
    write_rows(out / "final_candidate_evidence.tsv", rows)


def build_homology_components(
    out: Path,
    fate_filename: str = "candidate_copy_fate.tsv",
    output_filename: str = "candidate_synteny_components.tsv",
) -> None:
    """Build anchor-connected candidate components; these are synteny blocks, not orthogroups."""
    fate_path = out / "copy_fate" / fate_filename
    hit_path = out / "copy_fate" / "candidate_anchor_hits.tsv"
    if not fate_path.exists() or not hit_path.exists():
        return
    fate = list(csv.DictReader(fate_path.open(encoding="utf-8"), delimiter="\t"))
    targets = {r["target_id"] for r in fate}
    parent = {x: x for x in targets}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    hits = list(csv.DictReader(hit_path.open(encoding="utf-8"), delimiter="\t"))
    noncandidate = defaultdict(int)
    for h in hits:
        a, b = h["target_id"], h["partner_id"]
        if a in targets and b in targets:
            union(a, b)
        elif a in targets:
            noncandidate[a] += 1
    groups: dict[str, list[str]] = defaultdict(list)
    for target in sorted(targets):
        groups[find(target)].append(target)
    fate_by_gene = {r["target_id"]: r for r in fate}
    rows: list[dict[str, object]] = []
    for idx, members in enumerate(sorted(groups.values(), key=lambda x: x[0]), 1):
        phases = sorted({fate_by_gene[x]["phase"] for x in members})
        phase_pattern = "+".join(p for p in PHASES if p in phases)
        n_non = sum(noncandidate[x] for x in members)
        if set(phases) == set(PHASES):
            status = "four_phase_anchor_component"
        elif len(phases) > 1:
            status = "partial_phase_anchor_component"
        elif n_non:
            status = "candidate_anchor_to_non_candidate_partner"
        else:
            status = "singleton_or_unresolved"
        classes = set()
        for x in members:
            r = fate_by_gene[x]
            if r.get("PEBP") == "1": classes.add("PEBP")
            if r.get("MADS_both_domains") == "1": classes.add("MADS_both_domains")
            if r.get("C_LFY_FLO") == "1": classes.add("C_LFY_FLO")
            for family in TREE_FAMILIES:
                if family not in {"PEBP", "MADS", "C_LFY_FLO"} and r.get(family) == "1":
                    classes.add(family)
            if r.get("K_box") == "1" and r.get("MADS_both_domains") != "1": classes.add("MADS_single_domain")
            if r.get("SRF_TF") == "1" and r.get("MADS_both_domains") != "1": classes.add("MADS_single_domain")
        rows.append({"component_id": f"SYNCOMP_{idx:04d}", "members": ";".join(members), "n_members": len(members), "phase_pattern": phase_pattern, "family_classes": ";".join(sorted(classes)), "noncandidate_anchor_count": n_non, "copy_fate_status": status})
    write_rows(out / "copy_fate" / output_filename, rows)


def build_confirmed_ab_copy_fate(out: Path) -> None:
    """Summarize confirmed genes as A1/A2/B1/B2 synteny components.

    The table is intentionally conservative: a phase is reported as
    ``*_confirmed`` only when a confirmed candidate occurs in that phase.
    A missing phase is ``*_not_confirmed`` and is never interpreted as a
    biological gene loss without independent evidence.
    """
    fate_path = out / "copy_fate" / "confirmed_copy_fate.tsv"
    component_path = out / "copy_fate" / "confirmed_synteny_components.tsv"
    anchor_path = out / "copy_fate" / "confirmed_anchor_hits.tsv"
    if not fate_path.exists() or not component_path.exists():
        return

    fate_rows = list(csv.DictReader(fate_path.open(encoding="utf-8"), delimiter="\t"))
    components = list(csv.DictReader(component_path.open(encoding="utf-8"), delimiter="\t"))
    confirmed_ids = {r["target_id"] for r in fate_rows}
    fate_by_id = {r["target_id"]: r for r in fate_rows}
    # Tree placement and subfamily annotations live in confirmed_candidates,
    # whereas confirmed_copy_fate is deliberately a compact phase/anchor
    # table.  Join the two tables by target_id for the final matrix.
    metadata_by_id: dict[str, dict[str, str]] = {}
    metadata_path = out / "confirmed_candidates.tsv"
    if metadata_path.exists():
        metadata_by_id = {
            r["target_id"]: r
            for r in csv.DictReader(metadata_path.open(encoding="utf-8"), delimiter="\t")
            if r.get("target_id")
        }
    component_by_id: dict[str, str] = {}
    members_by_component: dict[str, list[str]] = {}
    for component in components:
        component_id = component["component_id"]
        members = [x for x in component.get("members", "").split(";") if x in confirmed_ids]
        members_by_component[component_id] = members
        for member in members:
            component_by_id[member] = component_id

    # Anchor files contain both directions for a confirmed-confirmed pair.
    # Deduplicate those edges and keep non-confirmed partners as unresolved
    # evidence rather than treating them as losses.
    internal_edges: dict[str, set[tuple[str, str]]] = defaultdict(set)
    unresolved_edges: dict[str, set[tuple[str, str]]] = defaultdict(set)
    phase_pair_values: dict[str, set[str]] = defaultdict(set)
    if anchor_path.exists():
        for hit in csv.DictReader(anchor_path.open(encoding="utf-8"), delimiter="\t"):
            target = hit.get("target_id", "")
            partner = hit.get("partner_id", "")
            component_id = component_by_id.get(target)
            if not component_id or not partner:
                continue
            phase_pair = hit.get("partner_phase_pair", "")
            if phase_pair:
                phase_pair_values[component_id].add(phase_pair)
            edge = tuple(sorted((target, partner)))
            if partner in confirmed_ids and component_by_id.get(partner) == component_id:
                internal_edges[component_id].add(edge)
            elif partner not in confirmed_ids:
                unresolved_edges[component_id].add(edge)

    output_rows: list[dict[str, object]] = []
    for component in components:
        component_id = component["component_id"]
        members = members_by_component.get(component_id, [])
        if not members:
            continue
        phase_ids = {
            phase: sorted(x for x in members if fate_by_id.get(x, {}).get("phase") == phase)
            for phase in PHASES
        }
        present = {phase: int(bool(ids)) for phase, ids in phase_ids.items()}
        has_a = bool(present["A1"] or present["A2"])
        has_b = bool(present["B1"] or present["B2"])
        complete_a = bool(present["A1"] and present["A2"])
        complete_b = bool(present["B1"] and present["B2"])
        if complete_a and complete_b:
            copy_fate = "A_and_B_confirmed"
        elif complete_a and not has_b:
            copy_fate = "A_only_confirmed"
        elif complete_b and not has_a:
            copy_fate = "B_only_confirmed"
        else:
            copy_fate = "no_complete_confirmed_copy_set"
        family_values = sorted({metadata_by_id.get(x, {}).get("tree_family", "") for x in members if metadata_by_id.get(x, {}).get("tree_family")})
        subfamily_values = sorted({metadata_by_id.get(x, {}).get("subfamily_group", "") for x in members if metadata_by_id.get(x, {}).get("subfamily_group")})
        output_rows.append({
            "copy_group_id": component_id,
            "tree_family": ";".join(x for x in family_values if x),
            "subfamily_group": ";".join(subfamily_values),
            "n_confirmed_members": len(members),
            "A1_ids": ";".join(phase_ids["A1"]),
            "A2_ids": ";".join(phase_ids["A2"]),
            "B1_ids": ";".join(phase_ids["B1"]),
            "B2_ids": ";".join(phase_ids["B2"]),
            "A1_present": present["A1"],
            "A2_present": present["A2"],
            "B1_present": present["B1"],
            "B2_present": present["B2"],
            "A_subgenome_status": "A_confirmed" if has_a else "A_not_confirmed",
            "B_subgenome_status": "B_confirmed" if has_b else "B_not_confirmed",
            "subgenome_copy_fate": copy_fate,
            "A_copy_count": len(phase_ids["A1"]) + len(phase_ids["A2"]),
            "B_copy_count": len(phase_ids["B1"]) + len(phase_ids["B2"]),
            "phase_pattern": "+".join(phase for phase in PHASES if present[phase]),
            "anchor_pairs": ";".join(sorted(phase_pair_values.get(component_id, set()))),
            "confirmed_anchor_edge_count": len(internal_edges.get(component_id, set())),
            "nonconfirmed_partner_edge_count": len(unresolved_edges.get(component_id, set())),
            "interpretation": "confirmed syntenic copies; missing phase is unresolved, not a gene-loss call",
        })
    write_rows(out / "copy_fate" / "confirmed_ab_copy_fate.tsv", output_rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--force-trees", action="store_true", help="recompute existing IQ-TREE outputs")
    args = ap.parse_args()
    out = args.outdir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    (out / "hmmsearch").mkdir(parents=True, exist_ok=True)
    (out / "sequences").mkdir(parents=True, exist_ok=True)
    (out / "alignments").mkdir(parents=True, exist_ok=True)
    (out / "trees").mkdir(parents=True, exist_ok=True)
    (out / "copy_fate").mkdir(parents=True, exist_ok=True)
    log = out / "pipeline.log"
    if not HMM_DB.exists():
        raise FileNotFoundError(HMM_DB)
    at_pep = REF_DIR / "Arabidopsis_thaliana.TAIR10.pep.all.fa.gz"
    if not at_pep.exists():
        raise FileNotFoundError(at_pep)
    rosaceae_meta, rosaceae_records = load_rosaceae_references()
    if not rosaceae_meta or not rosaceae_records:
        raise FileNotFoundError(
            "Curated Rosaceae references are missing; run scripts/prepare_rosaceae_references.py first"
        )

    profiles = extract_profiles(HMM_DB, out / "reference" / "selected_profiles.hmm")
    pzi_fa = out / "targets" / "pzi_A1A2B1B2.primary.pep.fa"
    phase_by_id = make_pzi_fasta(pzi_fa)
    at_symbols, combined = make_combined_target(
        pzi_fa,
        at_pep,
        out / "targets" / "combined.pep.fa",
        external_records=rosaceae_records,
    )
    at_records = {k: v for k, v in combined.items() if k in at_symbols}
    reference_metadata = build_reference_metadata(at_symbols, at_records, rosaceae_meta)

    all_rows: list[dict[str, object]] = []
    family_hits: dict[str, dict[str, dict[str, object]]] = {}
    for family, acc in PROFILES.items():
        domtbl = out / "hmmsearch" / f"{family}.domtblout"
        run([HMMSEARCH, "--cut_ga", "--cpu", str(args.threads), "--domtblout", str(domtbl), "--noali", str(out / "reference" / f"{family}.hmm"), str(out / "targets" / "combined.pep.fa")], log)
        rows = parse_domtbl(domtbl, family, at_symbols)
        best = best_hits(rows)
        family_hits[family] = best
        all_rows.extend(rows)

    candidate_rows: list[dict[str, object]] = []
    for family, hits in family_hits.items():
        for target, row in hits.items():
            is_pzi = target.startswith("pzi_hap")
            symbol = str(row.get("arabidopsis_symbol", ""))
            is_reference = (
                not is_pzi
                and target in reference_metadata
                and reference_metadata[target].get("tree_family") in PROFILE_TREE_FAMILIES.get(family, set())
            )
            if not is_pzi and not is_reference:
                continue
            row = dict(row)
            row["domain_quality"] = "full_domain" if float(row["profile_coverage"]) >= 0.5 else "partial_domain_review"
            ref_meta = reference_metadata.get(target, {})
            candidate_rows.append({
                "family": family,
                "target_id": target,
                "phase": phase_by_id.get(target, "Arabidopsis"),
                "arabidopsis_symbol": symbol,
                "reference_species": ref_meta.get("species", ""),
                "reference_gene_symbol": ref_meta.get("gene_symbol", ""),
                "reference_subfamily_group": ref_meta.get("subfamily_group", ""),
                "is_reference": int(is_reference),
                **row,
            })
    write_rows(out / "hmmsearch" / "all_domain_hits.tsv", all_rows)
    write_rows(out / "hmmsearch" / "candidate_domain_hits.tsv", candidate_rows)

    # Write a generic domain-integration table for all composite families and
    # retain the historical MADS-specific filename for compatibility.
    composite_rows: list[dict[str, object]] = []
    for family, requirements in FAMILY_REQUIREMENTS.items():
        if len(requirements) < 2:
            continue
        targets = set().union(*(set(family_hits[p]) for p in requirements))
        for target in sorted(targets):
            if not target.startswith("pzi_hap") and target not in reference_metadata:
                continue
            row: dict[str, object] = {
                "family": family,
                "target_id": target,
                "phase": phase_by_id.get(target, "Arabidopsis"),
                "reference_species": reference_metadata.get(target, {}).get("species", ""),
                "reference_gene_symbol": reference_metadata.get(target, {}).get("gene_symbol", ""),
                "reference_subfamily_group": reference_metadata.get(target, {}).get("subfamily_group", ""),
                "complete_required_domains": int(all(p in family_hits and target in family_hits[p] for p in requirements)),
                "complete_required_domains_strict": int(target in complete_pzi_targets(family, family_hits) or (target in reference_metadata and all(p in family_hits and target in family_hits[p] for p in requirements))),
            }
            for profile in requirements:
                hit = family_hits[profile].get(target)
                row[f"{profile}_hit"] = int(hit is not None)
                row[f"{profile}_coverage"] = hit["profile_coverage"] if hit else ""
                row[f"{profile}_iEvalue"] = hit["domain_i_evalue"] if hit else ""
            composite_rows.append(row)
    write_rows(out / "hmmsearch" / "domain_integration.tsv", composite_rows)
    mads_rows = [r for r in composite_rows if r.get("family") == "MADS"]
    write_rows(out / "hmmsearch" / "mads_domain_integration.tsv", mads_rows)

    # The strict analysis universe is one row per PZI protein/family with all
    # required domains complete.  Family identity is resolved only below by
    # the reference-anchored maximum-likelihood trees.
    strict_candidates: list[dict[str, object]] = []
    for family in TREE_FAMILIES:
        for target in complete_pzi_targets(family, family_hits):
            strict_candidates.append({
                "family": family,
                "target_id": target,
                "phase": phase_by_id.get(target, ""),
                "domain_quality": "full_domain",
                "is_reference": 0,
            })
    strict_candidates.sort(key=lambda r: (str(r["family"]), str(r["target_id"])))
    write_rows(out / "strict_candidates.tsv", strict_candidates)

    tree_targets = make_alignment_inputs(out, family_hits, combined, reference_metadata)
    run_align_trees(out, tree_targets, log, force=args.force_trees)
    tree_rows: list[dict[str, object]] = []
    for family in tree_targets:
        treefile = out / "trees" / f"{family}.treefile"
        tree_rows.append({"family": family, "n_sequences": len(tree_targets[family]), "treefile": str(treefile.relative_to(ROOT)) if treefile.exists() else "", "tree_status": "completed" if treefile.exists() else "missing"})
    write_rows(out / "trees" / "tree_manifest.tsv", tree_rows)
    tree_placements(out, tree_targets, reference_metadata)

    # Copy-fate uses only the strict 134-protein universe; absence of an anchor
    # is unresolved until the final confirmation filter is applied.
    build_copy_fate(out, strict_candidates, phase_by_id)
    merge_final_evidence(out)
    final_rows = list(csv.DictReader((out / "final_candidate_evidence.tsv").open(encoding="utf-8"), delimiter="\t"))
    confirmed = [r for r in final_rows if r.get("final_status") == "confirmed"]
    write_rows(out / "confirmed_candidates.tsv", confirmed)
    write_rows(out / "high_confidence_candidates.tsv", confirmed)

    confirmed_ids = {r["target_id"] for r in confirmed}
    fate_path = out / "copy_fate" / "candidate_copy_fate.tsv"
    if fate_path.exists():
        fate_rows = list(csv.DictReader(fate_path.open(encoding="utf-8"), delimiter="\t"))
        confirmed_fate = [r for r in fate_rows if r.get("target_id") in confirmed_ids]
        write_rows(out / "copy_fate" / "confirmed_copy_fate.tsv", confirmed_fate)
    anchor_path = out / "copy_fate" / "candidate_anchor_hits.tsv"
    if anchor_path.exists():
        anchor_rows = list(csv.DictReader(anchor_path.open(encoding="utf-8"), delimiter="\t"))
        confirmed_anchor_rows = [r for r in anchor_rows if r.get("target_id") in confirmed_ids]
        write_rows(out / "copy_fate" / "confirmed_anchor_hits.tsv", confirmed_anchor_rows)
    build_homology_components(
        out,
        fate_filename="confirmed_copy_fate.tsv",
        output_filename="confirmed_synteny_components.tsv",
    )
    build_confirmed_ab_copy_fate(out)

    # Phase-level copy-fate summary restricted to confirmed candidates.
    confirmed_fate_path = out / "copy_fate" / "confirmed_copy_fate.tsv"
    confirmed_anchor_path = out / "copy_fate" / "confirmed_anchor_hits.tsv"
    if confirmed_fate_path.exists():
        confirmed_fate_rows = list(csv.DictReader(confirmed_fate_path.open(encoding="utf-8"), delimiter="\t"))
        confirmed_anchor_rows = list(csv.DictReader(confirmed_anchor_path.open(encoding="utf-8"), delimiter="\t")) if confirmed_anchor_path.exists() else []
        confirmed_summary: list[dict[str, object]] = []
        summary_predicates = [
            ("PEBP", lambda r: r.get("PEBP") == "1"),
            ("MADS_both_domains", lambda r: r.get("MADS_both_domains") == "1"),
            ("C_LFY_FLO", lambda r: r.get("C_LFY_FLO") == "1"),
        ]
        summary_predicates.extend(
            (family, lambda r, family=family: r.get(family) == "1")
            for family in TREE_FAMILIES
            if family not in {"PEBP", "MADS", "C_LFY_FLO"}
        )
        for family, predicate in summary_predicates:
            for phase in PHASES:
                members = {r["target_id"] for r in confirmed_fate_rows if r.get("phase") == phase and predicate(r)}
                anchored = {r["target_id"] for r in confirmed_anchor_rows if r.get("target_id") in members}
                confirmed_summary.append({
                    "family": family,
                    "phase": phase,
                    "confirmed_count": len(members),
                    "confirmed_anchored_count": len(anchored),
                    "confirmed_unanchored_count": len(members - anchored),
                    "anchor_rate": f"{len(anchored) / len(members):.4f}" if members else "",
                })
        write_rows(out / "copy_fate" / "confirmed_copy_fate_family_summary.tsv", confirmed_summary)

    family_summary: list[dict[str, object]] = []
    for family in TREE_FAMILIES:
        members = [r for r in confirmed if r.get("tree_family") == family]
        groups: dict[str, int] = defaultdict(int)
        for row in members:
            groups[row.get("subfamily_group", "unresolved")] += 1
        for group, count in sorted(groups.items()):
            family_summary.append({"family": family, "subfamily_group": group, "confirmed_count": count})
    write_rows(out / "confirmed_family_summary.tsv", family_summary)

    summary = [
        {"metric": "PZI_proteins_scanned", "value": len(phase_by_id)},
        {"metric": "Arabidopsis_proteins_scanned", "value": len(at_records)},
        {"metric": "Rosaceae_reference_proteins", "value": len(rosaceae_records)},
        {"metric": "strict_candidate_records", "value": len(strict_candidates)},
        {"metric": "strict_candidate_unique_proteins", "value": len({r["target_id"] for r in strict_candidates})},
        {"metric": "confirmed_candidate_records", "value": len(confirmed)},
        {"metric": "confirmed_candidate_unique_proteins", "value": len({r["target_id"] for r in confirmed})},
        {"metric": "PZI_MADS_both_domain_proteins", "value": sum(1 for r in mads_rows if str(r["target_id"]).startswith("pzi_hap") and r["complete_required_domains"] == 1)},
        {"metric": "PZI_MADS_strict_both_domain_proteins", "value": sum(1 for r in mads_rows if str(r["target_id"]).startswith("pzi_hap") and r["complete_required_domains_strict"] == 1)},
        {"metric": "reference_anchor_proteins_in_trees", "value": sum(1 for family in TREE_FAMILIES for r in tree_targets[family] if r in reference_metadata)},
    ]
    write_rows(out / "pipeline_summary.tsv", summary)
    (out / "README.md").write_text(
        "# Inflorescence family evidence pipeline\n\n"
        "PZI proteins were screened with Pfam-A profiles for PEBP (PF01161), MADS K-box/SRF-TF (PF01486/PF00319), LFY/FLO (PF17538), AP2 (PF00847), SPL/SBP (PF03110), TCP (PF03634), KNOX1/KNOX2 (PF03790/PF03791), Homeobox (PF00046) and YABBY (PF04690) using HMMER hmmsearch --cut_ga. Composite MADS and KNOX families require both complete domains. All families are resolved with reference-anchored MAFFT/trimAl/IQ-TREE analyses (1,000 UFBoot and 1,000 SH-aLRT replicates).\n\n"
        f"The strict analysis universe contains {len(strict_candidates)} PZI protein-family records across {len(TREE_FAMILIES)} tree families. A candidate is labelled confirmed only when it has complete domain evidence, a reference-anchored placement with UFBoot and SH-aLRT at least 80, at least two independent references in the same annotated subfamily group, at least two of those references supporting the target-side clade, and at least one screened local JCVI anchor. All other strict candidates remain in final_candidate_evidence.tsv as excluded with explicit reasons.\n\n"
        "Run `python scripts/test_ab_bias.py --matrix <outdir>/copy_fate/confirmed_ab_copy_fate.tsv` for exact component-level A/B tests. Missing or unconfirmed phases remain unresolved and are not called gene losses.\n\n"
        "Key outputs include hmmsearch/all_domain_hits.tsv, hmmsearch/domain_integration.tsv, strict_candidates.tsv, trees/tree_placements.tsv, final_candidate_evidence.tsv, confirmed_candidates.tsv, copy_fate/confirmed_ab_copy_fate.tsv and formal A/B-bias reports.\n",
        encoding="utf-8",
    )
    print(out)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
