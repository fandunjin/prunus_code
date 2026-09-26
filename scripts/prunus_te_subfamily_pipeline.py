#!/usr/bin/env python3
"""Extract copy-level TE sequences and build supported subfamily trees.

This script is intended for the HPC project because it requires the frozen
pan-TE masking outputs, the 12 target genome FASTAs and external tools. It
never assigns a named subfamily from a consensus sequence alone.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


FAMILY_RE = re.compile(r"^(panTE_\d+(?:_(?:INT|LTR))?)$")
TARGET_ORDER = [
    "pziA1", "pziA2", "pziB1", "pziB2",
    "pspA1", "pspA2", "pspB1", "pspB2",
    "ppaA", "ppaB", "pcam", "ppe",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run(command: list[str], log: Path, check: bool = True) -> subprocess.CompletedProcess:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("a", encoding="utf-8") as handle:
        handle.write("$ " + " ".join(map(str, command)) + "\n")
        proc = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, text=True)
    if check and proc.returncode:
        raise RuntimeError(f"command failed ({proc.returncode}): {' '.join(command)}; see {log}")
    return proc


def check_tool(path_or_name: str, label: str) -> str:
    resolved = shutil.which(path_or_name) if not Path(path_or_name).is_file() else path_or_name
    if not resolved:
        raise SystemExit(f"ERROR: required tool unavailable: {label} ({path_or_name})")
    return resolved


def parse_fasta(path: Path):
    name = None
    chunks: list[str] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def read_library_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    for header, sequence in parse_fasta(path):
        match = re.search(r"(panTE_\d+(?:_(?:INT|LTR))?)", header)
        if match:
            family = match.group(1)
            lengths[family] = max(lengths.get(family, 0), len(sequence))
    return lengths


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_library(library: Path, manifest: Path) -> str:
    rows = read_tsv(manifest)
    expected = ""
    for row in rows:
        manifest_library = Path(row.get("library", "")).name
        if manifest_library == library.name:
            expected = row.get("sha256_file", "")
            break
    checksum_file = library.with_name(library.name + ".sha256")
    checksum_candidates = [checksum_file]
    if expected and Path(expected).is_file():
        checksum_candidates.insert(0, Path(expected))
    checksum_text = ""
    for candidate in checksum_candidates:
        if candidate.is_file():
            tokens = candidate.read_text(encoding="utf-8", errors="replace").split()
            if tokens:
                checksum_text = tokens[0]
                break
    if not checksum_text and re.fullmatch(r"[0-9a-fA-F]{64}", expected or ""):
        checksum_text = expected
    if not checksum_text:
        raise SystemExit(
            "ERROR: no local SHA256 checksum was available for the frozen pan-TE library. "
            f"Place {library.name}.sha256 beside the FASTA."
        )
    actual = sha256(library)
    if actual.lower() != checksum_text.lower():
        raise SystemExit(f"ERROR: frozen pan-TE SHA256 mismatch: {library}\nexpected={checksum_text}\nactual={actual}")
    return actual


def select_families(selection: Path) -> list[dict[str, str]]:
    rows = read_tsv(selection)
    if not rows:
        raise SystemExit(f"ERROR: empty selection table: {selection}")
    required = {"pan_family", "tree_mode"}
    missing = required - set(rows[0])
    if missing:
        raise SystemExit(f"ERROR: selection table missing columns: {sorted(missing)}")
    return rows


def parse_rmout(path: Path, selected: set[str]) -> list[dict[str, object]]:
    hits: list[dict[str, object]] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 11 or not fields[0].isdigit():
                continue
            family = fields[9]
            if family not in selected:
                continue
            try:
                start = int(fields[5]) - 1
                end = int(fields[6])
            except ValueError:
                continue
            if end <= start:
                continue
            strand = "-" if fields[8].upper().startswith("C") else "+"
            hits.append({
                "seqid": fields[4], "start": start, "end": end,
                "strand": strand, "family": family,
                "class_family": fields[10], "score": fields[0],
                "divergence": fields[1],
            })
    return hits


def merge_hits(hits: list[dict[str, object]], max_gap: int = 100) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for hit in hits:
        grouped[(str(hit["family"]), str(hit["seqid"]), str(hit["strand"]))].append(hit)
    merged: list[dict[str, object]] = []
    for (family, seqid, strand), group in grouped.items():
        group.sort(key=lambda x: (int(x["start"]), int(x["end"])))
        current = None
        for hit in group:
            if current is None:
                current = dict(hit)
                current["hit_count"] = 1
                continue
            if int(hit["start"]) <= int(current["end"]) + max_gap:
                current["end"] = max(int(current["end"]), int(hit["end"]))
                current["hit_count"] = int(current["hit_count"]) + 1
            else:
                merged.append(current)
                current = dict(hit)
                current["hit_count"] = 1
        if current is not None:
            merged.append(current)
    for index, record in enumerate(merged, 1):
        record["length_bp"] = int(record["end"]) - int(record["start"])
        record["copy_local_id"] = index
    return merged


def write_beds_and_catalogue(args, selected_rows, library_lengths):
    selected = {row["pan_family"] for row in selected_rows}
    by_target: dict[str, list[dict[str, object]]] = defaultdict(list)
    for target in TARGET_ORDER:
        rmout = args.mask_root / target / f"{target}.fa.out"
        if not rmout.is_file():
            raise SystemExit(f"ERROR: missing standardized RepeatMasker output: {rmout}")
        raw = parse_rmout(rmout, selected)
        merged = merge_hits(raw, max_gap=args.merge_gap)
        grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
        for record in merged:
            grouped[str(record["family"])].append(record)
        for family, records in grouped.items():
            ref_len = library_lengths.get(family, 0)
            min_len = max(args.min_copy_bp, int(ref_len * args.min_copy_fraction))
            max_len = int(ref_len * args.max_copy_fraction) if ref_len else 0
            records = [r for r in records if int(r["length_bp"]) >= min_len and (not max_len or int(r["length_bp"]) <= max_len)]
            records.sort(key=lambda r: (int(r["length_bp"]), int(r["hit_count"])), reverse=True)
            for record in records[:args.max_copies_per_target]:
                record["target"] = target
                by_target[target].append(record)

    bed_dir = args.outdir / "beds"
    bed_dir.mkdir(parents=True, exist_ok=True)
    catalogue: list[dict[str, object]] = []
    for target, records in by_target.items():
        bed_path = bed_dir / f"{target}.bed"
        with bed_path.open("w", encoding="utf-8") as handle:
            for ordinal, record in enumerate(records, 1):
                name = f"{record['family']}|{target}|{record['seqid']}|{record['start'] + 1}|{record['end']}|copy{ordinal}"
                record["bed_name"] = name
                handle.write(f"{record['seqid']}\t{record['start']}\t{record['end']}\t{name}\t0\t{record['strand']}\n")
                catalogue.append(record)
    fields = ["target", "family", "seqid", "start", "end", "strand", "length_bp", "hit_count", "class_family", "bed_name"]
    write_tsv(args.outdir / "copy_catalogue.tsv", catalogue, fields)
    return bed_dir, catalogue


def extract_sequences(args, bed_dir, catalogue):
    family_dir = args.outdir / "family_fastas"
    family_dir.mkdir(parents=True, exist_ok=True)
    for stale in family_dir.glob("*.fa"):
        stale.unlink()
    for target in TARGET_ORDER:
        bed = bed_dir / f"{target}.bed"
        if not bed.is_file() or bed.stat().st_size == 0:
            continue
        genome = args.mask_root / target / f"{target}.fa"
        fasta = args.outdir / "target_fastas" / f"{target}.fa"
        fasta.parent.mkdir(parents=True, exist_ok=True)
        run([args.bedtools, "getfasta", "-fi", str(genome), "-bed", str(bed), "-s", "-name", "-fo", str(fasta)], args.outdir / "logs" / f"{target}.getfasta.log")
        sequences = list(parse_fasta(fasta))
        by_family: dict[str, list[tuple[str, str]]] = defaultdict(list)
        for name, sequence in sequences:
            family = name.split("|", 1)[0]
            by_family[family].append((name, sequence))
        for family, records in by_family.items():
            path = family_dir / f"{family}.fa"
            with path.open("a", encoding="utf-8") as handle:
                for name, sequence in records:
                    handle.write(f">{name}\n{sequence}\n")
    return family_dir


def cluster_family(args, fasta: Path, family: str) -> str:
    count = sum(1 for _ in parse_fasta(fasta))
    if count < 3:
        return "too_few_copy_sequences"
    outdir = args.outdir / "clusters"
    outdir.mkdir(parents=True, exist_ok=True)
    if args.mmseqs:
        prefix = outdir / family
        tmp = outdir / f"{family}.tmp"
        run([args.mmseqs, "easy-cluster", str(fasta), str(prefix), str(tmp), "--min-seq-id", str(args.cluster_identity), "-c", str(args.cluster_coverage), "--cov-mode", "1", "--threads", str(args.threads)], args.outdir / "logs" / f"{family}.mmseqs.log")
        return "mmseqs2"
    output = outdir / f"{family}.cdhit.fa"
    run([args.cdhit, "-i", str(fasta), "-o", str(output), "-c", str(args.cluster_identity), "-n", "5", "-d", "0", "-T", str(args.threads), "-M", "0"], args.outdir / "logs" / f"{family}.cdhit.log")
    return "cd-hit-est"


def read_domain_records(path: Path, genes: set[str], expected_clade: str = "") -> list[tuple[str, str, str]]:
    records = []
    for name, sequence in parse_fasta(path):
        header = name
        gene = ""
        clade = ""
        for token in header.split(";"):
            if token.startswith("gene="):
                gene = token.split("=", 1)[1]
            if token.startswith("clade="):
                clade = token.split("=", 1)[1]
        if gene in genes and (not expected_clade or clade == expected_clade):
            records.append((name, sequence, gene))
    return records


def run_domain_trees(args, selected_rows, family_dir):
    status_rows = []
    domain_dir = args.outdir / "domains"
    alignment_dir = args.outdir / "alignments"
    tree_dir = args.outdir / "trees"
    for subdir in (domain_dir, alignment_dir, tree_dir):
        subdir.mkdir(parents=True, exist_ok=True)
    for row in selected_rows:
        family = row["pan_family"]
        mode = row["tree_mode"]
        fasta = family_dir / f"{family}.fa"
        result = dict(row)
        result["copy_fasta"] = str(fasta)
        if not fasta.is_file() or sum(1 for _ in parse_fasta(fasta)) < 3:
            result["status"] = "too_few_copy_sequences"
            status_rows.append(result)
            continue
        if mode == "Gypsy_RT_RH_INT":
            prefix = args.outdir / "tesorter" / family
            prefix.parent.mkdir(parents=True, exist_ok=True)
            run([args.tesorter, str(fasta), "-db", "rexdb-plant", "-p", str(args.threads)], args.outdir / "logs" / f"{family}.tesorter.log")
            candidates = list(fasta.parent.glob(fasta.name + "*.rexdb-plant.dom.faa"))
            candidates += list(prefix.parent.glob(prefix.name + "*.rexdb-plant.dom.faa"))
            if not candidates:
                candidates = list(args.outdir.rglob(f"{family}*.rexdb-plant.dom.faa"))
            if not candidates:
                result["status"] = "missing_rexdb_domain_output"
                status_rows.append(result)
                continue
            records = read_domain_records(candidates[-1], {"RT", "RH", "INT"}, str(row.get("gypsy_lineage", "")).replace("Gypsy/", ""))
            if len(records) < 3:
                result["status"] = "insufficient_RT_RH_INT_records"
                status_rows.append(result)
                continue
        else:
            if not args.mutator_domain_root:
                result["status"] = "missing_mutator_domain"
                status_rows.append(result)
                continue
            candidates = list(args.mutator_domain_root.glob(f"{family}*.faa")) + list(args.mutator_domain_root.glob(f"{family}*.fa"))
            if not candidates:
                result["status"] = "missing_mutator_domain"
                status_rows.append(result)
                continue
            records = [(name, sequence, "transposase") for name, sequence in parse_fasta(candidates[0])]
            if len(records) < 3:
                result["status"] = "insufficient_mutator_domain_records"
                status_rows.append(result)
                continue
        by_gene: dict[str, list[tuple[str, str]]] = defaultdict(list)
        for name, sequence, gene in records:
            copy_id = name.split("#", 1)[0]
            by_gene[gene].append((copy_id, sequence))
        domain_trees = []
        domain_counts = {}
        support_warnings = []
        for gene, gene_records in sorted(by_gene.items()):
            if len(gene_records) < 3:
                continue
            domain_fasta = domain_dir / f"{family}.{gene}.faa"
            with domain_fasta.open("w", encoding="utf-8") as handle:
                for copy_id, sequence in gene_records:
                    handle.write(f">{copy_id}\n{sequence}\n")
            alignment = alignment_dir / f"{family}.{gene}.aln.faa"
            with alignment.open("w", encoding="utf-8") as handle, (args.outdir / "logs" / f"{family}.{gene}.mafft.stderr").open("w", encoding="utf-8") as err:
                proc = subprocess.run([args.mafft, "--auto", "--thread", str(args.threads), "--quiet", str(domain_fasta)], stdout=handle, stderr=err, text=True)
            if proc.returncode:
                continue
            prefix = tree_dir / f"{family}.{gene}"
            run([args.iqtree, "-st", "AA", "-s", str(alignment), "-m", "MFP", "-B", str(args.bootstrap), "--alrt", str(args.alrt), "-T", str(args.threads), "--prefix", str(prefix)], args.outdir / "logs" / f"{family}.{gene}.iqtree.log")
            tree_file = Path(str(prefix) + ".treefile")
            report_file = Path(str(prefix) + ".iqtree")
            report = report_file.read_text(encoding="utf-8", errors="replace") if report_file.is_file() else ""
            support_present = "SH-aLRT" in report and "UFBoot" in report and tree_file.is_file()
            if tree_file.is_file():
                domain_trees.append(str(tree_file))
                domain_counts[gene] = len(gene_records)
                if not support_present:
                    support_warnings.append(gene)
        result["domain_trees"] = ",".join(domain_trees)
        result["domain_record_count"] = ";".join(f"{k}:{v}" for k, v in sorted(domain_counts.items()))
        result["support_metrics_status"] = "present" if domain_trees and not support_warnings else ("warning" if domain_trees else "missing")
        result["status"] = "tree_inferred" if domain_trees and not support_warnings else ("tree_without_support_metrics" if domain_trees else "no_domain_tree")
        status_rows.append(result)
    fields = sorted({key for row in status_rows for key in row})
    write_tsv(args.outdir / "subfamily_tree_status.tsv", status_rows, fields)
    return status_rows


def parse_newick(text: str):
    text = text.strip().rstrip(";")
    cursor = 0

    def token_until(stops: str) -> str:
        nonlocal cursor
        start = cursor
        while cursor < len(text) and text[cursor] not in stops:
            cursor += 1
        return text[start:cursor]

    def parse_node():
        nonlocal cursor
        if cursor < len(text) and text[cursor] == "(":
            cursor += 1
            children = []
            while True:
                children.append(parse_node())
                if cursor >= len(text):
                    break
                if text[cursor] == ",":
                    cursor += 1
                    continue
                if text[cursor] == ")":
                    cursor += 1
                    break
            label = token_until(",()")
            token = label.split(":", 1)[0].strip()
            ufboot = None
            alrt = None
            if "/" in token:
                left, right = token.split("/", 1)
                try:
                    ufboot, alrt = float(left), float(right)
                except ValueError:
                    pass
            return {"children": children, "name": "", "ufboot": ufboot, "alrt": alrt}
        label = token_until(",()")
        return {"children": [], "name": label.split(":", 1)[0].strip(), "ufboot": None, "alrt": None}

    return parse_node()


def node_leaves(node) -> list[str]:
    if not node["children"]:
        return [node["name"]]
    leaves = []
    for child in node["children"]:
        leaves.extend(node_leaves(child))
    return leaves


def collect_supported_clades(node, family: str, domain: str, rows: list[dict[str, object]], ufboot_min: float = 95.0, alrt_min: float = 80.0):
    if not node["children"]:
        return
    for child in node["children"]:
        collect_supported_clades(child, family, domain, rows, ufboot_min, alrt_min)
    leaves = sorted(set(node_leaves(node)))
    if len(leaves) >= 3 and node.get("ufboot") is not None and node.get("alrt") is not None:
        if node["ufboot"] >= ufboot_min and node["alrt"] >= alrt_min:
            rows.append({
                "pan_family": family,
                "domain": domain,
                "ufboot": node["ufboot"],
                "sh_alrt": node["alrt"],
                "copy_count": len(leaves),
                "copies": ",".join(leaves),
                "copy_set": frozenset(leaves),
            })


def write_subfamily_calls(args, status_rows) -> None:
    clade_rows = []
    for row in status_rows:
        tree_paths = [x for x in str(row.get("domain_trees", "")).split(",") if x]
        for tree_path in tree_paths:
            tree_file = Path(tree_path)
            if not tree_file.is_file():
                continue
            parts = tree_file.stem.split(".")
            family = row["pan_family"]
            domain = parts[-1]
            tree = parse_newick(tree_file.read_text(encoding="utf-8", errors="replace"))
            all_leaves = set(node_leaves(tree))
            before = len(clade_rows)
            collect_supported_clades(tree, family, domain, clade_rows)
            # A node containing every copy is the family root, not a
            # biologically informative subfamily.
            clade_rows[before:] = [r for r in clade_rows[before:] if r["copy_set"] != all_leaves]
    raw_fields = ["pan_family", "domain", "ufboot", "sh_alrt", "copy_count", "copies"]
    write_tsv(args.outdir / "tree_supported_clades.tsv", [{k: v for k, v in r.items() if k != "copy_set"} for r in clade_rows], raw_fields)
    calls = []
    for family in sorted({r["pan_family"] for r in clade_rows}):
        family_clades = [r for r in clade_rows if r["pan_family"] == family]
        # Retain maximal supported clades; nested supported clades are useful
        # in the raw table but should not be counted twice as subfamilies.
        maximal = []
        for candidate in family_clades:
            if any(
                candidate["copy_set"] < other["copy_set"]
                and other["ufboot"] >= candidate["ufboot"]
                and other["sh_alrt"] >= candidate["sh_alrt"]
                for other in family_clades
            ):
                continue
            maximal.append(candidate)
        for index, candidate in enumerate(sorted(maximal, key=lambda r: (len(r["copy_set"]), r["domain"]), reverse=True), 1):
            matching_domains = set()
            for other in family_clades:
                if other["domain"] == candidate["domain"]:
                    continue
                intersection = len(candidate["copy_set"] & other["copy_set"])
                union = len(candidate["copy_set"] | other["copy_set"])
                if union and intersection / union >= 0.80:
                    matching_domains.add(other["domain"])
            domain_count = 1 + len(matching_domains)
            mode = next((str(r.get("tree_mode", "")) for r in status_rows if r["pan_family"] == family), "")
            calls.append({
                "pan_family": family,
                "subfamily_id": f"{family}_SF{index:02d}",
                "tree_mode": mode,
                "representative_domain": candidate["domain"],
                "supporting_domains": ",".join(sorted({candidate["domain"]} | matching_domains)),
                "supporting_domain_count": domain_count,
                "ufboot": candidate["ufboot"],
                "sh_alrt": candidate["sh_alrt"],
                "copy_count": candidate["copy_count"],
                "copies": candidate["copies"],
                "subfamily_status": "supported_multidomain" if domain_count >= 2 else "supported_single_domain",
            })
    call_fields = ["pan_family", "subfamily_id", "tree_mode", "representative_domain", "supporting_domains", "supporting_domain_count", "ufboot", "sh_alrt", "copy_count", "copies", "subfamily_status"]
    write_tsv(args.outdir / "subfamily_calls.tsv", calls, call_fields)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--pan-library", type=Path, required=True)
    parser.add_argument("--library-manifest", type=Path, required=True)
    parser.add_argument("--mask-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--bedtools", default="bedtools")
    parser.add_argument("--mmseqs", default="")
    parser.add_argument("--cdhit", default="cd-hit-est")
    parser.add_argument("--tesorter", default="TEsorter")
    parser.add_argument("--mafft", default="mafft")
    parser.add_argument("--iqtree", default="iqtree2")
    parser.add_argument("--mutator-domain-root", type=Path, default=None)
    parser.add_argument("--threads", type=int, default=24)
    parser.add_argument("--max-copies-per-target", type=int, default=50)
    parser.add_argument("--min-copy-bp", type=int, default=200)
    parser.add_argument("--min-copy-fraction", type=float, default=0.45)
    parser.add_argument("--max-copy-fraction", type=float, default=1.75)
    parser.add_argument("--merge-gap", type=int, default=100)
    parser.add_argument("--cluster-identity", type=float, default=0.80)
    parser.add_argument("--cluster-coverage", type=float, default=0.80)
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--alrt", type=int, default=1000)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    check_tool(args.bedtools, "bedtools")
    if args.mmseqs:
        check_tool(args.mmseqs, "MMseqs2")
    else:
        check_tool(args.cdhit, "CD-HIT-est or MMseqs2")
    check_tool(args.tesorter, "TEsorter")
    check_tool(args.mafft, "MAFFT")
    check_tool(args.iqtree, "IQ-TREE2")
    library_sha = verify_library(args.pan_library, args.library_manifest)
    library_lengths = read_library_lengths(args.pan_library)
    selected_rows = select_families(args.selection)
    missing = sorted({row["pan_family"] for row in selected_rows} - set(library_lengths))
    if missing:
        raise SystemExit(f"ERROR: selected families absent from frozen pan-TE FASTA: {', '.join(missing[:20])}")
    write_tsv(args.outdir / "run_manifest.tsv", [{"key": "library_sha256", "value": library_sha}, {"key": "selection_rows", "value": len(selected_rows)}], ["key", "value"])
    bed_dir, catalogue = write_beds_and_catalogue(args, selected_rows, library_lengths)
    family_dir = extract_sequences(args, bed_dir, catalogue)
    cluster_status = []
    for row in selected_rows:
        family = row["pan_family"]
        fasta = family_dir / f"{family}.fa"
        if fasta.is_file():
            method = cluster_family(args, fasta, family)
            cluster_status.append({"pan_family": family, "cluster_method": method})
    write_tsv(args.outdir / "cluster_status.tsv", cluster_status, ["pan_family", "cluster_method"])
    statuses = run_domain_trees(args, selected_rows, family_dir)
    write_subfamily_calls(args, statuses)
    inferred = sum(row.get("status") == "tree_inferred" for row in statuses)
    print(f"Completed copy extraction and tree attempts: {len(statuses)} families; trees inferred: {inferred}")


if __name__ == "__main__":
    main()
