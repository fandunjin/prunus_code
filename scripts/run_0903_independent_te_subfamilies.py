#!/usr/bin/env python3
"""Build Gypsy/Mutator cluster trees from independent EDTA outputs.

This is deliberately separate from the pan-TE workflow. It reads only the
independent input manifest made by ``run_0903_independent_edta_subgenome.py``.
Consensus sequences are clustered independently for Gypsy and Mutator-like
classes. For each selected cluster, intact EDTA intervals are preferred for
copy extraction; ordinary EDTA annotation intervals are used only when an
intact GFF3 is unavailable. A tree is a subfamily candidate until its domain
and structural evidence is reviewed.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
import threading
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


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


def parse_attrs(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for item in text.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = value
    return out


def check_tool(name: str) -> str:
    resolved = shutil.which(name)
    if not resolved:
        raise SystemExit(f"ERROR: required tool not found: {name}")
    return resolved


def family_id(header: str) -> str:
    return header.split("#", 1)[0]


def family_class(header: str) -> str:
    return header.split("#", 1)[1] if "#" in header else "unknown"


MUTATOR_CLASSES = {
    "DNA/DTM", "MITE/DTM", "DNA/MULE-MuDR",
    # newer EDTA classification vocabulary used by ppe/pcam/rch libraries
    "DNAnona/MULE", "DNAnona/MULEtir", "DNAauto/MULE",
}


def class_group(classification: str) -> str | None:
    if classification.startswith("LTR/Gypsy"):
        return "Gypsy"
    if classification in MUTATOR_CLASSES:
        return "Mutator"
    return None


_LOG_LOCK = threading.Lock()


def run_cmd(command: list[str], log: Path) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with _LOG_LOCK:
        with log.open("a", encoding="utf-8") as handle:
            handle.write("$ " + " ".join(map(str, command)) + "\n")
    with log.open("a", encoding="utf-8") as handle:
        proc = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, text=True)
    if proc.returncode:
        raise RuntimeError(f"command failed ({proc.returncode}); see {log}")


def write_fasta(path: Path, records: dict[str, str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, seq in records.items():
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i + 80] + "\n")


def load_consensus(rows: list[dict[str, str]], out: Path, min_len: int) -> tuple[dict[str, str], list[dict[str, object]]]:
    records: dict[str, str] = {}
    catalogue: list[dict[str, object]] = []
    seen_libraries: set[str] = set()
    for row in rows:
        library = Path(row["edta_library_fasta"])
        if not library.is_file():
            continue
        library_key = str(library.resolve())
        if library_key in seen_libraries:
            continue
        seen_libraries.add(library_key)
        owner = row.get("library_owner") or row["target"]
        for header, sequence in parse_fasta(library):
            fid = family_id(header)
            classification = family_class(header)
            group = class_group(classification)
            if group is None or len(sequence) < min_len:
                continue
            # A combined phased EDTA library is loaded once.  Its owner ID
            # prevents duplicate consensus records for PSP/PPA rows that
            # share the same physical library.
            record_id = f"{owner}__{fid}"
            records[record_id] = sequence.upper()
            catalogue.append({
                "record_id": record_id,
                "target": owner,
                "species": row["species"],
                "subgenome": row["subgenome"],
                "independent_family_id": fid,
                "classification": classification,
                "group": group,
                "length_bp": len(sequence),
                "library_fasta": str(library),
                "library_owner": owner,
            })
    write_fasta(out / "independent_gypsy_mutator_consensus.fa", records)
    write_tsv(out / "consensus_catalogue.tsv", catalogue, list(catalogue[0]) if catalogue else [
        "record_id", "target", "species", "subgenome", "independent_family_id",
        "classification", "group", "length_bp", "library_fasta",
    ])
    return records, catalogue


def parse_clusters(path: Path) -> dict[str, list[str]]:
    clusters: dict[str, list[str]] = defaultdict(list)
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                clusters[fields[0]].append(fields[1])
    return dict(clusters)


def _scan_gff_for_families(
    row: dict[str, str], gff: Path, wanted_families: set[str], max_copies: int,
) -> list[dict[str, object]]:
    seq_regex = re.compile(row["seq_regex"]) if row.get("seq_regex") else None
    seen = 0
    intervals: list[dict[str, object]] = []
    for line in gff.open(encoding="utf-8", errors="replace"):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            continue
        seqid, _src, _feature, start_s, end_s, _score, strand, _phase, attr_text = fields
        if seq_regex and not seq_regex.search(seqid):
            continue
        attrs = parse_attrs(attr_text)
        fid = attrs.get("Name", "")
        if fid not in wanted_families:
            continue
        try:
            start, end = int(start_s), int(end_s)
        except ValueError:
            continue
        if end < start:
            start, end = end, start
        seen += 1
        if seen > max_copies:
            break
        intervals.append({
            "seqid": seqid,
            "start0": start - 1,
            "end": end,
            "strand": strand if strand in {"+", "-"} else "+",
            "family_id": fid,
            "target": row["target"],
            "genome_fasta": row["annotation_genome_fasta"],
        })
    return intervals


def extract_intervals(
    row: dict[str, str], wanted_families: set[str], max_copies: int,
) -> tuple[list[dict[str, object]], str]:
    """Return BED records.

    The intact GFF3 is preferred because it contains structurally validated
    elements.  EDTA versions differ between targets: ppe/pcam intact GFF3 use
    coordinate strings instead of library family IDs, so families without an
    intact match fall back to the annotation GFF3 (searched per family, no
    double counting).
    """
    candidates: list[tuple[str, Path]] = []
    intact = Path(row["edta_intact_gff3"])
    ordinary = Path(row["edta_gff3"])
    if intact.is_file():
        candidates.append(("intact", intact))
    if ordinary.is_file():
        candidates.append(("annotation_fallback", ordinary))
    if not candidates:
        return [], "missing_gff3"

    intervals: list[dict[str, object]] = []
    remaining = set(wanted_families)
    sources: list[str] = []
    for source, gff in candidates:
        if not remaining:
            break
        found = _scan_gff_for_families(row, gff, remaining, max_copies)
        if found:
            sources.append(source)
            intervals.extend(found)
            remaining -= {str(item["family_id"]) for item in found}
    return intervals, ";".join(sources) if sources else "no_match"


def build_tree_for_cluster(
    group: str,
    cluster_id: str,
    members: list[str],
    records: dict[str, str],
    catalogue_by_record: dict[str, dict[str, object]],
    rows_by_target: dict[str, dict[str, str]],
    rows_by_library_owner: dict[str, list[dict[str, str]]],
    out: Path,
    mafft: str,
    iqtree: str,
    bedtools: str,
    max_copies: int,
    log: Path,
    iqtree_threads: str = "AUTO",
    resume: bool = False,
) -> dict[str, object]:
    cluster_dir = out / "trees" / group / cluster_id
    cluster_dir.mkdir(parents=True, exist_ok=True)
    tree_prefix = cluster_dir / "iqtree"
    if resume and Path(str(tree_prefix) + ".treefile").is_file():
        copies_fa_existing = cluster_dir / "copies.fa"
        n_copies = sum(1 for _ in parse_fasta(copies_fa_existing)) if copies_fa_existing.is_file() else 0
        return {
            "group": group,
            "cluster_id": cluster_id,
            "n_consensus_members": len(members),
            "n_copy_sequences": n_copies,
            "copy_source": "resumed",
            "status": "tree_built_resumed",
            "cluster_fasta": str(cluster_dir / "consensus.fa"),
            "copy_fasta": str(copies_fa_existing),
            "alignment": str(cluster_dir / "copies.aln.fa"),
            "tree": str(tree_prefix) + ".treefile",
        }
    wanted_by_owner: dict[str, set[str]] = defaultdict(set)
    for member in members:
        meta = catalogue_by_record[member]
        owner = str(meta.get("library_owner") or meta["target"])
        wanted_by_owner[owner].add(str(meta["independent_family_id"]))

    intervals: list[dict[str, object]] = []
    sources: set[str] = set()
    # For a combined PSP/PPA library, apply each consensus family to every
    # manifest row sharing that library owner so all subgenome partitions are
    # represented despite loading the physical library only once.
    for owner, wanted in wanted_by_owner.items():
        for row in rows_by_library_owner.get(owner, [rows_by_target.get(owner)]):
            if row is None:
                continue
            extracted, source = extract_intervals(row, wanted, max_copies)
            intervals.extend(extracted)
            sources.add(source)
    bed = cluster_dir / "copies.bed"
    with bed.open("w", encoding="utf-8") as handle:
        for i, item in enumerate(intervals, 1):
            copy_id = f"{cluster_id}__copy{i}__{item['target']}__{item['family_id']}"
            handle.write("\t".join(map(str, [
                item["seqid"], item["start0"], item["end"], copy_id, ".", item["strand"]
            ])) + "\n")
    copies_fa = cluster_dir / "copies.fa"
    if intervals:
        # One bedtools call per cluster keeps the extraction reproducible and
        # supports PPA, whose annotation and denominator FASTAs differ.
        by_genome: dict[str, list[dict[str, object]]] = defaultdict(list)
        for item in intervals:
            by_genome[str(item["genome_fasta"])].append(item)
        pieces: list[str] = []
        for idx, (genome, items) in enumerate(by_genome.items(), 1):
            subbed = cluster_dir / f"copies_{idx}.bed"
            with subbed.open("w", encoding="utf-8") as handle:
                for i, item in enumerate(items, 1):
                    copy_id = f"{cluster_id}__copy{i}__{item['target']}__{item['family_id']}"
                    handle.write("\t".join(map(str, [
                        item["seqid"], item["start0"], item["end"], copy_id, ".", item["strand"]
                    ])) + "\n")
            piece = cluster_dir / f"copies_{idx}.fa"
            run_cmd([bedtools, "getfasta", "-s", "-name", "-fi", genome, "-bed", str(subbed), "-fo", str(piece)], log)
            pieces.append(piece.read_text(encoding="utf-8"))
        copies_fa.write_text("".join(pieces), encoding="utf-8")
    n_copies = sum(1 for _ in parse_fasta(copies_fa)) if copies_fa.is_file() else 0
    alignment = cluster_dir / "copies.aln.fa"
    tree_prefix = cluster_dir / "iqtree"
    status = "too_few_copies"
    if n_copies >= 4:
        with alignment.open("w", encoding="utf-8") as handle:
            proc = subprocess.run([mafft, "--auto", str(copies_fa)], stdout=handle, stderr=subprocess.PIPE, text=True)
        if proc.returncode:
            raise RuntimeError(f"MAFFT failed for {cluster_id}: {proc.stderr[-500:]}")
        run_cmd([iqtree, "-s", str(alignment), "-m", "MFP", "-B", "1000", "--alrt", "1000", "-T", str(iqtree_threads), "--prefix", str(tree_prefix)], log)
        status = "tree_built"
    return {
        "group": group,
        "cluster_id": cluster_id,
        "n_consensus_members": len(members),
        "n_copy_sequences": n_copies,
        "copy_source": ";".join(sorted(sources)),
        "status": status,
        "cluster_fasta": str(cluster_dir / "consensus.fa"),
        "copy_fasta": str(copies_fa),
        "alignment": str(alignment),
        "tree": str(tree_prefix) + ".treefile",
    }


def run(args: argparse.Namespace) -> int:
    root = Path(args.project_root).resolve()
    out = Path(args.output_dir).resolve() if args.output_dir else root / "0903" / "TE" / "independent_edta" / "subfamilies"
    out.mkdir(parents=True, exist_ok=True)
    manifest = Path(args.manifest).resolve() if args.manifest else root / "0903" / "TE" / "independent_edta" / "input_manifest.tsv"
    if not manifest.is_file():
        raise SystemExit(f"ERROR: missing manifest: {manifest}")
    rows = read_tsv(manifest)
    availability = []
    for row in rows:
        availability.append({
            "target": row["target"],
            "group": row["species"],
            "library_ok": int(Path(row["edta_library_fasta"]).is_file()),
            "intact_gff_ok": int(Path(row["edta_intact_gff3"]).is_file()),
            "annotation_gff_ok": int(Path(row["edta_gff3"]).is_file()),
            "library": row["edta_library_fasta"],
            "intact_gff": row["edta_intact_gff3"],
        })
    write_tsv(out / "tree_input_availability.tsv", availability, list(availability[0]))
    if args.require_all_libraries and any(not x["library_ok"] for x in availability):
        raise SystemExit("ERROR: --require-all-libraries requested, but one or more independent TElib.fa files are missing")

    check_tool(args.mmseqs)
    bedtools = check_tool(args.bedtools)
    mafft = check_tool(args.mafft)
    iqtree = check_tool(args.iqtree)
    log = out / "logs" / "subfamily_trees.log"
    records, catalogue = load_consensus(rows, out, args.min_consensus_length)
    if not records:
        raise SystemExit("ERROR: no Gypsy or Mutator consensus records were found in independent libraries")
    catalogue_by_record = {str(r["record_id"]): r for r in catalogue}
    rows_by_target = {r["target"]: r for r in rows}
    rows_by_library_owner: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        rows_by_library_owner[row.get("library_owner") or row["target"]].append(row)
    all_status: list[dict[str, object]] = []
    for group in ("Gypsy", "Mutator"):
        group_records = {k: v for k, v in records.items() if catalogue_by_record[k]["group"] == group}
        if not group_records:
            continue
        group_fa = out / f"{group.lower()}_consensus.fa"
        write_fasta(group_fa, group_records)
        prefix = out / "mmseqs" / group.lower()
        prefix.parent.mkdir(parents=True, exist_ok=True)
        tmp = out / "mmseqs" / f"tmp_{group.lower()}"
        run_cmd([
            args.mmseqs, "easy-cluster", str(group_fa), str(prefix), str(tmp),
            "--min-seq-id", str(args.min_seq_id), "-c", str(args.min_coverage),
            "--cov-mode", "1", "--threads", str(args.threads),
        ], log)
        cluster_tsv = Path(str(prefix) + "_cluster.tsv")
        clusters = parse_clusters(cluster_tsv)
        map_rows = []
        for cluster_id, members in clusters.items():
            for member in members:
                map_rows.append({
                    "group": group,
                    "cluster_id": cluster_id,
                    "record_id": member,
                    "target": catalogue_by_record[member]["target"],
                    "subgenome": catalogue_by_record[member]["subgenome"],
                    "independent_family_id": catalogue_by_record[member]["independent_family_id"],
                    "classification": catalogue_by_record[member]["classification"],
                })
        write_tsv(out / f"{group.lower()}_clusters.tsv", map_rows, list(map_rows[0]) if map_rows else [
            "group", "cluster_id", "record_id", "target", "subgenome",
            "independent_family_id", "classification",
        ])
        ordered_clusters = [
            (cid, members) for cid, members in
            sorted(clusters.items(), key=lambda kv: (-len(kv[1]), kv[0]))
            if len(members) >= args.min_cluster_members
        ][:args.max_clusters]
        def run_cluster(cluster_id: str, members: list[str]) -> dict[str, object]:
            cluster_fa = out / "trees" / group / cluster_id / "consensus.fa"
            cluster_fa.parent.mkdir(parents=True, exist_ok=True)
            write_fasta(cluster_fa, {m: records[m] for m in members})
            return build_tree_for_cluster(
                group, cluster_id, members, records, catalogue_by_record,
                rows_by_target, rows_by_library_owner, out, mafft, iqtree,
                bedtools, args.max_copies, log,
                iqtree_threads=args.iqtree_threads, resume=args.resume,
            )

        def run_cluster_safe(cluster_id: str, members: list[str]) -> dict[str, object]:
            try:
                return run_cluster(cluster_id, members)
            except Exception as exc:
                print(f"[warn] {group} cluster {cluster_id} failed: {exc}", flush=True)
                return {
                    "group": group,
                    "cluster_id": cluster_id,
                    "n_consensus_members": len(members),
                    "n_copy_sequences": "",
                    "copy_source": "",
                    "status": f"failed: {exc}"[:300],
                    "cluster_fasta": "",
                    "copy_fasta": "",
                    "alignment": "",
                    "tree": "",
                }

        if args.workers > 1 and len(ordered_clusters) > 1:
            print(f"[{group}] {len(ordered_clusters)} trees, {args.workers} concurrent workers, "
                  f"{args.iqtree_threads} IQ-TREE threads each", flush=True)
            with ThreadPoolExecutor(max_workers=args.workers) as pool:
                futures = [pool.submit(run_cluster_safe, cid, mem) for cid, mem in ordered_clusters]
                for fut in as_completed(futures):
                    all_status.append(fut.result())
        else:
            for cluster_id, members in ordered_clusters:
                all_status.append(run_cluster_safe(cluster_id, members))
    write_tsv(out / "subfamily_tree_status.tsv", all_status, list(all_status[0]) if all_status else [
        "group", "cluster_id", "n_consensus_members", "n_copy_sequences", "copy_source",
        "status", "cluster_fasta", "copy_fasta", "alignment", "tree",
    ])
    print(f"Independent consensus and tree outputs: {out}")
    print(f"Tree candidates: {len(all_status)}")
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--manifest", default="")
    parser.add_argument("--output-dir", default="")
    parser.add_argument("--mmseqs", default="mmseqs")
    parser.add_argument("--bedtools", default="bedtools")
    parser.add_argument("--mafft", default="mafft")
    parser.add_argument("--iqtree", default="iqtree2")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--workers", type=int, default=1,
                        help="Number of clusters processed concurrently (parallel independent IQ-TREE runs)")
    parser.add_argument("--iqtree-threads", default="AUTO",
                        help="Threads per IQ-TREE run (e.g. 6 when using 4 workers); default AUTO")
    parser.add_argument("--resume", action="store_true",
                        help="Skip clusters that already have iqtree.treefile (resume after a restart)")
    parser.add_argument("--min-seq-id", type=float, default=0.80)
    parser.add_argument("--min-coverage", type=float, default=0.80)
    parser.add_argument("--min-consensus-length", type=int, default=200)
    parser.add_argument("--min-cluster-members", type=int, default=1,
                        help="Only build trees for clusters with at least this many consensus members")
    parser.add_argument("--max-copies", type=int, default=200)
    parser.add_argument("--max-clusters", type=int, default=200)
    parser.add_argument("--require-all-libraries", action="store_true")
    args = parser.parse_args()
    raise SystemExit(run(args))


if __name__ == "__main__":
    main()
