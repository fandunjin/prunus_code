#!/usr/bin/env python3
"""Run synteny-based Ks/4DTV/KaKs analysis from JCVI gene pairs.

This script uses JCVI syntenic gene-pair tables as the pair source, aligns each
protein pair with MAFFT, back-translates with PAL2NAL, computes pairwise 4DTV,
and runs KaKs_Calculator on concatenated AXT files per comparison.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shlex
import shutil
import statistics
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path


FOURFOLD_PREFIXES = {"GC", "GG", "CC", "AC", "GT", "TC", "CG", "CT"}
VALID_BASES = set("ACGT")
STOP_CODONS = {"TAA", "TAG", "TGA"}


@dataclass(frozen=True)
class PairRecord:
    pair_id: str
    comparison: str
    comparison_type: str
    query_label: str
    ref_label: str
    block_id: str
    query_gene: str
    ref_gene: str
    query_seqid: str
    ref_seqid: str
    block_support: int
    chromosome_pair_support: int
    chromosome_filter_rule: str


@dataclass
class AlignmentResult:
    pair: PairRecord
    status: str
    message: str = ""
    aligned_codons: int = 0
    ungapped_codons: int = 0
    fourfold_sites: int = 0
    fourfold_differences: int = 0
    p_4dtv: str = ""
    axt_block: str = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", required=True)
    parser.add_argument("--jcvi-dir", required=True)
    parser.add_argument(
        "--pair-table",
        default="",
        help="Optional JCVI-syntenic-pair TSV. Defaults to <jcvi-dir>/summary/jcvi_syntenic_pairs.tsv.",
    )
    parser.add_argument(
        "--include-comparison",
        action="append",
        default=[],
        help="Restrict to named comparisons from the pair table; repeatable.",
    )
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--mafft-cmd", default="mafft")
    parser.add_argument("--pal2nal-cmd", default="pal2nal.pl")
    parser.add_argument("--kaks-cmd", default="KaKs_Calculator")
    parser.add_argument("--kaks-method", default="YN")
    parser.add_argument("--min-aa", type=int, default=50)
    parser.add_argument("--min-codons", type=int, default=50)
    parser.add_argument("--min-4d-sites", type=int, default=4)
    parser.add_argument(
        "--pair-mode",
        choices=["block_representative", "all_pairs", "collinear_all_pairs"],
        default="block_representative",
        help=(
            "Use one pair per JCVI block, every raw pair, or the publication-oriented "
            "collinear filter followed by all retained pairs."
        ),
    )
    parser.add_argument(
        "--min-chrom-pair-anchors",
        type=int,
        default=75,
        help="Minimum aggregate anchors for an A-B chromosome-segment relationship.",
    )
    parser.add_argument(
        "--min-chrom-pair-fraction",
        type=float,
        default=0.005,
        help="Minimum fraction of comparison-wide anchors for an A-B chromosome relationship.",
    )
    parser.add_argument(
        "--min-reference-block-anchors",
        type=int,
        default=20,
        help="Minimum anchors in a syntenic block for pcam/ppe/gtr reference comparisons.",
    )
    parser.add_argument("--max-pairs-per-comparison", type=int, default=0)
    parser.add_argument("--keep-temp", type=int, default=0)
    return parser.parse_args()


def safe_name(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_.:-]+", "_", value)
    return value.strip("_") or "pair"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    current: str | None = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current is not None:
                    records[current] = "".join(chunks).upper()
                current = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if current is not None:
            records[current] = "".join(chunks).upper()
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for idx in range(0, len(seq), 80):
                handle.write(seq[idx : idx + 80] + "\n")


def read_fasta_ordered(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    current: str | None = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current is not None:
                    records.append((current, "".join(chunks).upper()))
                current = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if current is not None:
            records.append((current, "".join(chunks).upper()))
    return records


def normalize_pep(seq: str) -> str:
    return seq.upper().replace("*", "").replace(" ", "").replace("\t", "")


def normalize_cds(seq: str) -> tuple[str, str]:
    seq = seq.upper().replace("U", "T")
    seq = re.sub(r"[^ACGTN-]", "", seq).replace("-", "")
    if len(seq) < 3:
        return "", "cds_too_short"
    if len(seq) % 3 != 0:
        return "", "cds_length_not_multiple_of_3"
    codons = [seq[i : i + 3] for i in range(0, len(seq), 3)]
    if codons and codons[-1] in STOP_CODONS:
        codons = codons[:-1]
    if any(codon in STOP_CODONS for codon in codons):
        return "", "internal_stop_codon"
    seq = "".join(codons)
    if len(seq) < 3:
        return "", "cds_too_short_after_stop_trim"
    return seq, "ok"


def comparison_type(row: dict[str, str]) -> str:
    if row.get("comparison_type"):
        # The Figure 2 all-taxon anchor table uses ``gtr_reference`` for the
        # same external-reference class named ``cross_species_reference`` in
        # the original JCVI summary.
        return "cross_species_reference" if row["comparison_type"] == "gtr_reference" else row["comparison_type"]
    q = row.get("query_label", "")
    r = row.get("ref_label", "")
    if q.startswith("pzi_hap") and r.startswith("pzi_hap"):
        return "internal_haplotype"
    return "external_reference"


def chromosome_key(seqid: str) -> str:
    match = re.match(r"^C(\d{2})(?:_|$)", seqid or "")
    return match.group(1) if match else ""


def load_pairs(
    pair_table: Path,
    include_comparisons: set[str],
    pair_mode: str,
    max_pairs_per_comparison: int,
    min_chrom_pair_anchors: int,
    min_chrom_pair_fraction: float,
    min_reference_block_anchors: int,
) -> tuple[list[PairRecord], list[dict[str, object]]]:
    if not pair_table.is_file():
        raise SystemExit(f"Missing JCVI syntenic pair table: {pair_table}")
    rows = read_tsv(pair_table)
    if include_comparisons:
        rows = [row for row in rows if row.get("comparison", "") in include_comparisons]
    comparison_totals: dict[str, int] = {}
    chromosome_pair_counts: dict[tuple[str, str, str], int] = {}
    chromosome_pair_blocks: dict[tuple[str, str, str], set[str]] = {}
    block_counts: dict[tuple[str, str], int] = {}
    comparison_types: dict[str, str] = {}
    for row in rows:
        comparison = row.get("comparison", "")
        if not comparison:
            continue
        comparison_totals[comparison] = comparison_totals.get(comparison, 0) + 1
        pair_key = (comparison, row.get("query_seqid", ""), row.get("ref_seqid", ""))
        chromosome_pair_counts[pair_key] = chromosome_pair_counts.get(pair_key, 0) + 1
        chromosome_pair_blocks.setdefault(pair_key, set()).add(row.get("block_id", ""))
        block_key = (comparison, row.get("block_id", ""))
        block_counts[block_key] = block_counts.get(block_key, 0) + 1
        comparison_types[comparison] = comparison_type(row)

    chromosome_filter_rows: list[dict[str, object]] = []
    accepted_chromosome_pairs: set[tuple[str, str, str]] = set()
    for (comparison, query_seqid, ref_seqid), support in sorted(chromosome_pair_counts.items()):
        comp_type = comparison_types.get(comparison, "")
        threshold = max(
            min_chrom_pair_anchors,
            int(math.ceil(comparison_totals.get(comparison, 0) * min_chrom_pair_fraction)),
        )
        query_key = chromosome_key(query_seqid)
        ref_key = chromosome_key(ref_seqid)
        if comp_type == "cross_species_reference":
            accepted = any(
                block_counts.get((comparison, block_id), 0) >= min_reference_block_anchors
                for block_id in chromosome_pair_blocks.get((comparison, query_seqid, ref_seqid), set())
            )
            rule = "syntenic_block_minimum_support"
            audit_threshold: object = min_reference_block_anchors
        elif comp_type == "self_synteny":
            # These pairs have already been restricted by the i-ADHoRe
            # adapter to real anchorpoints between different chromosomes.
            accepted = True
            rule = "iadhore_cross_chromosome_real_anchorpoint"
            audit_threshold = "pre_filtered"
        elif comp_type.startswith("cross_species_"):
            accepted = bool(query_key and query_key == ref_key)
            rule = "same_chromosome_key"
            audit_threshold = "NA"
        elif comp_type == "within_species_homeolog":
            accepted = support >= threshold
            rule = "major_collinear_chromosome_segment"
            audit_threshold = threshold
        else:
            accepted = bool(query_key and query_key == ref_key)
            rule = "same_chromosome_key_fallback"
            audit_threshold = "NA"
        if accepted:
            accepted_chromosome_pairs.add((comparison, query_seqid, ref_seqid))
        chromosome_filter_rows.append(
            {
                "comparison": comparison,
                "comparison_type": comp_type,
                "query_seqid": query_seqid,
                "ref_seqid": ref_seqid,
                "query_chromosome_key": query_key,
                "ref_chromosome_key": ref_key,
                "anchor_support": support,
                "support_threshold": audit_threshold,
                "filter_rule": rule,
                "accepted": "PASS" if accepted else "FAIL",
            }
        )

    seen: set[tuple[str, str, str]] = set()
    chosen_blocks: set[tuple[str, str]] = set()
    counts: dict[str, int] = {}
    pairs: list[PairRecord] = []
    for row in rows:
        comparison = row.get("comparison", "")
        qgene = row.get("query_gene", "")
        rgene = row.get("ref_gene", "")
        block_id = row.get("block_id", "")
        if not comparison or not qgene or not rgene:
            continue
        chromosome_pair = (comparison, row.get("query_seqid", ""), row.get("ref_seqid", ""))
        block_support = block_counts.get((comparison, block_id), 0)
        comp_type = comparison_type(row)
        if pair_mode == "collinear_all_pairs":
            if comp_type == "cross_species_reference":
                if block_support < min_reference_block_anchors:
                    continue
            elif comp_type == "self_synteny":
                # The source table is an already filtered i-ADHoRe result.
                pass
            elif chromosome_pair not in accepted_chromosome_pairs:
                continue
        key = (comparison, qgene, rgene)
        if key in seen:
            continue
        if pair_mode == "block_representative":
            block_key = (comparison, block_id)
            if block_key in chosen_blocks:
                continue
            chosen_blocks.add(block_key)
        seen.add(key)
        counts.setdefault(comparison, 0)
        if max_pairs_per_comparison > 0 and counts[comparison] >= max_pairs_per_comparison:
            continue
        counts[comparison] += 1
        pair_id = f"{safe_name(comparison)}__pair{counts[comparison]:06d}"
        pairs.append(
            PairRecord(
                pair_id=pair_id,
                comparison=comparison,
                comparison_type=comparison_type(row),
                query_label=row.get("query_label", ""),
                ref_label=row.get("ref_label", ""),
                block_id=block_id,
                query_gene=qgene,
                ref_gene=rgene,
                query_seqid=row.get("query_seqid", ""),
                ref_seqid=row.get("ref_seqid", ""),
                block_support=block_support,
                chromosome_pair_support=chromosome_pair_counts.get(chromosome_pair, 0),
                chromosome_filter_rule=(
                    "not_applied"
                    if pair_mode != "collinear_all_pairs"
                    else (
                        "syntenic_block_minimum_support"
                        if comp_type == "cross_species_reference"
                        else (
                            "iadhore_cross_chromosome_real_anchorpoint"
                            if comp_type == "self_synteny"
                            else (
                            "same_chromosome_key"
                            if comp_type.startswith("cross_species_")
                            else "major_collinear_chromosome_segment"
                            )
                        )
                    )
                ),
            )
        )
    return pairs, chromosome_filter_rows


def label_fasta_paths(analysis_dir: Path, jcvi_dir: Path, label: str) -> tuple[Path, Path]:
    jcvi_input = jcvi_dir / "inputs" / label
    pep = jcvi_input / f"{label}.pep"
    cds = jcvi_input / f"{label}.cds"
    if pep.is_file() and cds.is_file():
        return pep, cds
    canonical = analysis_dir / "canonical" / label
    pep = canonical / "primary.pep.fa"
    cds = canonical / "primary.cds.fa"
    return pep, cds


def load_sequences(analysis_dir: Path, jcvi_dir: Path, labels: set[str]) -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    pep_by_label: dict[str, dict[str, str]] = {}
    cds_by_label: dict[str, dict[str, str]] = {}
    for label in sorted(labels):
        pep_path, cds_path = label_fasta_paths(analysis_dir, jcvi_dir, label)
        if not pep_path.is_file() or not cds_path.is_file():
            raise SystemExit(f"Missing pep/cds for {label}: {pep_path} {cds_path}")
        pep_by_label[label] = read_fasta(pep_path)
        cds_by_label[label] = read_fasta(cds_path)
    return pep_by_label, cds_by_label


def run_command(cmd_text: str, cwd: Path, stdout_path: Path | None = None) -> tuple[int, str]:
    cmd = shlex.split(cmd_text)
    if not cmd:
        return 127, "empty command"
    if stdout_path is None:
        proc = subprocess.run(cmd, cwd=str(cwd), text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        with stdout_path.open("w", encoding="utf-8") as out:
            proc = subprocess.run(cmd, cwd=str(cwd), text=True, stdout=out, stderr=subprocess.PIPE)
    return proc.returncode, proc.stderr.strip()


def parse_codon_alignment(path: Path) -> tuple[str, str] | None:
    records = read_fasta_ordered(path)
    if len(records) != 2:
        return None
    seq_a = records[0][1].upper()
    seq_b = records[1][1].upper()
    if len(seq_a) != len(seq_b) or len(seq_a) % 3 != 0:
        return None
    return seq_a, seq_b


def pairwise_4dtv(seq_a: str, seq_b: str) -> tuple[int, int, str]:
    sites = 0
    diffs = 0
    for idx in range(0, len(seq_a), 3):
        codon_a = seq_a[idx : idx + 3]
        codon_b = seq_b[idx : idx + 3]
        if "-" in codon_a or "-" in codon_b or "N" in codon_a or "N" in codon_b:
            continue
        if any(base not in VALID_BASES for base in codon_a + codon_b):
            continue
        if codon_a[:2] not in FOURFOLD_PREFIXES or codon_b[:2] not in FOURFOLD_PREFIXES:
            continue
        sites += 1
        if codon_a[2] != codon_b[2]:
            diffs += 1
    if sites == 0:
        return sites, diffs, ""
    return sites, diffs, f"{diffs / sites:.8f}"


def process_pair(
    pair: PairRecord,
    pep_by_label: dict[str, dict[str, str]],
    cds_by_label: dict[str, dict[str, str]],
    mafft_cmd: str,
    pal2nal_cmd: str,
    min_aa: int,
    min_codons: int,
    min_4d_sites: int,
    tmp_root: Path,
    keep_temp: bool,
) -> AlignmentResult:
    qpep = pep_by_label.get(pair.query_label, {}).get(pair.query_gene)
    rpep = pep_by_label.get(pair.ref_label, {}).get(pair.ref_gene)
    qcds_raw = cds_by_label.get(pair.query_label, {}).get(pair.query_gene)
    rcds_raw = cds_by_label.get(pair.ref_label, {}).get(pair.ref_gene)
    if qpep is None or rpep is None:
        return AlignmentResult(pair, "skip_missing_pep")
    if qcds_raw is None or rcds_raw is None:
        return AlignmentResult(pair, "skip_missing_cds")
    qpep = normalize_pep(qpep)
    rpep = normalize_pep(rpep)
    if len(qpep) < min_aa or len(rpep) < min_aa:
        return AlignmentResult(pair, "skip_pep_too_short")
    qcds, qstatus = normalize_cds(qcds_raw)
    rcds, rstatus = normalize_cds(rcds_raw)
    if qstatus != "ok":
        return AlignmentResult(pair, f"skip_query_{qstatus}")
    if rstatus != "ok":
        return AlignmentResult(pair, f"skip_ref_{rstatus}")

    work_dir = Path(tempfile.mkdtemp(prefix=pair.pair_id + "_", dir=str(tmp_root)))
    try:
        pep_in = work_dir / "pair.pep.fa"
        cds_in = work_dir / "pair.cds.fa"
        pep_aln = work_dir / "pair.pep.aln.fa"
        codon_aln = work_dir / "pair.codon.aln.fa"
        write_fasta(pep_in, [(pair.query_gene, qpep), (pair.ref_gene, rpep)])
        write_fasta(cds_in, [(pair.query_gene, qcds), (pair.ref_gene, rcds)])

        code, err = run_command(f"{mafft_cmd} --anysymbol --quiet {pep_in.name}", work_dir, stdout_path=pep_aln)
        if code != 0:
            return AlignmentResult(pair, "skip_mafft_failed", err[:300])

        code, err = run_command(
            f"{pal2nal_cmd} {pep_aln.name} {cds_in.name} -output fasta -nogap",
            work_dir,
            stdout_path=codon_aln,
        )
        if code != 0:
            return AlignmentResult(pair, "skip_pal2nal_failed", err[:300])

        parsed = parse_codon_alignment(codon_aln)
        if parsed is None:
            return AlignmentResult(pair, "skip_bad_codon_alignment")
        seq_a, seq_b = parsed
        aligned_codons = len(seq_a) // 3
        if aligned_codons < min_codons:
            return AlignmentResult(pair, "skip_aligned_codon_too_short", aligned_codons=aligned_codons)
        ungapped_codons = 0
        for idx in range(0, len(seq_a), 3):
            if "-" not in seq_a[idx : idx + 3] and "-" not in seq_b[idx : idx + 3]:
                ungapped_codons += 1
        fourfold_sites, fourfold_diffs, p_4dtv = pairwise_4dtv(seq_a, seq_b)
        if fourfold_sites < min_4d_sites:
            return AlignmentResult(
                pair,
                "skip_too_few_4d_sites",
                aligned_codons=aligned_codons,
                ungapped_codons=ungapped_codons,
                fourfold_sites=fourfold_sites,
                fourfold_differences=fourfold_diffs,
                p_4dtv=p_4dtv,
            )
        axt_block = f"{pair.pair_id}\n{seq_a}\n{seq_b}\n\n"
        return AlignmentResult(
            pair,
            "used",
            aligned_codons=aligned_codons,
            ungapped_codons=ungapped_codons,
            fourfold_sites=fourfold_sites,
            fourfold_differences=fourfold_diffs,
            p_4dtv=p_4dtv,
            axt_block=axt_block,
        )
    finally:
        if not keep_temp:
            shutil.rmtree(work_dir, ignore_errors=True)


def result_row(result: AlignmentResult) -> dict[str, object]:
    p = result.pair
    return {
        "pair_id": p.pair_id,
        "comparison": p.comparison,
        "comparison_type": p.comparison_type,
        "query_label": p.query_label,
        "ref_label": p.ref_label,
        "block_id": p.block_id,
        "query_gene": p.query_gene,
        "ref_gene": p.ref_gene,
        "query_seqid": p.query_seqid,
        "ref_seqid": p.ref_seqid,
        "block_support": p.block_support,
        "chromosome_pair_support": p.chromosome_pair_support,
        "chromosome_filter_rule": p.chromosome_filter_rule,
        "status": result.status,
        "message": result.message,
        "aligned_codons": result.aligned_codons,
        "ungapped_codons": result.ungapped_codons,
        "fourfold_sites": result.fourfold_sites,
        "fourfold_differences": result.fourfold_differences,
        "p_4dtv": result.p_4dtv,
    }


def run_kaks(kaks_cmd: str, method: str, axt: Path, out: Path) -> tuple[int, str]:
    cmd = f"{kaks_cmd} -i {axt} -o {out} -m {method}"
    return run_command(cmd, cwd=axt.parent)


def parse_float(value: str | None) -> float | None:
    if value is None:
        return None
    value = value.strip()
    if not value or value.upper() in {"NA", "NAN", "NULL", "NONE"}:
        return None
    try:
        val = float(value)
    except ValueError:
        return None
    if math.isnan(val) or math.isinf(val):
        return None
    return val


def parse_kaks_table(path: Path) -> dict[str, dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return {}
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return {}
        first = reader.fieldnames[0]
        out: dict[str, dict[str, str]] = {}
        for row in reader:
            seq = row.get("Sequence") or row.get("sequence") or row.get(first) or ""
            if seq:
                out[seq] = row
        return out


def merge_kaks(
    results: list[AlignmentResult],
    kaks_by_pair: dict[str, dict[str, str]],
    out_path: Path,
) -> list[dict[str, object]]:
    fieldnames = [
        "pair_id",
        "comparison",
        "comparison_type",
        "query_label",
        "ref_label",
        "block_id",
        "query_gene",
        "ref_gene",
        "query_seqid",
        "ref_seqid",
        "block_support",
        "chromosome_pair_support",
        "chromosome_filter_rule",
        "aligned_codons",
        "ungapped_codons",
        "fourfold_sites",
        "fourfold_differences",
        "p_4dtv",
        "kaks_status",
        "method",
        "Ka",
        "Ks",
        "Ka_Ks",
        "P_value",
    ]
    rows: list[dict[str, object]] = []
    for result in results:
        if result.status != "used":
            continue
        p = result.pair
        krow = kaks_by_pair.get(p.pair_id, {})
        ka = krow.get("Ka", "")
        ks = krow.get("Ks", "")
        ratio = krow.get("Ka/Ks", krow.get("Ka_Ks", ""))
        p_value = krow.get("P-Value(Fisher)", krow.get("P_value", ""))
        rows.append(
            {
                "pair_id": p.pair_id,
                "comparison": p.comparison,
                "comparison_type": p.comparison_type,
                "query_label": p.query_label,
                "ref_label": p.ref_label,
                "block_id": p.block_id,
                "query_gene": p.query_gene,
                "ref_gene": p.ref_gene,
                "query_seqid": p.query_seqid,
                "ref_seqid": p.ref_seqid,
                "block_support": p.block_support,
                "chromosome_pair_support": p.chromosome_pair_support,
                "chromosome_filter_rule": p.chromosome_filter_rule,
                "aligned_codons": result.aligned_codons,
                "ungapped_codons": result.ungapped_codons,
                "fourfold_sites": result.fourfold_sites,
                "fourfold_differences": result.fourfold_differences,
                "p_4dtv": result.p_4dtv,
                "kaks_status": "ok" if krow else "missing_kaks",
                "method": krow.get("Method", ""),
                "Ka": ka,
                "Ks": ks,
                "Ka_Ks": ratio,
                "P_value": p_value,
            }
        )
    write_tsv(out_path, rows, fieldnames)
    return rows


def median(values: list[float]) -> str:
    if not values:
        return ""
    return f"{statistics.median(values):.8f}"


def mean(values: list[float]) -> str:
    if not values:
        return ""
    return f"{statistics.mean(values):.8f}"


def write_summary(outdir: Path, pair_counts: dict[str, int], align_rows: list[dict[str, object]], merged_rows: list[dict[str, object]]) -> None:
    align_by_comp: dict[str, list[dict[str, object]]] = {}
    merged_by_comp: dict[str, list[dict[str, object]]] = {}
    for row in align_rows:
        align_by_comp.setdefault(str(row["comparison"]), []).append(row)
    for row in merged_rows:
        merged_by_comp.setdefault(str(row["comparison"]), []).append(row)
    comps = sorted(set(pair_counts) | set(align_by_comp) | set(merged_by_comp))
    rows: list[dict[str, object]] = []
    for comp in comps:
        arows = align_by_comp.get(comp, [])
        mrows = merged_by_comp.get(comp, [])
        used = [r for r in arows if r["status"] == "used"]
        ks = [x for x in (parse_float(str(r.get("Ks", ""))) for r in mrows) if x is not None]
        ka = [x for x in (parse_float(str(r.get("Ka", ""))) for r in mrows) if x is not None]
        ratios = [x for x in (parse_float(str(r.get("Ka_Ks", ""))) for r in mrows) if x is not None]
        tv = [x for x in (parse_float(str(r.get("p_4dtv", ""))) for r in mrows) if x is not None]
        rows.append(
            {
                "comparison": comp,
                "input_pairs": pair_counts.get(comp, 0),
                "used_alignments": len(used),
                "kaks_rows": len(mrows),
                "ks_median": median(ks),
                "ks_mean": mean(ks),
                "ka_median": median(ka),
                "kaks_ratio_median": median(ratios),
                "p_4dtv_median": median(tv),
                "ks_gt_3_count": sum(1 for x in ks if x > 3),
                "p_4dtv_gt_0_75_count": sum(1 for x in tv if x > 0.75),
            }
        )
    write_tsv(
        outdir / "summary" / "wgd_kaks_comparison_summary.tsv",
        rows,
        [
            "comparison",
            "input_pairs",
            "used_alignments",
            "kaks_rows",
            "ks_median",
            "ks_mean",
            "ka_median",
            "kaks_ratio_median",
            "p_4dtv_median",
            "ks_gt_3_count",
            "p_4dtv_gt_0_75_count",
        ],
    )


def histogram(values: list[float], max_value: float, bin_width: float) -> list[int]:
    bins = int(max_value / bin_width)
    counts = [0] * bins
    for value in values:
        if value < 0 or value >= max_value:
            continue
        idx = min(int(value / bin_width), bins - 1)
        counts[idx] += 1
    return counts


def write_svg_histogram(path: Path, title: str, values_by_comp: dict[str, list[float]], max_value: float, bin_width: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    comps = sorted(values_by_comp)
    width = 1200
    row_h = 120
    left = 160
    top = 60
    plot_w = 960
    height = top + row_h * max(1, len(comps)) + 50
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;font-size:12px}.title{font-size:18px;font-weight:700}.axis{stroke:#333;stroke-width:1}.bar{fill:#4477AA;opacity:.82}</style>',
        f'<text class="title" x="{left}" y="32">{title}</text>',
    ]
    for row_idx, comp in enumerate(comps):
        y0 = top + row_idx * row_h
        counts = histogram(values_by_comp[comp], max_value, bin_width)
        ymax = max(counts) if counts else 1
        lines.append(f'<text x="{left - 10}" y="{y0 + 55}" text-anchor="end">{comp}</text>')
        lines.append(f'<line class="axis" x1="{left}" y1="{y0 + 90}" x2="{left + plot_w}" y2="{y0 + 90}"/>')
        for idx, count in enumerate(counts):
            if count == 0:
                continue
            x = left + plot_w * idx / len(counts)
            bw = max(1, plot_w / len(counts) - 1)
            bh = 80 * count / ymax
            lines.append(f'<rect class="bar" x="{x:.2f}" y="{y0 + 90 - bh:.2f}" width="{bw:.2f}" height="{bh:.2f}"><title>{comp} bin={idx} count={count}</title></rect>')
        lines.append(f'<text x="{left}" y="{y0 + 108}">0</text>')
        lines.append(f'<text x="{left + plot_w}" y="{y0 + 108}" text-anchor="end">{max_value:g}</text>')
    lines.append("</svg>")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_histograms(outdir: Path, merged_rows: list[dict[str, object]]) -> None:
    ks_by_comp: dict[str, list[float]] = {}
    tv_by_comp: dict[str, list[float]] = {}
    ratio_by_comp: dict[str, list[float]] = {}
    for row in merged_rows:
        comp = str(row["comparison"])
        ks = parse_float(str(row.get("Ks", "")))
        tv = parse_float(str(row.get("p_4dtv", "")))
        ratio = parse_float(str(row.get("Ka_Ks", "")))
        if ks is not None:
            ks_by_comp.setdefault(comp, []).append(ks)
        if tv is not None:
            tv_by_comp.setdefault(comp, []).append(tv)
        if ratio is not None and ratio <= 3:
            ratio_by_comp.setdefault(comp, []).append(ratio)
    write_svg_histogram(outdir / "figures" / "ks_distribution_by_comparison.svg", "Syntenic pair Ks distribution", ks_by_comp, 3.0, 0.05)
    write_svg_histogram(outdir / "figures" / "4dtv_distribution_by_comparison.svg", "Syntenic pair 4DTV distribution", tv_by_comp, 1.0, 0.025)
    write_svg_histogram(outdir / "figures" / "kaks_ratio_distribution_by_comparison.svg", "Syntenic pair Ka/Ks distribution", ratio_by_comp, 3.0, 0.05)


def main() -> int:
    args = parse_args()
    analysis_dir = Path(args.analysis_dir)
    jcvi_dir = Path(args.jcvi_dir)
    pair_table = Path(args.pair_table) if args.pair_table else jcvi_dir / "summary" / "jcvi_syntenic_pairs.tsv"
    outdir = Path(args.outdir)
    pair_dir = outdir / "pairs"
    axt_dir = outdir / "axt"
    kaks_dir = outdir / "kaks"
    tmp_root = outdir / "tmp_pair_alignments"
    for path in (pair_dir, axt_dir, kaks_dir, outdir / "summary", outdir / "figures", tmp_root):
        path.mkdir(parents=True, exist_ok=True)

    pairs, chromosome_filter_rows = load_pairs(
        pair_table,
        set(args.include_comparison),
        args.pair_mode,
        args.max_pairs_per_comparison,
        args.min_chrom_pair_anchors,
        args.min_chrom_pair_fraction,
        args.min_reference_block_anchors,
    )
    if not pairs:
        raise SystemExit("No syntenic gene pairs found for WGD/KaKs analysis.")
    labels = {p.query_label for p in pairs} | {p.ref_label for p in pairs}
    pep_by_label, cds_by_label = load_sequences(analysis_dir, jcvi_dir, labels)

    pairs_by_comp: dict[str, list[PairRecord]] = {}
    for pair in pairs:
        pairs_by_comp.setdefault(pair.comparison, []).append(pair)
    pair_counts = {comp: len(items) for comp, items in pairs_by_comp.items()}

    pair_plan_rows = [
        {
            "pair_id": p.pair_id,
            "comparison": p.comparison,
            "comparison_type": p.comparison_type,
            "query_label": p.query_label,
            "ref_label": p.ref_label,
            "block_id": p.block_id,
            "query_gene": p.query_gene,
            "ref_gene": p.ref_gene,
            "query_seqid": p.query_seqid,
            "ref_seqid": p.ref_seqid,
            "block_support": p.block_support,
            "chromosome_pair_support": p.chromosome_pair_support,
            "chromosome_filter_rule": p.chromosome_filter_rule,
        }
        for p in pairs
    ]
    write_tsv(
        outdir / "summary" / "syntenic_pair_plan.tsv",
        pair_plan_rows,
        [
            "pair_id", "comparison", "comparison_type", "query_label", "ref_label",
            "block_id", "query_gene", "ref_gene", "query_seqid", "ref_seqid",
            "block_support", "chromosome_pair_support", "chromosome_filter_rule",
        ],
    )
    write_tsv(
        outdir / "summary" / "chromosome_pair_filter.tsv",
        chromosome_filter_rows,
        [
            "comparison", "comparison_type", "query_seqid", "ref_seqid",
            "query_chromosome_key", "ref_chromosome_key", "anchor_support",
            "support_threshold", "filter_rule", "accepted",
        ],
    )

    all_results: list[AlignmentResult] = []
    all_merged_rows: list[dict[str, object]] = []
    for comparison in sorted(pairs_by_comp):
        comp_pairs = pairs_by_comp[comparison]
        comp_safe = safe_name(comparison)
        axt_path = axt_dir / f"{comp_safe}.axt"
        results: list[AlignmentResult] = []
        with ThreadPoolExecutor(max_workers=max(1, args.threads)) as executor:
            futures = [
                executor.submit(
                    process_pair,
                    pair,
                    pep_by_label,
                    cds_by_label,
                    args.mafft_cmd,
                    args.pal2nal_cmd,
                    args.min_aa,
                    args.min_codons,
                    args.min_4d_sites,
                    tmp_root,
                    bool(args.keep_temp),
                )
                for pair in comp_pairs
            ]
            for future in as_completed(futures):
                results.append(future.result())
        results.sort(key=lambda r: r.pair.pair_id)
        all_results.extend(results)
        with axt_path.open("w", encoding="utf-8") as axt_handle:
            for result in results:
                if result.status == "used":
                    axt_handle.write(result.axt_block)

        align_rows = [result_row(result) for result in results]
        write_tsv(
            pair_dir / f"{comp_safe}.alignment_stats.tsv",
            align_rows,
            [
                "pair_id",
                "comparison",
                "comparison_type",
                "query_label",
                "ref_label",
                "block_id",
                "query_gene",
                "ref_gene",
                "query_seqid",
                "ref_seqid",
                "block_support",
                "chromosome_pair_support",
                "chromosome_filter_rule",
                "status",
                "message",
                "aligned_codons",
                "ungapped_codons",
                "fourfold_sites",
                "fourfold_differences",
                "p_4dtv",
            ],
        )
        used = [r for r in results if r.status == "used"]
        kaks_table: dict[str, dict[str, str]] = {}
        kaks_out = kaks_dir / f"{comp_safe}.KaKs_Calculator.{args.kaks_method}.tsv"
        if used:
            code, err = run_kaks(args.kaks_cmd, args.kaks_method, axt_path, kaks_out)
            if code != 0:
                (kaks_dir / f"{comp_safe}.KaKs_Calculator.stderr.txt").write_text(err + "\n", encoding="utf-8")
                raise SystemExit(f"KaKs_Calculator failed for {comparison}; see {comp_safe}.KaKs_Calculator.stderr.txt")
            kaks_table = parse_kaks_table(kaks_out)
            if len(kaks_table) != len(used):
                expected = {result.pair.pair_id for result in used}
                missing = sorted(expected - set(kaks_table))
                missing_path = kaks_dir / f"{comp_safe}.missing_kaks_pair_ids.txt"
                missing_path.write_text("\n".join(missing) + "\n", encoding="utf-8")
                raise SystemExit(
                    f"KaKs_Calculator returned {len(kaks_table)}/{len(used)} rows for {comparison}; "
                    f"missing IDs: {missing_path}"
                )
        merged = merge_kaks(results, kaks_table, kaks_dir / f"{comp_safe}.merged_4dtv_kaks.tsv")
        all_merged_rows.extend(merged)

    all_align_rows = [result_row(result) for result in all_results]
    write_tsv(
        outdir / "summary" / "all_pair_alignment_stats.tsv",
        all_align_rows,
        [
            "pair_id",
            "comparison",
            "comparison_type",
            "query_label",
            "ref_label",
            "block_id",
            "query_gene",
            "ref_gene",
            "query_seqid",
            "ref_seqid",
            "block_support",
            "chromosome_pair_support",
            "chromosome_filter_rule",
            "status",
            "message",
            "aligned_codons",
            "ungapped_codons",
            "fourfold_sites",
            "fourfold_differences",
            "p_4dtv",
        ],
    )
    write_tsv(
        outdir / "summary" / "all_merged_4dtv_kaks.tsv",
        all_merged_rows,
        [
            "pair_id",
            "comparison",
            "comparison_type",
            "query_label",
            "ref_label",
            "block_id",
            "query_gene",
            "ref_gene",
            "query_seqid",
            "ref_seqid",
            "block_support",
            "chromosome_pair_support",
            "chromosome_filter_rule",
            "aligned_codons",
            "ungapped_codons",
            "fourfold_sites",
            "fourfold_differences",
            "p_4dtv",
            "kaks_status",
            "method",
            "Ka",
            "Ks",
            "Ka_Ks",
            "P_value",
        ],
    )
    write_summary(outdir, pair_counts, all_align_rows, all_merged_rows)
    write_histograms(outdir, all_merged_rows)

    if not args.keep_temp:
        shutil.rmtree(tmp_root, ignore_errors=True)
    print(f"comparisons={len(pairs_by_comp)} input_pairs={len(pairs)} aligned_pairs={sum(1 for r in all_results if r.status == 'used')} kaks_rows={len(all_merged_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
