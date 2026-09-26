#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Figure 2 stage 5: parameterized HiFi support review for stage-4 candidates.
# The caller must explicitly supply HIFI_BAM because the project contains more
# than one read technology and not every alignment uses copy1 coordinates.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BLOCK_ROOT="${BLOCK_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/neutral_block_ancestry}"
OUT_ROOT="${OUT_ROOT:-${BLOCK_ROOT}/05_hifi_junction_review}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
HIFI_BAM="${HIFI_BAM:-}"
HIFI_BAI="${HIFI_BAI:-}"
MIN_MAPQ="${MIN_MAPQ:-20}"
MIN_FLANK_BP="${MIN_FLANK_BP:-1000}"
MIN_ALIGNMENT_BP="${MIN_ALIGNMENT_BP:-5000}"
MIN_SPANNING_READS="${MIN_SPANNING_READS:-3}"

die() { echo "[ERROR] $*" >&2; exit 2; }
req() { [[ -s "$1" ]] || die "Missing input: $1"; }

req "${BLOCK_ROOT}/04_candidate_validation/junction_100kb_review_windows.bed"
[[ -n "${HIFI_BAM}" ]] || die "Set HIFI_BAM to a coordinate-sorted HiFi BAM aligned to the exact copy1 reference."
req "${HIFI_BAM}"
if [[ -n "${HIFI_BAI}" ]]; then
    req "${HIFI_BAI}"
elif [[ ! -s "${HIFI_BAM}.bai" ]]; then
    die "Set HIFI_BAI because ${HIFI_BAM}.bai is unavailable."
fi
"${PYTHON_BIN}" -c 'import pysam' >/dev/null 2>&1 || die "pysam is unavailable in ${PYTHON_BIN}"

if [[ -e "${OUT_ROOT}/COMPLETE" ]]; then
    echo "[OK] Existing HiFi junction review: ${OUT_ROOT}"
    exit 0
fi
mkdir -p "${OUT_ROOT}"

args=(
    --windows "${BLOCK_ROOT}/04_candidate_validation/junction_100kb_review_windows.bed"
    --bam "${HIFI_BAM}"
    --outdir "${OUT_ROOT}"
    --min-mapq "${MIN_MAPQ}"
    --min-flank-bp "${MIN_FLANK_BP}"
    --min-alignment-bp "${MIN_ALIGNMENT_BP}"
    --min-spanning-reads "${MIN_SPANNING_READS}"
)
[[ -z "${HIFI_BAI}" ]] || args+=(--bai "${HIFI_BAI}")
"${PYTHON_BIN}" "${ROOT_DIR}/scripts/summarize_junction_hifi_support.py" "${args[@]}"

{
    echo "python=$(${PYTHON_BIN} --version 2>&1)"
    echo "hifi_bam=${HIFI_BAM}"
    echo "hifi_bai=${HIFI_BAI:-${HIFI_BAM}.bai}"
    echo "min_mapq=${MIN_MAPQ}"
    echo "min_flank_bp=${MIN_FLANK_BP}"
    echo "min_alignment_bp=${MIN_ALIGNMENT_BP}"
    echo "min_spanning_reads=${MIN_SPANNING_READS}"
    echo "claim_boundary=HiFi support is a review criterion; confirm coordinate compatibility and inspect Hi-C before calling a junction."
} > "${OUT_ROOT}/software_and_parameters.txt"
touch "${OUT_ROOT}/COMPLETE"
echo "[OK] HiFi junction support review completed: ${OUT_ROOT}"
