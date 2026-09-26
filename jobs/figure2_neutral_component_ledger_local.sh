#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Figure 2 stage 6: construct neutral components only after manual confirmation.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BLOCK_ROOT="${BLOCK_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/neutral_block_ancestry}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
REVIEW_SHEET="${REVIEW_SHEET:-${BLOCK_ROOT}/05_junction_manual_review.tsv}"
COPY1_FAI="${COPY1_FAI:-${ROOT_DIR}/final_genome/pzi_copy1.fa.fai}"
COPY_MAP="${COPY_MAP:-${ROOT_DIR}/final_assembly/haplotype_provenance/pzi_2n32.corrected_homolog_pairs.tsv}"
OUT_ROOT="${OUT_ROOT:-${BLOCK_ROOT}/06_neutral_component_ledger}"

if [[ ! -s "${COPY1_FAI}" ]]; then
    COPY1_FASTA="${COPY1_FAI%.fai}"
    command -v samtools >/dev/null 2>&1 || { echo "[ERROR] samtools is required to index ${COPY1_FASTA}" >&2; exit 2; }
    [[ -s "${COPY1_FASTA}" ]] || { echo "[ERROR] Missing copy1 FASTA: ${COPY1_FASTA}" >&2; exit 2; }
    samtools faidx "${COPY1_FASTA}"
fi

for f in "${REVIEW_SHEET}" "${BLOCK_ROOT}/03_syntenic_fragments/syntenic_fragments.tsv" "${COPY1_FAI}" "${COPY_MAP}"; do
    [[ -s "${f}" ]] || { echo "[ERROR] Missing input: ${f}" >&2; exit 2; }
done
if [[ -e "${OUT_ROOT}/COMPLETE" ]]; then
    echo "[OK] Existing neutral component ledger: ${OUT_ROOT}"
    exit 0
fi
mkdir -p "${OUT_ROOT}"
"${PYTHON_BIN}" "${ROOT_DIR}/scripts/build_neutral_component_ledger.py" \
    --review-sheet "${REVIEW_SHEET}" \
    --fragments "${BLOCK_ROOT}/03_syntenic_fragments/syntenic_fragments.tsv" \
    --query-fai "${COPY1_FAI}" \
    --copy-map "${COPY_MAP}" \
    --outdir "${OUT_ROOT}"
touch "${OUT_ROOT}/COMPLETE"
echo "[OK] Neutral component ledger completed: ${OUT_ROOT}"
