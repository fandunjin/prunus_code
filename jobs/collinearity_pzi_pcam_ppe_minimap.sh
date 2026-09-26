#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

# FASTA-level (annotation-free) macro-synteny between PZI copy1 and the
# corymbose P. campanulata (pcam) and solitary-flowered P. persica (ppe)
# genomes under ${ROOT_DIR}/data/. Same visualisation and filters as the
# pun/oce job; minimap2 reads .fna(.gz) directly. Pairs use "query:target".

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/analyse/collinearity_pzi_pcam_ppe}"
PZI_FASTA="${PZI_FASTA:-${PROJECT_ROOT}/final_genome/pzi_copy1.fa}"
PCAM_FASTA="${PCAM_FASTA:-${ROOT_DIR}/data/pcam/GCA_050084535.1_PCAv2.0_genomic.fna.gz}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
MINIMAP_ENV="${MINIMAP_ENV:-minimap}"
REPORT_SCRIPT="${REPORT_SCRIPT:-${ROOT_DIR}/jobs/minimap_collinearity_report.py}"

THREADS="${THREADS:-36}"
MINIMAP_PRESET="${MINIMAP_PRESET:-asm10}"
MIN_MAPQ="${MIN_MAPQ:-20}"
MIN_ALN_LEN="${MIN_ALN_LEN:-10000}"
PAIRS="${PAIRS:-pcam:pzi ppe:pzi}"
FORCE_RERUN="${FORCE_RERUN:-0}"

first_existing() {
    local p
    for p in "$@"; do
        if [[ -s "${p}" ]]; then echo "${p}"; return 0; fi
    done
    return 1
}

PPE_FASTA="${PPE_FASTA:-$(first_existing \
    "${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fna" \
    "${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fna.gz" || true)}"
GTR_FASTA="${GTR_FASTA:-$(first_existing \
    "${ROOT_DIR}/data/gtr/GCA_018257905.1_V1_genomic.fna" \
    "${ROOT_DIR}/data/gtr/GCA_018257905.1_V1_genomic.fna.gz" || true)}"

require_file() {
    local path="$1"; local label="$2"
    [[ -s "${path}" ]] || { echo "[ERROR] Missing ${label}: ${path}" >&2; exit 2; }
}

require_file "${PZI_FASTA}" "PZI copy1 FASTA"
require_file "${PCAM_FASTA}" "PCAM FASTA"
require_file "${CONDA_SH}" "conda initialization script"
require_file "${REPORT_SCRIPT}" "collinearity report script"

mkdir -p "${WORK_DIR}/00_logs" "${WORK_DIR}/02_paf" "${WORK_DIR}/03_report" "${WORK_DIR}/04_meta"
mkdir -p "${ROOT_DIR}/logs"

source "${CONDA_SH}"
conda activate "${MINIMAP_ENV}"
command -v minimap2 >/dev/null 2>&1 || { echo "[ERROR] minimap2 is unavailable in ${MINIMAP_ENV}" >&2; exit 2; }
if command -v python3 >/dev/null 2>&1; then PYTHON_BIN=python3; else PYTHON_BIN=python; fi
"${PYTHON_BIN}" -c 'import matplotlib, numpy' 2>/dev/null || {
    echo "[ERROR] ${MINIMAP_ENV} python lacks matplotlib/numpy" >&2; exit 2; }

declare -A FA=( [pzi]="${PZI_FASTA}" [pcam]="${PCAM_FASTA}" [ppe]="${PPE_FASTA:-}" [gtr]="${GTR_FASTA:-}" )
for pair in ${PAIRS}; do
    for sp in "${pair%%:*}" "${pair##*:}"; do
        [[ -n "${FA[$sp]:-}" ]] || { echo "[ERROR] No FASTA configured for species '${sp}' (pairs: ${PAIRS})" >&2; exit 2; }
        require_file "${FA[$sp]}" "${sp} FASTA"
    done
done

NP=$(wc -w <<<"${PAIRS}")
TP=$(( THREADS / NP )); (( TP < 4 )) && TP=4

align_pair() {
    local pair="$1" q="${1%%:*}" s="${1##*:}"
    local raw="${WORK_DIR}/02_paf/${q}_vs_${s}.${MINIMAP_PRESET}.paf"
    local log="${WORK_DIR}/00_logs/${q}_vs_${s}.minimap2.log"
    if [[ "${FORCE_RERUN}" != "1" && -s "${raw}" ]]; then
        echo "[INFO] Reusing existing PAF: ${raw}"
        return 0
    fi
    rm -f "${raw}"
    echo "[RUN] minimap2 ${q} (query) vs ${s} (target), preset=${MINIMAP_PRESET}, threads=${TP}"
    minimap2 -x "${MINIMAP_PRESET}" -t "${TP}" "${FA[$s]}" "${FA[$q]}" > "${raw}" 2> "${log}"
    local n; n=$(wc -l < "${raw}")
    echo "[OK] ${q}_vs_${s}: ${n} raw PAF lines"
}

report_pair() {
    local pair="$1" q="${1%%:*}" s="${1##*:}"
    local raw="${WORK_DIR}/02_paf/${q}_vs_${s}.${MINIMAP_PRESET}.paf"
    local out="${WORK_DIR}/03_report/${q}_vs_${s}"
    local log="${WORK_DIR}/00_logs/${q}_vs_${s}.report.log"
    require_file "${raw}" "raw PAF for ${pair}"
    "${PYTHON_BIN}" "${REPORT_SCRIPT}" \
        --paf "${raw}" \
        --outdir "${out}" \
        --label "${q}_vs_${s}" \
        --title "${q} vs ${s} minimap2 macro-synteny (${MINIMAP_PRESET})" \
        --min-mapq "${MIN_MAPQ}" \
        --min-alignment-len "${MIN_ALN_LEN}" \
        > "${log}" 2>&1
    echo "[OK] report ${q}_vs_${s}: ${out}/dotplot.svg, ${out}/blocks.tsv, ${out}/summary.tsv"
}

echo "=== PZI copy1 vs PCAM / PPE macro-synteny ==="
echo "root    : ${ROOT_DIR}"
echo "output  : ${WORK_DIR}"
echo "pairs   : ${PAIRS}"
echo "preset  : ${MINIMAP_PRESET}; MAPQ >= ${MIN_MAPQ}; length >= ${MIN_ALN_LEN} bp"

pids=()
for pair in ${PAIRS}; do
    align_pair "${pair}" &
    pids+=("$!")
done
for pid in "${pids[@]}"; do
    wait "${pid}" || { echo "[ERROR] one alignment job failed" >&2; exit 3; }
done

for pair in ${PAIRS}; do
    report_pair "${pair}"
done

{
    printf 'pair\tn_anchors\tcovered_target_Mb\tpct_target_covered\tcovered_query_Mb\tpct_query_covered\tn_blocks\tn_blocks_ge5\tanchor_in_block_frac\tmax_block_Mb\n'
    for pair in ${PAIRS}; do
        q="${pair%%:*}"; s="${pair##*:}"
        tail -n +2 "${WORK_DIR}/03_report/${q}_vs_${s}/summary.tsv"
    done
} > "${WORK_DIR}/04_meta/ALL_PAIRS_summary.tsv"

{
    echo "workflow: minimap2 macro-synteny (PAF) for pzi_copy1 / pcam / ppe"
    echo "date: $(date -Is)"
    echo "host: $(hostname)"
    echo "work_dir: ${WORK_DIR}"
    echo "pairs: ${PAIRS}"
    echo "preset: ${MINIMAP_PRESET}; MAPQ >= ${MIN_MAPQ}; min length: ${MIN_ALN_LEN}; threads: ${THREADS}"
    echo "minimap2: $(minimap2 --version)"
    "${PYTHON_BIN}" -c 'import sys, numpy, matplotlib; print("python", sys.version.split()[0], "| numpy", numpy.__version__, "| matplotlib", matplotlib.__version__)'
    for sp in pzi pcam ppe gtr; do
        [[ -n "${FA[$sp]:-}" ]] || continue
        f="${FA[$sp]}"; printf 'input %s: %s (%s bytes)\n' "${sp}" "${f}" "$(stat -c '%s' "${f}" 2>/dev/null || echo NA)"
    done
} > "${WORK_DIR}/04_meta/tool_versions_and_inputs.txt" 2>&1

cat > "${WORK_DIR}/README.txt" <<EOF
PZI copy1 vs PCAM / PPE FASTA-level macro-synteny
=================================================

Inputs
------
- pzi : ${PZI_FASTA}
- pcam: ${PCAM_FASTA}
- ppe : ${PPE_FASTA:-<not configured>}
- pairs (query:target): ${PAIRS}
- optional: gtr=${GTR_FASTA:-<not configured>} can be added via
  PAIRS="pcam:pzi ppe:pzi gtr:pzi" (all pairwise combinations are allowed).

Method and filters
------------------
- minimap2 -x ${MINIMAP_PRESET} per pair (intrageneric Prunus comparisons;
  use MINIMAP_PRESET=asm5 for speed or asm20 for maximum sensitivity).
- Retained anchors: MAPQ >= ${MIN_MAPQ}, alignment length >= ${MIN_ALN_LEN} bp.
- Sequence lengths are read from the PAF (columns 2/7); no .fai is required.

Outputs
-------
- 02_paf/<query>_vs_<target>.<preset>.paf   raw minimap2 alignments
- 03_report/<query>_vs_<target>/filtered.paf
- 03_report/<query>_vs_<target>/dotplot.svg|.png   black-point anchor dotplot
- 03_report/<query>_vs_<target>/blocks.tsv   collinear blocks (>= 3 anchors)
- 03_report/<query>_vs_<target>/summary.tsv  one-line pair statistics
- 03_report/<query>_vs_<target>/stats.txt
- 04_meta/ALL_PAIRS_summary.tsv, tool_versions_and_inputs.txt

Re-runs
-------
- Default: existing PAF files are reused (FORCE_RERUN=0).
- Force realign: FORCE_RERUN=1.

Boundary
--------
FASTA-level structural correspondence only; no orthology or parental-origin
inference. Gene-level synteny, Ks and phylogenetics remain separate steps.
EOF

echo "[DONE] ${WORK_DIR}"
echo "[OK] combined summary: ${WORK_DIR}/04_meta/ALL_PAIRS_summary.tsv"
