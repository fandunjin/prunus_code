#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/haploid_repartition/inputs}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/haploid_repartition/figure}"
WORK_ROOT="${WORK_ROOT:-${OUT_ROOT}/work}"
SCRIPT="${ROOT_DIR}/scripts/plot_repartitioned_pairwise_heatmap.py"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
MINIMAP_ENV="${MINIMAP_ENV:-minimap}"
THREADS="${THREADS:-32}"
MINIMAP2_PRESET="${MINIMAP2_PRESET:-asm10}"

die() { echo "[ERROR] $*" >&2; exit 2; }
req() { [[ -s "$1" ]] || die "Missing or empty $2: $1"; }
req "${SCRIPT}" "heatmap plotting script"
req "${CONDA_SH}" "conda initialization"
mkdir -p "${OUT_ROOT}" "${WORK_ROOT}" "${ROOT_DIR}/logs"

source "${CONDA_SH}"
conda activate "${MINIMAP_ENV}"
command -v minimap2 >/dev/null 2>&1 || die "minimap2 is unavailable in ${MINIMAP_ENV}"
command -v python >/dev/null 2>&1 || die "python is unavailable in ${MINIMAP_ENV}"

cat "${INPUT_ROOT}/pzi/genome/pzi_hapA1.fa" "${INPUT_ROOT}/pzi/genome/pzi_hapB1.fa" > "${WORK_ROOT}/pzi_A1B1.fa"
cat "${INPUT_ROOT}/psp/genome/psp_hapA1.fa" "${INPUT_ROOT}/psp/genome/psp_hapB1.fa" > "${WORK_ROOT}/psp_A1B1.fa"
cat "${INPUT_ROOT}/ppa/genome/ppa_hapA.fa" "${INPUT_ROOT}/ppa/genome/ppa_hapB.fa" > "${WORK_ROOT}/ppa_AB.fa"

run_pair() {
    local tag="$1" target="$2" query="$3"
    local paf="${WORK_ROOT}/${tag}.${MINIMAP2_PRESET}.paf"
    if [[ ! -s "${paf}" ]]; then
        minimap2 -x "${MINIMAP2_PRESET}" -c --secondary=yes -N 50 -p 0.50 \
            -t "${THREADS}" "${target}" "${query}" > "${paf}"
    else
        echo "[SKIP] Existing PAF: ${paf}"
    fi
    python "${SCRIPT}" \
        --paf "${paf}" \
        --query-fasta "${query}" \
        --target-fasta "${target}" \
        --out-stem "${OUT_ROOT}/${tag}.${MINIMAP2_PRESET}.heatmap" \
        --title "${tag} (${MINIMAP2_PRESET})" \
        --min-mapq 0 --min-aln-len 10000
}

run_pair pzi_vs_psp "${WORK_ROOT}/pzi_A1B1.fa" "${WORK_ROOT}/psp_A1B1.fa"
run_pair pzi_vs_ppa "${WORK_ROOT}/pzi_A1B1.fa" "${WORK_ROOT}/ppa_AB.fa"
run_pair psp_vs_ppa "${WORK_ROOT}/psp_A1B1.fa" "${WORK_ROOT}/ppa_AB.fa"

echo "[OK] Pairwise heatmaps: ${OUT_ROOT}"
