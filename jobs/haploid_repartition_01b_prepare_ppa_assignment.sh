#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/haploid_repartition/ppa_copy1_minimap/ppa_vs_pzi_copy1.asm5.best_pzi_copy1_hit_per_ppa.tsv}"
OUT_ASSIGNMENT="${OUT_ASSIGNMENT:-${ROOT_DIR}/haploid_repartition/assignments/ppa_A_B.rename_map.tsv}"
SCRIPT="${ROOT_DIR}/scripts/prepare_repartitioned_ppa_assignment.py"
REVIEW_STATUS="${REVIEW_STATUS:-approved}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
PY_ENV="${PY_ENV:-haphic}"

die() { echo "[ERROR] $*" >&2; exit 2; }
[[ -s "${BEST_HITS}" ]] || die "Missing PPA-PZI copy1 best-hit table: ${BEST_HITS}"
[[ -s "${SCRIPT}" ]] || die "Missing PPA assignment script: ${SCRIPT}"
[[ -s "${CONDA_SH}" ]] || die "Missing conda initialization: ${CONDA_SH}"

source "${CONDA_SH}"
conda activate "${PY_ENV}"
command -v python >/dev/null 2>&1 || die "python is unavailable in ${PY_ENV}"

python "${SCRIPT}" \
    --best-hits "${BEST_HITS}" \
    --out "${OUT_ASSIGNMENT}" \
    --review-status "${REVIEW_STATUS}"

echo "[OK] Revised PPA assignment: ${OUT_ASSIGNMENT}"
