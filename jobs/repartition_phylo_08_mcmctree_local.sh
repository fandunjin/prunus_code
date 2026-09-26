#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Dating wrappers using the existing MCMCtree parameterization.
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/analysisphylogenetics}"
TRIMMED_DIR="${TRIMMED_DIR:-${BASE}/04_alignments/main/trimmed}"
PARTITION_FILE="${PARTITION_FILE:-${BASE}/04_alignments/main/qc/concat_partitions.nex}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
IQTREE_BIN="${IQTREE_BIN:-iqtree3}"
MCMCTREE_BIN="${MCMCTREE_BIN:-mcmctree}"
THREADS="${THREADS:-8}"

[[ -e "${BASE}/08_final_summary/COMPLETE" ]] || { echo "[ERROR] Complete ASTRAL/support workflow first" >&2; exit 2; }
[[ -d "${TRIMMED_DIR}" && -s "${PARTITION_FILE}" ]] || { echo "[ERROR] Missing new trimmed alignments or partitions" >&2; exit 2; }
command -v "${IQTREE_BIN}" >/dev/null 2>&1 || { echo "[ERROR] IQ-TREE unavailable: ${IQTREE_BIN}" >&2; exit 2; }
command -v "${MCMCTREE_BIN}" >/dev/null 2>&1 || { echo "[ERROR] MCMCtree unavailable: ${MCMCTREE_BIN}" >&2; exit 2; }

OUT14="${BASE}/09_mcmctree_species14_plant_outgroup_fossil55_35"
OUT17="${BASE}/09_mcmctree_subgenome17_outgroup_fossil55_35"
ROOT_DIR="${ROOT_DIR}" OUT_ROOT="${OUT14}" TRIMMED_DIR="${TRIMMED_DIR}" PARTITION_FILE="${PARTITION_FILE}" \
  PYTHON_BIN="${PYTHON_BIN}" IQTREE_BIN="${IQTREE_BIN}" MCMCTREE_BIN="${MCMCTREE_BIN}" THREADS="${THREADS}" \
  bash "${ROOT_DIR}/jobs/figure2_mcmctree_species14_plant_outgroup_fossil55_35_local.sh"
ROOT_DIR="${ROOT_DIR}" OUT_ROOT="${OUT17}" TRIMMED_DIR="${TRIMMED_DIR}" PARTITION_FILE="${PARTITION_FILE}" \
  PYTHON_BIN="${PYTHON_BIN}" IQTREE_BIN="${IQTREE_BIN}" MCMCTREE_BIN="${MCMCTREE_BIN}" THREADS="${THREADS}" \
  bash "${ROOT_DIR}/jobs/figure2_mcmctree_subgenome17_outgroup_fossil55_35_local.sh"
printf 'MCMCtree dating runs completed from new repartitioned alignments.\n' > "${BASE}/MCMCtree_COMPLETE"
echo "[OK] New repartitioned MCMCtree species14 and subgenome17 runs completed"
