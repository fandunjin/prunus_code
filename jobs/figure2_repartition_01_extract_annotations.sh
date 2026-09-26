#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/scripts/extract_repartitioned_haplotype_annotations.py"
require_file "${ROOT_DIR}/haploid_repartition/maps/pzi_repartition_map.tsv"
require_file "${ROOT_DIR}/haploid_repartition/maps/psp_repartition_map.tsv"
require_file "${ROOT_DIR}/haploid_repartition/maps/ppa_repartition_map.tsv"
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || { echo "[ERROR] Python is unavailable: ${PYTHON_BIN}" >&2; exit 2; }

args=("${ROOT_DIR}/scripts/extract_repartitioned_haplotype_annotations.py" --project-root "${ROOT_DIR}")
[[ "${FORCE}" == "1" ]] && args+=(--force)
"${PYTHON_BIN}" "${args[@]}"

echo "[OK] Revised annotation extraction completed."
