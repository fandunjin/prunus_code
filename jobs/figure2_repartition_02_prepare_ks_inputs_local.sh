#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Standardize new haplotype annotations and canonical references for five Ks panels.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/downstream_analysis/figure2_repartition_ks}"
REPARTITION_ROOT="${REPARTITION_ROOT:-${ROOT_DIR}/haploid_repartition}"
CANONICAL_ROOT="${CANONICAL_ROOT:-${ROOT_DIR}/downstream_analysis/ks/stage/canonical}"
GTR_ROOT="${GTR_ROOT:-${ROOT_DIR}/data/gtr}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
FORCE="${FORCE:-0}"

die() { echo "[ERROR] $*" >&2; exit 2; }
[[ -s "${REPARTITION_ROOT}/inputs/ANNOTATIONS_COMPLETE" ]] || die "Missing corrected annotation completion marker: ${REPARTITION_ROOT}/inputs/ANNOTATIONS_COMPLETE"
[[ -s "${ROOT_DIR}/scripts/prepare_figure2_repartition_ks_inputs.py" ]] || die "Missing input preparation script"
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || die "Python unavailable: ${PYTHON_BIN}"

if [[ -e "${OUT_DIR}/COMPLETE" && "${FORCE}" != "1" ]]; then
  echo "[OK] Reusing corrected Figure 2 Ks inputs: ${OUT_DIR}"
  cat "${OUT_DIR}/input_validation.tsv"
  exit 0
fi

args=("${ROOT_DIR}/scripts/prepare_figure2_repartition_ks_inputs.py" \
  --project-root "${ROOT_DIR}" \
  --repartition-root "${REPARTITION_ROOT}" \
  --canonical-root "${CANONICAL_ROOT}" \
  --gtr-root "${GTR_ROOT}" \
  --outdir "${OUT_DIR}")
[[ "${FORCE}" == "1" ]] && args+=(--force)
"${PYTHON_BIN}" "${args[@]}"

[[ -s "${OUT_DIR}/COMPLETE" ]] || die "Input preparation did not create COMPLETE"
echo "[OK] Corrected Figure 2 Ks inputs: ${OUT_DIR}"
