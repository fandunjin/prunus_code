#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
PRIMARY_DIR="${PRIMARY_DIR:-${ROOT_DIR}/annotation/genes/primary_isoform_busco/primary}"
SCRIPT="${ROOT_DIR}/scripts/package_pzi_downstream_inputs.py"
CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
PY_ENV="${PY_ENV:-haphic}"

require_file() {
    local file_path="$1"
    local label="$2"
    if [ ! -s "${file_path}" ]; then
        echo "ERROR: ${label} does not exist or is empty: ${file_path}" >&2
        exit 1
    fi
}

require_file "${CONDA_SH}" "Conda initialization script"
require_file "${SCRIPT}" "downstream input packaging script"
require_file "${PRIMARY_DIR}/primary_summary.tsv" "primary isoform summary"

mkdir -p "${ROOT_DIR}/logs"

source "${CONDA_SH}"
conda activate "${PY_ENV}"

if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi

echo "=== Package PZI downstream analysis inputs ==="
echo "Root       : ${ROOT_DIR}"
echo "Primary dir: ${PRIMARY_DIR}"
echo "Output     : ${OUT_DIR}"
echo

"${PYTHON_BIN}" "${SCRIPT}" \
    --root-dir "${ROOT_DIR}" \
    --primary-dir "${PRIMARY_DIR}" \
    --out-dir "${OUT_DIR}"

echo "=== Downstream input package created ==="
echo "${OUT_DIR}"
echo "${OUT_DIR}/README_DOWNSTREAM_INPUTS.md"
echo "${OUT_DIR}/SUMMARY.tsv"
echo "${OUT_DIR}/MANIFEST.tsv"
