#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
PREFIX="pzi_2n32"
SCRIPT="${ROOT_DIR}/scripts/package_downstream_genome.py"
CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
PY_ENV="haphic"

require_file() {
    local file_path="$1"
    local label="$2"
    if [ ! -s "${file_path}" ]; then
        echo "ERROR: ${label} does not exist or is empty: ${file_path}" >&2
        exit 1
    fi
}

require_file "${SCRIPT}" "packaging Python script"
require_file "${CONDA_SH}" "Conda initialization script"
mkdir -p "${ROOT_DIR}/logs"

source "${CONDA_SH}"
conda activate "${PY_ENV}"

if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi

echo "=== Build downstream package ==="
echo "Source : ${SOURCE_DIR}"
echo "Output : ${FINAL_DIR}"
echo "Prefix : ${PREFIX}"
echo

"${PYTHON_BIN}" "${SCRIPT}" \
    --source-dir "${SOURCE_DIR}" \
    --final-dir "${FINAL_DIR}" \
    --prefix "${PREFIX}"

echo "=== Package created ==="
echo "${FINAL_DIR}"
echo "${FINAL_DIR}/MANIFEST.tsv"
echo "${FINAL_DIR}/CHECKSUMS.md5"
echo "${FINAL_DIR}/README_FINAL_PACKAGE.md"
