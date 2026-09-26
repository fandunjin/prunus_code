#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
SCRIPT="${ROOT_DIR}/scripts/collect_final_haplotype_genomes.py"
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

require_file "${SCRIPT}" "final haplotype genome collection script"
require_file "${CONDA_SH}" "Conda initialization script"
mkdir -p "${OUT_DIR}" "${ROOT_DIR}/logs"

source "${CONDA_SH}"
conda activate "${PY_ENV}"

if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi

echo "=== Collect final phased haplotype genomes ==="
echo "Root: ${ROOT_DIR}"
echo "Out : ${OUT_DIR}"
echo

"${PYTHON_BIN}" "${SCRIPT}" \
    --root-dir "${ROOT_DIR}" \
    --out-dir "${OUT_DIR}"

echo "=== Final haplotype genome package ==="
echo "${OUT_DIR}/pzi_hapA1.fa"
echo "${OUT_DIR}/pzi_hapA1.gff3"
echo "${OUT_DIR}/pzi_hapA1.cds.fa"
echo "${OUT_DIR}/pzi_hapA1.pep.fa"
echo "${OUT_DIR}/pzi_hapA2.fa"
echo "${OUT_DIR}/pzi_hapB1.fa"
echo "${OUT_DIR}/pzi_hapB2.fa"
echo "${OUT_DIR}/MANIFEST.tsv"
echo "${OUT_DIR}/SUMMARY.tsv"
echo "${OUT_DIR}/CHECKSUMS.md5"
echo "${OUT_DIR}/README_FINAL_GENOMES.md"
