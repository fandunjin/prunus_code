#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
OUT_DIR="${ROOT_DIR}/annotation/chr32_lai"
PREFIX="pzi_chr32"
THREADS=24
CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
LAI_ENV="LAI"
PREP_SCRIPT="${ROOT_DIR}/scripts/prepare_lai_inputs.py"
GFF2BED_SCRIPT="${GFF2BED_SCRIPT:-}"

require_file() {
    local file_path="$1"
    local label="$2"
    if [ ! -s "${file_path}" ]; then
        echo "ERROR: ${label} does not exist or is empty: ${file_path}" >&2
        exit 1
    fi
}

require_file "${PREP_SCRIPT}" "LAI preparation Python script"
require_file "${CONDA_SH}" "Conda initialization script"
mkdir -p "${OUT_DIR}" "${ROOT_DIR}/logs"

source "${CONDA_SH}"
conda activate "${LAI_ENV}"

if [ -z "${GFF2BED_SCRIPT}" ]; then
    if [ -s "${CONDA_ROOT}/envs/EDTA/share/EDTA/bin/gff2bed.pl" ]; then
        GFF2BED_SCRIPT="${CONDA_ROOT}/envs/EDTA/share/EDTA/bin/gff2bed.pl"
    elif [ -s "${CONDA_ROOT}/envs/EDTA/share/LTR_retriever/bin/gff2bed.pl" ]; then
        GFF2BED_SCRIPT="${CONDA_ROOT}/envs/EDTA/share/LTR_retriever/bin/gff2bed.pl"
    elif command -v gff2bed.pl >/dev/null 2>&1; then
        GFF2BED_SCRIPT="gff2bed.pl"
    else
        GFF2BED_SCRIPT="gff2bed.pl"
        echo "WARNING: gff2bed.pl was not found in PATH. This is OK only if native LAI inputs already exist." >&2
        echo "If EDTA-GFF3 fallback is needed, set GFF2BED_SCRIPT=/path/to/gff2bed.pl before submitting." >&2
    fi
fi

if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi

"${PYTHON_BIN}" "${PREP_SCRIPT}" --edta-dir "${EDTA_DIR}" --out-dir "${OUT_DIR}" --prefix "${PREFIX}" --gff2bed "${GFF2BED_SCRIPT}" --method auto

source "${OUT_DIR}/${PREFIX}.lai_inputs.env"

require_file "${GENOME}" "LAI genome"
require_file "${INTACT}" "LAI intact LTR pass.list"
require_file "${ALL}" "LAI all LTR RepeatMasker-like output"

cd "${OUT_DIR}"

echo "=== LAI: 2n=32 chromosomes ==="
echo "Genome : ${GENOME}"
echo "Intact : ${INTACT}"
echo "All    : ${ALL}"
echo "Result : ${LAI_RESULT}"
echo "Method : ${LAI_INPUT_METHOD}"
echo

if [ ! -s "${LAI_RESULT}" ]; then
    LAI -genome "${GENOME}" -intact "${INTACT}" -all "${ALL}" -q -t "${THREADS}"
else
    echo "Existing LAI output found, skipping LAI: ${LAI_RESULT}"
fi

require_file "${LAI_RESULT}" "LAI output"
echo "LAI chr32 completed: ${LAI_RESULT}"
