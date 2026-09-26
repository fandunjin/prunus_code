#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
GFF3="${ROOT_DIR}/annotation/genes/final_${HAP_SET}_helixer/pzi_${HAP_SET}.helixer.gff3"
TE_GFF3="${TE_GFF3:-}"
OUT_DIR="${ROOT_DIR}/annotation/genes/final_${HAP_SET}_helixer/postprocess"
PREFIX="pzi_${HAP_SET}.helixer"
SCRIPT="${ROOT_DIR}/scripts/helixer_postprocess_qc.py"
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

case "${HAP_SET}" in
    hapA1|hapA2|hapB1|hapB2) ;;
    *)
        echo "ERROR: HAP_SET must be hapA1, hapA2, hapB1, or hapB2. Observed: ${HAP_SET}" >&2
        exit 1
        ;;
esac

require_file "${GENOME}" "final phased hap FASTA"
require_file "${GFF3}" "Helixer GFF3"
require_file "${SCRIPT}" "Helixer postprocess/QC Python script"
require_file "${CONDA_SH}" "Conda initialization script"
mkdir -p "${OUT_DIR}" "${ROOT_DIR}/logs"

source "${CONDA_SH}"
conda activate "${PY_ENV}"

if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi

echo "=== Helixer postprocess/QC: final ${HAP_SET} ==="
echo "Genome : ${GENOME}"
echo "GFF3   : ${GFF3}"
if [ -n "${TE_GFF3}" ]; then
    require_file "${TE_GFF3}" "optional TE GFF3"
    echo "TE GFF3: ${TE_GFF3}"
else
    echo "TE GFF3: not used"
fi
echo "Out    : ${OUT_DIR}"
echo

if [ -n "${TE_GFF3}" ]; then
    "${PYTHON_BIN}" "${SCRIPT}" \
        --genome "${GENOME}" \
        --gff3 "${GFF3}" \
        --te-gff3 "${TE_GFF3}" \
        --out-dir "${OUT_DIR}" \
        --prefix "${PREFIX}" \
        --short-cds 300 \
        --gene-te-fraction 0.80 \
        --cds-te-fraction 0.50
else
    "${PYTHON_BIN}" "${SCRIPT}" \
        --genome "${GENOME}" \
        --gff3 "${GFF3}" \
        --out-dir "${OUT_DIR}" \
        --prefix "${PREFIX}" \
        --short-cds 300 \
        --gene-te-fraction 0.80 \
        --cds-te-fraction 0.50
fi

echo "=== Expected outputs ==="
echo "${OUT_DIR}/${PREFIX}.genes.gff3"
echo "${OUT_DIR}/${PREFIX}.cds.fa"
echo "${OUT_DIR}/${PREFIX}.protein.fa"
echo "${OUT_DIR}/${PREFIX}.mrna.fa"
echo "${OUT_DIR}/${PREFIX}.gene_model_qc.tsv"
echo "${OUT_DIR}/${PREFIX}.postprocess_summary.tsv"
