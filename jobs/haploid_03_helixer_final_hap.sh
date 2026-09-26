#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
OUT_DIR="${ROOT_DIR}/annotation/genes/final_${HAP_SET}_helixer"
OUT_GFF="${OUT_DIR}/pzi_${HAP_SET}.helixer.gff3"
STATS_DIR="${OUT_DIR}/stats"
CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
SINGULARITY_ENV="singularity"
HELIXER_CONTAINER_DIR="${ROOT_DIR}/software/helixer_container"
HELIXER_SIF="${HELIXER_SIF:-}"
MODEL_DIR="${ROOT_DIR}/annotation/genes/helixer_models"
MODEL_FILE="${MODEL_DIR}/land_plant/land_plant_v0.3_a_0080.h5"
SINGULARITY_TMPDIR="${ROOT_DIR}/annotation/tmp/singularity"
LINEAGE="land_plant"
SPECIES="pzi_${HAP_SET}"
GFF_STATS_SCRIPT="${ROOT_DIR}/scripts/gff3_basic_stats.py"

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

require_file "${INPUT_FASTA}" "final phased hap FASTA"
require_file "${CONDA_SH}" "Conda initialization script"
require_file "${GFF_STATS_SCRIPT}" "GFF3 stats script"
mkdir -p "${OUT_DIR}" "${STATS_DIR}" "${ROOT_DIR}/logs" "${MODEL_DIR}" "${SINGULARITY_TMPDIR}"

if [ -z "${HELIXER_SIF}" ]; then
    if [ -s "${HELIXER_CONTAINER_DIR}/helixer-docker_latest.sif" ]; then
        HELIXER_SIF="${HELIXER_CONTAINER_DIR}/helixer-docker_latest.sif"
    elif [ -s "${HELIXER_CONTAINER_DIR}/helixer-docker_helixer_v0.3.6_cuda_12.2.2-cudnn8.sif" ]; then
        HELIXER_SIF="${HELIXER_CONTAINER_DIR}/helixer-docker_helixer_v0.3.6_cuda_12.2.2-cudnn8.sif"
    else
        echo "ERROR: could not find Helixer .sif in ${HELIXER_CONTAINER_DIR}" >&2
        echo "Set HELIXER_SIF=/path/to/helixer.sif before submitting." >&2
        exit 1
    fi
fi
require_file "${HELIXER_SIF}" "Helixer Singularity image"

source "${CONDA_SH}"
conda activate "${SINGULARITY_ENV}"
export SINGULARITY_TMPDIR

if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi

echo "=== Helixer gene prediction: final ${HAP_SET} ==="
echo "Input  : ${INPUT_FASTA}"
echo "Output : ${OUT_GFF}"
echo "Lineage: ${LINEAGE}"
echo "Species: ${SPECIES}"
echo "SIF    : ${HELIXER_SIF}"
if [ -s "${MODEL_FILE}" ]; then
    echo "Model  : ${MODEL_FILE}"
    echo "Mode   : project downloaded-model-path"
else
    echo "Model  : Helixer default model path"
    echo "Mode   : no downloaded-model-path"
fi
nvidia-smi
singularity --version
echo

if [ ! -s "${OUT_GFF}" ]; then
    if [ -s "${MODEL_FILE}" ]; then
        singularity run --nv \
            -B "${ROOT_DIR}:/work" \
            "${HELIXER_SIF}" \
            Helixer.py \
            --lineage "${LINEAGE}" \
            --fasta-path "/work/haploid/final_haplotypes/pzi_${HAP_SET}.fa" \
            --species "${SPECIES}" \
            --downloaded-model-path "/work/annotation/genes/helixer_models" \
            --gff-output-path "/work/annotation/genes/final_${HAP_SET}_helixer/pzi_${HAP_SET}.helixer.gff3"
    else
        singularity run --nv \
            -B "${ROOT_DIR}:/work" \
            "${HELIXER_SIF}" \
            Helixer.py \
            --lineage "${LINEAGE}" \
            --fasta-path "/work/haploid/final_haplotypes/pzi_${HAP_SET}.fa" \
            --species "${SPECIES}" \
            --gff-output-path "/work/annotation/genes/final_${HAP_SET}_helixer/pzi_${HAP_SET}.helixer.gff3"
    fi
else
    echo "Existing Helixer GFF3 found, skipping prediction: ${OUT_GFF}"
fi

require_file "${OUT_GFF}" "Helixer final hap GFF3"
"${PYTHON_BIN}" "${GFF_STATS_SCRIPT}" --gff3 "${OUT_GFF}" --prefix "${SPECIES}" --out-dir "${STATS_DIR}"

echo "Helixer final ${HAP_SET} completed."
