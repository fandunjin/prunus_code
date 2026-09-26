#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${PROJECT_ROOT}/data/pmal"
INPUT_FASTA_GZ="${DATA_DIR}/GCA_048569315.1_ASM4856931v1_genomic.fna.gz"
ANNOTATION_DIR="${DATA_DIR}/annotation"
WORK_DIR="${ANNOTATION_DIR}/helixer/work"
OUT_DIR="${ANNOTATION_DIR}/helixer"
GENOME="${WORK_DIR}/pmal.genomic.fa"
OUT_GFF="${OUT_DIR}/pmal.helixer.gff3"
STATS_DIR="${OUT_DIR}/stats"

CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
SINGULARITY_ENV="singularity"
HELIXER_CONTAINER_DIR="${ROOT_DIR}/software/helixer_container"
HELIXER_SIF="${HELIXER_SIF:-}"
MODEL_DIR="${ROOT_DIR}/annotation/genes/helixer_models"
MODEL_FILE="${MODEL_DIR}/land_plant/land_plant_v0.3_a_0080.h5"
SINGULARITY_TMPDIR="${ANNOTATION_DIR}/tmp/singularity"
GFF_STATS_SCRIPT="${ROOT_DIR}/scripts/gff3_basic_stats.py"
LINEAGE="land_plant"
SPECIES="Prunus_mahaleb"

require_file() {
    local path="$1"
    local label="$2"
    if [ ! -s "${path}" ]; then
        echo "ERROR: ${label} does not exist or is empty: ${path}" >&2
        exit 1
    fi
}

require_file "${INPUT_FASTA_GZ}" "P. mahaleb genome"
require_file "${CONDA_SH}" "Conda initialization script"
require_file "${GFF_STATS_SCRIPT}" "GFF3 statistics script"
mkdir -p "${WORK_DIR}" "${OUT_DIR}" "${STATS_DIR}" "${MODEL_DIR}" "${SINGULARITY_TMPDIR}"

if [ ! -s "${GENOME}" ]; then
    gzip -t "${INPUT_FASTA_GZ}"
    gzip -dc "${INPUT_FASTA_GZ}" > "${GENOME}.tmp"
    mv "${GENOME}.tmp" "${GENOME}"
fi
require_file "${GENOME}" "uncompressed P. mahaleb genome"

if [ -z "${HELIXER_SIF}" ]; then
    if [ -s "${HELIXER_CONTAINER_DIR}/helixer-docker_latest.sif" ]; then
        HELIXER_SIF="${HELIXER_CONTAINER_DIR}/helixer-docker_latest.sif"
    elif [ -s "${HELIXER_CONTAINER_DIR}/helixer-docker_helixer_v0.3.6_cuda_12.2.2-cudnn8.sif" ]; then
        HELIXER_SIF="${HELIXER_CONTAINER_DIR}/helixer-docker_helixer_v0.3.6_cuda_12.2.2-cudnn8.sif"
    else
        echo "ERROR: no Helixer SIF found in ${HELIXER_CONTAINER_DIR}" >&2
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

echo "=== Helixer: ${SPECIES} ==="
echo "Genome : ${GENOME}"
echo "GFF3   : ${OUT_GFF}"
echo "SIF    : ${HELIXER_SIF}"
nvidia-smi

HELIXER_ARGS=(
    Helixer.py
    --lineage "${LINEAGE}"
    --fasta-path "/work/data/pmal/annotation/helixer/work/pmal.genomic.fa"
    --species "${SPECIES}"
    --gff-output-path "/work/data/pmal/annotation/helixer/pmal.helixer.gff3"
)
if [ -s "${MODEL_FILE}" ]; then
    HELIXER_ARGS+=(--downloaded-model-path "/work/annotation/genes/helixer_models")
fi

if [ ! -s "${OUT_GFF}" ]; then
    singularity run --nv -B "${ROOT_DIR}:/work" "${HELIXER_SIF}" "${HELIXER_ARGS[@]}"
else
    echo "Existing GFF3 found; skipping Helixer: ${OUT_GFF}"
fi

require_file "${OUT_GFF}" "Helixer GFF3"
"${PYTHON_BIN}" "${GFF_STATS_SCRIPT}" --gff3 "${OUT_GFF}" --prefix pmal_helixer --out-dir "${STATS_DIR}"
echo "Helixer completed: ${OUT_GFF}"
