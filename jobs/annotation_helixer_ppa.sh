#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail
# Helixer is genome-only. PPA transcript evidence is incorporated by the
# separate jobs/ppa_annotation_eviann.sh refinement job.

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/data/ppa}"
INPUT_FASTA_GZ="${INPUT_FASTA_GZ:-${DATA_DIR}/GCA_964059315.1_drPruPadu1.1_genomic.fna.gz}"
ROOT_PARENT_DIR="${ROOT_PARENT_DIR:-$(dirname "${ROOT_DIR}")}"
PPA_EDTA_DIR="${PPA_EDTA_DIR:-${ROOT_PARENT_DIR}/prunus_old/analysis/prunus_pan_te/discovery/ppa}"
PPA_EDTA_MASKED="${PPA_EDTA_MASKED:-${PPA_EDTA_DIR}/ppa.fa.mod.MAKER.masked}"
ANNOTATION_DIR="${ANNOTATION_DIR:-${DATA_DIR}/annotation}"
WORK_DIR="${WORK_DIR:-${ANNOTATION_DIR}/helixer/work}"
OUT_DIR="${OUT_DIR:-${ANNOTATION_DIR}/helixer}"
GENOME="${GENOME:-${PPA_EDTA_MASKED}}"
OUT_GFF="${OUT_GFF:-${OUT_DIR}/ppa.helixer.gff3}"
STATS_DIR="${STATS_DIR:-${OUT_DIR}/stats}"

CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
SINGULARITY_ENV="${SINGULARITY_ENV:-singularity}"
HELIXER_CONTAINER_DIR="${HELIXER_CONTAINER_DIR:-${ROOT_DIR}/software/helixer_container}"
HELIXER_SIF="${HELIXER_SIF:-}"
MODEL_DIR="${MODEL_DIR:-${ROOT_DIR}/annotation/genes/helixer_models}"
MODEL_FILE="${MODEL_FILE:-${MODEL_DIR}/land_plant/land_plant_v0.3_a_0080.h5}"
SINGULARITY_TMPDIR="${SINGULARITY_TMPDIR:-${ANNOTATION_DIR}/tmp/singularity}"
GFF_STATS_SCRIPT="${GFF_STATS_SCRIPT:-${ROOT_DIR}/scripts/gff3_basic_stats.py}"
LINEAGE="${LINEAGE:-land_plant}"
SPECIES="${SPECIES:-Prunus_padus}"

require_file() {
    local path="$1"
    local label="$2"
    if [ ! -s "${path}" ]; then
        echo "ERROR: ${label} does not exist or is empty: ${path}" >&2
        exit 1
    fi
}

require_file "${CONDA_SH}" "Conda initialization script"
require_file "${GFF_STATS_SCRIPT}" "GFF3 statistics script"
mkdir -p "${WORK_DIR}" "${OUT_DIR}" "${STATS_DIR}" "${MODEL_DIR}" "${SINGULARITY_TMPDIR}"

if [ ! -s "${GENOME}" ]; then
    if [[ "${GENOME}" != "${WORK_DIR}/"* ]]; then
        echo "ERROR: PPA EDTA-masked genome is missing: ${GENOME}" >&2
        echo "Set PPA_EDTA_MASKED/GENOME to a valid *.MAKER.masked FASTA." >&2
        exit 1
    fi
    require_file "${INPUT_FASTA_GZ}" "P. padus genome"
    gzip -t "${INPUT_FASTA_GZ}"
    gzip -dc "${INPUT_FASTA_GZ}" > "${GENOME}.tmp"
    mv "${GENOME}.tmp" "${GENOME}"
fi
require_file "${GENOME}" "P. padus EDTA-masked genome"

# ROOT_DIR/../prunus_old is outside ROOT_DIR.  Canonicalize before deciding
# whether the file can use the /work bind; otherwise the literal path prefix
# makes /work/../prunus_old look like an in-container path that does not exist.
if command -v readlink >/dev/null 2>&1; then
    GENOME_CANONICAL="$(readlink -f -- "${GENOME}")"
else
    GENOME_CANONICAL="${GENOME}"
fi
require_file "${GENOME_CANONICAL}" "canonical P. padus EDTA-masked genome"
GENOME="${GENOME_CANONICAL}"

SINGULARITY_BINDS=(-B "${ROOT_DIR}:/work")
if [[ "${GENOME}" == "${ROOT_DIR}/"* ]]; then
    GENOME_REL="${GENOME#${ROOT_DIR}/}"
    GENOME_CONTAINER="/work/${GENOME_REL}"
else
    GENOME_BASENAME="$(basename "${GENOME}")"
    GENOME_DIR="$(dirname "${GENOME}")"
    SINGULARITY_BINDS+=(-B "${GENOME_DIR}:/ppa_input")
    GENOME_CONTAINER="/ppa_input/${GENOME_BASENAME}"
fi
OUT_REL="${OUT_GFF#${ROOT_DIR}/}"
[[ "${OUT_REL}" != "${OUT_GFF}" ]] || { echo "OUT_GFF must be under ROOT_DIR for Singularity binding" >&2; exit 1; }

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
    --fasta-path "${GENOME_CONTAINER}"
    --species "${SPECIES}"
    --gff-output-path "/work/${OUT_REL}"
)
if [ -s "${MODEL_FILE}" ]; then
    HELIXER_ARGS+=(--downloaded-model-path "/work/annotation/genes/helixer_models")
fi

if [ ! -s "${OUT_GFF}" ]; then
    singularity run --nv "${SINGULARITY_BINDS[@]}" "${HELIXER_SIF}" "${HELIXER_ARGS[@]}"
else
    echo "Existing GFF3 found; skipping Helixer: ${OUT_GFF}"
fi

require_file "${OUT_GFF}" "Helixer GFF3"
"${PYTHON_BIN}" "${GFF_STATS_SCRIPT}" --gff3 "${OUT_GFF}" --prefix ppa_helixer --out-dir "${STATS_DIR}"
echo "Helixer completed: ${OUT_GFF}"
