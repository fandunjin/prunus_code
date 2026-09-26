#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/haploid_repartition/inputs}"
GENE_ROOT="${GENE_ROOT:-${ROOT_DIR}/annotation/repartition/genes}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
SINGULARITY_ENV="${SINGULARITY_ENV:-singularity}"
HELIXER_SIF="${HELIXER_SIF:-${ROOT_DIR}/software/helixer_container/helixer-docker_helixer_v0.3.6_cuda_12.2.2-cudnn8.sif}"
MODEL_DIR="${MODEL_DIR:-${GENE_ROOT}/helixer_models}"
LINEAGE="${LINEAGE:-land_plant}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing or empty $2: $1"; }
case "${SPECIES}:${LABEL}" in
    pzi:hapA1|pzi:hapA2|pzi:hapB1|pzi:hapB2|psp:hapA1|psp:hapA2|psp:hapB1|psp:hapB2|ppa:hapA|ppa:hapB) ;;
    *) die "Unsupported SPECIES/LABEL=${SPECIES}/${LABEL}" ;;
esac

EDTA_DIR="${EDTA_DIR:-${EDTA_INPUT_ROOT}/${SPECIES}/edta/${LABEL}}"
INPUT_FASTA="${EDTA_DIR}/${SPECIES}_${LABEL}.fa.mod.MAKER.masked"
OUT_DIR="${GENE_ROOT}/${SPECIES}/${LABEL}_helixer"
OUT_GFF="${OUT_DIR}/${SPECIES}_${LABEL}.edta_masked.helixer.gff3"
EDTA_REL="${EDTA_DIR#${ROOT_DIR}/}"
OUT_REL="${OUT_DIR#${ROOT_DIR}/}"
MODEL_REL="${MODEL_DIR#${ROOT_DIR}/}"
[[ "${EDTA_REL}" != "${EDTA_DIR}" && "${OUT_REL}" != "${OUT_DIR}" ]] || die "EDTA_DIR and GENE_ROOT must be under ROOT_DIR for Singularity binding"
require_file "${INPUT_FASTA}" "EDTA masked genome"
require_file "${CONDA_SH}" "conda initialization"
require_file "${HELIXER_SIF}" "Helixer SIF"
mkdir -p "${OUT_DIR}" "${MODEL_DIR}" "${ROOT_DIR}/annotation/repartition/tmp/singularity/${SPECIES}_${LABEL}"

source "${CONDA_SH}"
conda activate "${SINGULARITY_ENV}"
export SINGULARITY_TMPDIR="${ROOT_DIR}/annotation/repartition/tmp/singularity/${SPECIES}_${LABEL}"
command -v singularity >/dev/null 2>&1 || die "singularity is unavailable"

echo "=== Helixer: ${SPECIES} ${LABEL} ==="
if [[ ! -s "${OUT_GFF}" ]]; then
    cmd=(singularity run --nv -B "${ROOT_DIR}:/work" "${HELIXER_SIF}" Helixer.py --lineage "${LINEAGE}" --fasta-path "/work/${EDTA_REL}/${SPECIES}_${LABEL}.fa.mod.MAKER.masked" --species "${SPECIES}_${LABEL}_edta_masked" --gff-output-path "/work/${OUT_REL}/${SPECIES}_${LABEL}.edta_masked.helixer.gff3")
    [[ -s "${MODEL_DIR}/${LINEAGE}/land_plant_v0.3_a_0080.h5" ]] && cmd+=(--downloaded-model-path "/work/${MODEL_REL}")
    "${cmd[@]}"
else
    echo "[SKIP] Existing Helixer GFF3: ${OUT_GFF}"
fi
require_file "${OUT_GFF}" "Helixer GFF3"
echo "[OK] Helixer output: ${OUT_GFF}"
