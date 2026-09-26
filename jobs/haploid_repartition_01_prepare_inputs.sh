#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/haploid_repartition}"
PZI_COPY1="${PZI_COPY1:-${ROOT_DIR}/final_assembly/haplotype_provenance/pzi_2n32.copy1.hap_provenance.chr16.fa}"
PZI_COPY2="${PZI_COPY2:-${ROOT_DIR}/final_assembly/haplotype_provenance/pzi_2n32.copy2.hap_provenance.chr16.fa}"
# The complete-assembly mode is the default because current PZI EviAnn
# coordinates are defined on pzi32.  Set PZI_SOURCE_MODE=copy only for
# phase-specific genome analyses whose coordinates use the copy FASTAs.
PZI_SOURCE_MODE="${PZI_SOURCE_MODE:-full32}"
PZI_FULL32="${PZI_FULL32:-${ROOT_DIR}/final_assembly/pzi_2n32.chrom/pzi_2n32.chromosomes.fa}"
if [[ ! -s "${PZI_FULL32}" && -s "${ROOT_DIR}/final_assembly/pzi_2n32.chromosomes.fa" ]]; then
    PZI_FULL32="${ROOT_DIR}/final_assembly/pzi_2n32.chromosomes.fa"
fi
# Current copy1/copy2 provenance map. It links the normalized chr01--chr16
# IDs in the two FASTAs to their source chromosomes in the 2n=32 assembly.
PZI_MAP="${PZI_MAP:-${ROOT_DIR}/final_assembly/haplotype_provenance/pzi_2n32.corrected_homolog_pairs.tsv}"
PPA_GENOME="${PPA_GENOME:-${ROOT_DIR}/data/ppa/GCA_964059315.1_drPruPadu1.1_genomic.fna.gz}"
PPA_ASSIGNMENT="${PPA_ASSIGNMENT:-${ROOT_DIR}/haploid_repartition/assignments/ppa_A_B.rename_map.tsv}"
PSP_GENOME="${PSP_GENOME:-${ROOT_DIR}/data/psp/haploid/genome/psp_A1A2B1B2.renamed.chr32.fa}"
PSP_ASSIGNMENT="${PSP_ASSIGNMENT:-${ROOT_DIR}/data/psp/haploid/psp_A1A2B1B2.rename_map.tsv}"
PZI_ONLY="${PZI_ONLY:-1}"
SCRIPT="${ROOT_DIR}/scripts/prepare_repartitioned_subgenomes.py"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
PY_ENV="${PY_ENV:-haphic}"
FORCE="${FORCE:-0}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing or empty $2: $1"; }

case "${PZI_SOURCE_MODE}" in
    copy)
        require_file "${PZI_COPY1}" "PZI copy1 chr16 FASTA"
        require_file "${PZI_COPY2}" "PZI copy2 chr16 FASTA"
        ;;
    full32)
        require_file "${PZI_FULL32}" "complete PZI 2n=32 FASTA"
        ;;
    *)
        die "PZI_SOURCE_MODE must be copy or full32 (found: ${PZI_SOURCE_MODE})"
        ;;
esac
require_file "${PZI_MAP}" "current PZI corrected homolog provenance map"
if [[ "$(basename "${PZI_MAP}")" == "chromosome_rename_map.tsv" ]]; then
    die "The retired haploid_rename/chromosome_rename_map.tsv cannot be used; provide final_assembly/haplotype_provenance/pzi_2n32.corrected_homolog_pairs.tsv"
fi
if [[ "${PZI_ONLY}" != "1" ]]; then
    require_file "${PPA_GENOME}" "PPA primary FASTA"
    require_file "${PPA_ASSIGNMENT}" "PPA assignment table"
    require_file "${PSP_GENOME}" "PSP renamed chr32 FASTA"
    require_file "${PSP_ASSIGNMENT}" "PSP assignment table"
fi
require_file "${SCRIPT}" "preparation script"
require_file "${CONDA_SH}" "conda initialization"

mkdir -p "${ROOT_DIR}/logs"
source "${CONDA_SH}"
conda activate "${PY_ENV}"
command -v python >/dev/null 2>&1 || die "python is unavailable in ${PY_ENV}"

args=(
    --pzi-map "${PZI_MAP}"
    --out-root "${OUT_ROOT}"
)
if [[ "${PZI_SOURCE_MODE}" == "full32" ]]; then
    args+=(--pzi-full32 "${PZI_FULL32}")
else
    args+=(--pzi-copy1 "${PZI_COPY1}" --pzi-copy2 "${PZI_COPY2}")
fi
if [[ "${PZI_ONLY}" != "1" ]]; then
    args+=(
        --ppa-genome "${PPA_GENOME}"
        --ppa-assignment "${PPA_ASSIGNMENT}"
        --psp-genome "${PSP_GENOME}"
        --psp-assignment "${PSP_ASSIGNMENT}"
    )
fi
[[ "${FORCE}" == "1" ]] && args+=(--force)

echo "=== Revised PZI/PPA/PSP A/B partition ==="
echo "A copy1: chr01 chr02 chr05 chr06 chr07 chr10 chr12 chr14"
echo "B copy1: chr03 chr04 chr08 chr09 chr11 chr13 chr15 chr16"
echo "Chromosome IDs: C1A-C8A and C1B-C8B"
echo "PZI source mode: ${PZI_SOURCE_MODE}"
if [[ "${PZI_SOURCE_MODE}" == "full32" ]]; then
    echo "Complete 2n=32 FASTA: ${PZI_FULL32}"
else
    echo "Copy1 FASTA: ${PZI_COPY1}"
    echo "Copy2 FASTA: ${PZI_COPY2}"
fi
echo "Provenance map: ${PZI_MAP}"
echo "Mode: $([[ \"${PZI_ONLY}\" == \"1\" ]] && echo PZI-only || echo PZI+PPA+PSP)"
echo "Output: ${OUT_ROOT}"
python "${SCRIPT}" "${args[@]}"

echo "[OK] Inputs and maps written: ${OUT_ROOT}"
