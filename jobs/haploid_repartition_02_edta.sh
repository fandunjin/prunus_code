#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/haploid_repartition/inputs}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/annotation/repartition/repeats}"
SPECIES="${SPECIES:-pzi}"
LABEL="${LABEL:-hapA1}"
THREADS="${THREADS:-24}"
EDTA_ENV="${EDTA_ENV:-EDTA}"
EDTA_FORCE="${EDTA_FORCE:-1}"
EDTA_OVERWRITE="${EDTA_OVERWRITE:-0}"
EDTA_CDS="${EDTA_CDS:-${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_cds_from_genomic.fna}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing or empty $2: $1"; }

case "${SPECIES}:${LABEL}" in
    pzi:hapA1|pzi:hapA2|pzi:hapB1|pzi:hapB2|psp:hapA1|psp:hapA2|psp:hapB1|psp:hapB2|ppa:hapA|ppa:hapB) ;;
    *) die "Unsupported SPECIES/LABEL=${SPECIES}/${LABEL}" ;;
esac

INPUT="${INPUT_ROOT}/${SPECIES}/genome/${SPECIES}_${LABEL}.fa"
OUT_DIR="${OUT_ROOT}/${SPECIES}/${LABEL}_edta"
WORK_FASTA="${OUT_DIR}/${SPECIES}_${LABEL}.fa"
require_file "${INPUT}" "${SPECIES} ${LABEL} genome input"
require_file "${CONDA_SH}" "conda initialization"
mkdir -p "${OUT_DIR}" "${ROOT_DIR}/logs"

if [[ ! -s "${WORK_FASTA}" ]]; then
    cp -p "${INPUT}" "${WORK_FASTA}"
elif ! cmp -s "${INPUT}" "${WORK_FASTA}"; then
    die "Existing EDTA work FASTA differs from input; move ${OUT_DIR} aside before rerun"
fi

source "${CONDA_SH}"
conda activate "${EDTA_ENV}"
cd "${OUT_DIR}"

echo "=== Independent EDTA: ${SPECIES} ${LABEL} ==="
echo "Input : ${WORK_FASTA}"
echo "Mode  : EDTA.pl --step all (no panEDTA/uniform library)"
echo "CDS   : ${EDTA_CDS:-not used}"

if [[ ! -s "${WORK_FASTA}.mod.EDTA.TEanno.gff3" ]]; then
    cmd=(EDTA.pl --genome "${WORK_FASTA}" --species others --step all --anno 1 --sensitive 1 --force "${EDTA_FORCE}" --overwrite "${EDTA_OVERWRITE}" --threads "${THREADS}")
    [[ -s "${EDTA_CDS}" ]] && cmd+=(--cds "${EDTA_CDS}")
    "${cmd[@]}"
else
    echo "[SKIP] Existing EDTA result: ${WORK_FASTA}.mod.EDTA.TEanno.gff3"
fi

for output in \
    "${WORK_FASTA}.mod.EDTA.TEanno.gff3" \
    "${WORK_FASTA}.mod.EDTA.TEanno.sum" \
    "${WORK_FASTA}.mod.EDTA.TElib.fa" \
    "${WORK_FASTA}.mod.MAKER.masked"; do
    require_file "${output}" "EDTA output"
done

printf 'species\tlabel\tgenome\tte_gff3\tte_summary\tte_library\tmasked_genome\tstatus\n' > "${OUT_DIR}/edta.completed.tsv"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\tpass\n' "${SPECIES}" "${LABEL}" "${WORK_FASTA}" "${WORK_FASTA}.mod.EDTA.TEanno.gff3" "${WORK_FASTA}.mod.EDTA.TEanno.sum" "${WORK_FASTA}.mod.EDTA.TElib.fa" "${WORK_FASTA}.mod.MAKER.masked" >> "${OUT_DIR}/edta.completed.tsv"
echo "[OK] EDTA outputs: ${OUT_DIR}"
