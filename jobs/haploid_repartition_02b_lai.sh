#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

# Run once per revised genome, for example:
#   env SPECIES=pzi LABEL=hapA1 jsub < jobs/haploid_repartition_02b_lai.sh
# EDTA repartition must have completed in haploid_repartition/inputs/*/edta first.

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/haploid_repartition/inputs}"
LAI_ROOT="${LAI_ROOT:-${ROOT_DIR}/annotation/repartition/lai}"
SPECIES="${SPECIES:-pzi}"
LABEL="${LABEL:-hapA1}"
LAI_THREADS="${LAI_THREADS:-24}"
LAI_METHOD="${LAI_METHOD:-auto}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
LAI_ENV="${LAI_ENV:-LAI}"
PREP_SCRIPT="${PREP_SCRIPT:-${ROOT_DIR}/scripts/prepare_lai_inputs.py}"
GFF2BED_SCRIPT="${GFF2BED_SCRIPT:-}"
FORCE_PREP="${FORCE_PREP:-0}"
FORCE_LAI="${FORCE_LAI:-0}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() {
    local path="$1" label="$2"
    [[ -s "${path}" ]] || die "Missing or empty ${label}: ${path}"
}

case "${SPECIES}:${LABEL}" in
    pzi:hapA1|pzi:hapA2|pzi:hapB1|pzi:hapB2|psp:hapA1|psp:hapA2|psp:hapB1|psp:hapB2|ppa:hapA|ppa:hapB) ;;
    *) die "Unsupported SPECIES/LABEL=${SPECIES}/${LABEL}" ;;
esac

EDTA_DIR="${EDTA_DIR:-${EDTA_INPUT_ROOT}/${SPECIES}/edta/${LABEL}}"
LAI_DIR="${LAI_ROOT}/${SPECIES}/${LABEL}_lai"
PREFIX="${SPECIES}_${LABEL}"
TE_GFF3="${EDTA_DIR}/${PREFIX}.fa.mod.EDTA.TEanno.gff3"

require_file "${CONDA_SH}" "conda initialization script"
require_file "${PREP_SCRIPT}" "LAI input preparation script"
require_file "${TE_GFF3}" "EDTA TE annotation GFF3"
require_file "${EDTA_DIR}/${PREFIX}.fa.mod.EDTA.TEanno.sum" "EDTA TE summary"
[[ "${FORCE_PREP}" == "0" || "${FORCE_PREP}" == "1" ]] || die "FORCE_PREP must be 0 or 1"
[[ "${FORCE_LAI}" == "0" || "${FORCE_LAI}" == "1" ]] || die "FORCE_LAI must be 0 or 1"
mkdir -p "${LAI_DIR}" "${ROOT_DIR}/logs"

source "${CONDA_SH}"
conda activate "${LAI_ENV}"
command -v python3 >/dev/null 2>&1 && PYTHON_BIN=python3 || PYTHON_BIN=python
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || die "Python is unavailable in ${LAI_ENV}"
PYTHON_VERSION="$("${PYTHON_BIN}" -c 'import sys; print(".".join(map(str, sys.version_info[:3])))')" \
    || die "Unable to query Python version: ${PYTHON_BIN}"
"${PYTHON_BIN}" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 6) else 1)' \
    || die "Python >= 3.6 is required; found ${PYTHON_VERSION}"
command -v LAI >/dev/null 2>&1 || die "LAI is unavailable in ${LAI_ENV}"

if [[ -n "${GFF2BED_SCRIPT}" && -s "${GFF2BED_SCRIPT}" ]]; then
    GFF2BED_RESOLVED="${GFF2BED_SCRIPT}"
elif [[ -s "${CONDA_ROOT}/envs/EDTA/share/EDTA/bin/gff2bed.pl" ]]; then
    GFF2BED_RESOLVED="${CONDA_ROOT}/envs/EDTA/share/EDTA/bin/gff2bed.pl"
elif [[ -s "${CONDA_ROOT}/envs/EDTA/share/LTR_retriever/bin/gff2bed.pl" ]]; then
    GFF2BED_RESOLVED="${CONDA_ROOT}/envs/EDTA/share/LTR_retriever/bin/gff2bed.pl"
elif command -v gff2bed.pl >/dev/null 2>&1; then
    GFF2BED_RESOLVED="$(command -v gff2bed.pl)"
else
    die "gff2bed.pl was not found; set GFF2BED_SCRIPT=/path/to/gff2bed.pl"
fi

echo "=== LAI: ${SPECIES} ${LABEL} ==="
echo "EDTA dir : ${EDTA_DIR}"
echo "LAI dir  : ${LAI_DIR}"
echo "Prefix   : ${PREFIX}"
echo "Method   : ${LAI_METHOD}"
echo "Threads  : ${LAI_THREADS}"
echo "gff2bed  : ${GFF2BED_RESOLVED}"
echo "Force prep: ${FORCE_PREP}"
echo "Force LAI : ${FORCE_LAI}"

if [[ "${FORCE_PREP}" == "1" ]]; then
    rm -f "${LAI_DIR}/${PREFIX}.lai_inputs.env" "${LAI_DIR}/${PREFIX}.lai_inputs.tsv"
fi
if [[ "${FORCE_LAI}" == "1" ]]; then
    rm -f "${LAI_DIR}/${PREFIX}.all_ltr.out.LAI" \
        "${LAI_DIR}/${PREFIX}.all_ltr.out.LAI."*
fi

"${PYTHON_BIN}" "${PREP_SCRIPT}" \
    --edta-dir "${EDTA_DIR}" \
    --out-dir "${LAI_DIR}" \
    --prefix "${PREFIX}" \
    --gff2bed "${GFF2BED_RESOLVED}" \
    --method "${LAI_METHOD}"

require_file "${LAI_DIR}/${PREFIX}.lai_inputs.env" "LAI input environment file"
# shellcheck disable=SC1090
source "${LAI_DIR}/${PREFIX}.lai_inputs.env"
require_file "${GENOME}" "LAI genome"
require_file "${INTACT}" "LAI intact LTR pass list"
require_file "${ALL}" "LAI all-LTR RepeatMasker output"

cd "${LAI_DIR}"
echo "Genome : ${GENOME}"
echo "Intact : ${INTACT}"
echo "All    : ${ALL}"
echo "Result : ${LAI_RESULT}"

if [[ ! -s "${LAI_RESULT}" ]]; then
    LAI -genome "${GENOME}" -intact "${INTACT}" -all "${ALL}" -q -t "${LAI_THREADS}"
else
    echo "[SKIP] Existing LAI result: ${LAI_RESULT}"
fi
require_file "${LAI_RESULT}" "LAI result"
if [[ -e "${LAI_DIR}/${PREFIX}.all_ltr.out.LAI.LTR.fa" &&
      ! -s "${LAI_DIR}/${PREFIX}.all_ltr.out.LAI.LTR.fa" ]]; then
    die "${SPECIES}/${LABEL}: LAI produced an empty LTR FASTA; check genome/LAI IDs"
fi

SUMMARY="${LAI_ROOT}/${SPECIES}/${LABEL}.lai.completed.tsv"
{
    printf 'species\tlabel\tedta_dir\tte_gff3\tlai_dir\tgenome\tintact\tall_ltr\tlai_result\tmethod\tstatus\n'
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tpass\n' \
        "${SPECIES}" "${LABEL}" "${EDTA_DIR}" "${TE_GFF3}" "${LAI_DIR}" \
        "${GENOME}" "${INTACT}" "${ALL}" "${LAI_RESULT}" "${LAI_INPUT_METHOD}"
} > "${SUMMARY}.tmp"
mv "${SUMMARY}.tmp" "${SUMMARY}"

echo "[OK] LAI result: ${LAI_RESULT}"
echo "[OK] LAI summary: ${SUMMARY}"
