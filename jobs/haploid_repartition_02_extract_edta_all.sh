#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

# Split the already completed EDTA annotation for each complete assembly.
# This job does not invoke EDTA and therefore preserves the original TE calls.
ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/haploid_repartition/inputs}"
SPECIES_SET="${SPECIES_SET:-pzi psp ppa}"
FORCE="${FORCE:-0}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
PYTHON_ENV="${PYTHON_ENV:-haphic}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
SCRIPT="${SCRIPT:-${ROOT_DIR}/scripts/repartition_edta_results.py}"
PZI_SOURCE_DIR="${PZI_SOURCE_DIR:-${ROOT_DIR}/annotation/chr32_edta}"
PSP_SOURCE_DIR="${PSP_SOURCE_DIR:-${ROOT_DIR}/data/psp/annotation/repeats/chr32_edta}"
PPA_SOURCE_DIR="${PPA_SOURCE_DIR:-${ROOT_DIR}/../prunus_old/analysis/prunus_pan_te/discovery/ppa}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing or empty $2: $1"; }

require_file "${SCRIPT}" "EDTA repartition script"
require_file "${CONDA_SH}" "Conda initialization script"
source "${CONDA_SH}"
command -v conda >/dev/null 2>&1 || die "Conda is unavailable after sourcing: ${CONDA_SH}"
conda activate "${PYTHON_ENV}"
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || die "Python is unavailable: ${PYTHON_BIN}"
PYTHON_VERSION="$("${PYTHON_BIN}" -c 'import sys; print(".".join(map(str, sys.version_info[:3])))')" \
    || die "Unable to query Python version: ${PYTHON_BIN}"
"${PYTHON_BIN}" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 6) else 1)' \
    || die "Python >= 3.6 is required; found ${PYTHON_VERSION}"
mkdir -p "${INPUTS}" "${ROOT_DIR}/logs"

read -r -a SPECIES_LIST <<< "${SPECIES_SET}"
[[ "${#SPECIES_LIST[@]}" -gt 0 ]] || die "SPECIES_SET is empty"
for species in "${SPECIES_LIST[@]}"; do
    case "${species}" in
        pzi) require_file "${PZI_SOURCE_DIR}/pzi_2n32.chromosomes.fa" "PZI complete EDTA FASTA" ;;
        psp) require_file "${PSP_SOURCE_DIR}/psp_2n32.chromosomes.fa" "PSP complete EDTA FASTA" ;;
        ppa) require_file "${PPA_SOURCE_DIR}/ppa.fa" "PPA complete EDTA FASTA" ;;
        *) die "Unknown species in SPECIES_SET: ${species} (choose pzi, psp or ppa)" ;;
    esac
    if [[ "${species}" == "ppa" ]]; then
        haps=(hapA hapB)
    else
        haps=(hapA1 hapA2 hapB1 hapB2)
    fi
    for hap in "${haps[@]}"; do
        require_file "${INPUTS}/${species}/genome/${species}_${hap}.fa" "${species}/${hap} partitioned genome"
    done
done

echo "=== Extract haplotype EDTA annotations ==="
echo "Species       : ${SPECIES_SET}"
echo "PZI EDTA      : ${PZI_SOURCE_DIR}"
echo "PSP EDTA      : ${PSP_SOURCE_DIR}"
echo "PPA EDTA      : ${PPA_SOURCE_DIR}"
echo "Input/output  : ${INPUTS}/{species}/{genome,edta}"
echo "Mode          : exact forward/reverse-complement MD5 verified chromosome subset; no EDTA rerun"
echo "Python        : $(command -v "${PYTHON_BIN}") (${PYTHON_VERSION}); Conda env: ${PYTHON_ENV}"

args=(
    "${SCRIPT}"
    --species "${SPECIES_LIST[@]}"
    --root "${ROOT_DIR}"
    --inputs "${INPUTS}"
    --pzi-source-dir "${PZI_SOURCE_DIR}"
    --psp-source-dir "${PSP_SOURCE_DIR}"
    --ppa-source-dir "${PPA_SOURCE_DIR}"
)
[[ "${FORCE}" == "1" ]] && args+=(--force)
[[ "${FORCE}" == "0" || "${FORCE}" == "1" ]] || die "FORCE must be 0 or 1"

"${PYTHON_BIN}" "${args[@]}"
echo "[OK] Haplotype EDTA extraction completed under ${INPUTS}/*/edta/"
