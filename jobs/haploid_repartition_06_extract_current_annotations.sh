#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}${CONDA_ROOT}/etc/profile.d/conda.sh}"
PYTHON_ENV="${PYTHON_ENV:-haphic}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
FORCE="${FORCE:-0}"
SPECIES_SET="${SPECIES_SET:-pzi psp ppa}"
SCRIPT="${ROOT_DIR}/scripts/extract_repartitioned_haplotype_annotations.py"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing or empty $2: $1"; }

require_file "${SCRIPT}" "current annotation extraction script"
require_file "${CONDA_SH}" "Conda initialization script"
source "${CONDA_SH}"
command -v conda >/dev/null 2>&1 || die "Conda is unavailable after sourcing: ${CONDA_SH}"
conda activate "${PYTHON_ENV}"
read -r -a SPECIES_LIST <<< "${SPECIES_SET}"
[[ "${#SPECIES_LIST[@]}" -gt 0 ]] || die "SPECIES_SET is empty; choose one or more of pzi psp ppa"
for species in "${SPECIES_LIST[@]}"; do
  case "${species}" in
    pzi) require_file "${ROOT_DIR}/haploid_repartition/maps/pzi_repartition_map.tsv" "PZI repartition map" ;;
    psp) require_file "${ROOT_DIR}/haploid_repartition/maps/psp_repartition_map.tsv" "PSP repartition map" ;;
    ppa) require_file "${ROOT_DIR}/haploid_repartition/maps/ppa_repartition_map.tsv" "PPA repartition map" ;;
    *) die "Unknown species in SPECIES_SET: ${species} (choose pzi, psp or ppa)" ;;
  esac
done
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || die "Python is unavailable: ${PYTHON_BIN}"
PYTHON_VERSION="$("${PYTHON_BIN}" -c 'import sys; print(".".join(map(str, sys.version_info[:3])))')" \
  || die "Unable to query Python version: ${PYTHON_BIN}"
"${PYTHON_BIN}" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 6) else 1)' \
  || die "Python >= 3.6 is required; found ${PYTHON_VERSION}"

args=("${SCRIPT}" --project-root "${ROOT_DIR}")
args+=(--species "${SPECIES_LIST[@]}")
[[ "${FORCE}" == "1" ]] && args+=(--force)
[[ "${FORCE}" == "0" || "${FORCE}" == "1" ]] || die "FORCE must be 0 or 1"

echo "=== Extract current complete-assembly annotations ==="
echo "Species: ${SPECIES_SET}"
echo "PZI: annotation/genes/pzi_eviann/postprocess/pzi32.eviann.*"
echo "PSP: data/psp/annotation/genes/psp32_eviann/postprocess/psp32.eviann.*"
echo "PPA: data/ppa/annotation/eviann/postprocess/ppa.eviann.* (fallback: data/ppa/postprocess)"
echo "Output: ${ROOT_DIR}/haploid_repartition/inputs/{pzi,psp,ppa}/{gff,cds,pep}"
echo "Force overwrite: ${FORCE}"
echo "Python: $(command -v "${PYTHON_BIN}") (${PYTHON_VERSION}); Conda env: ${PYTHON_ENV}"

"${PYTHON_BIN}" "${args[@]}"
echo "[OK] Current annotation extraction completed."
