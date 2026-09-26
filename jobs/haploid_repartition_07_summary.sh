#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}${CONDA_ROOT}/etc/profile.d/conda.sh}"
PYTHON_ENV="${PYTHON_ENV:-haphic}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
SCRIPT="${SCRIPT:-${ROOT_DIR}/scripts/summarize_repartitioned_haplotypes.py}"
OUTPUT_DIR="${OUTPUT_DIR:-${ROOT_DIR}/annotation/repartition/summary}"
INPUT_ROOT="${INPUT_ROOT:-${ROOT_DIR}/haploid_repartition/inputs}"
GENE_ROOT="${GENE_ROOT:-${ROOT_DIR}/annotation/repartition/genes}"
LAI_ROOT="${LAI_ROOT:-${ROOT_DIR}/annotation/repartition/lai}"
NCRNA_ROOT="${NCRNA_ROOT:-${ROOT_DIR}/annotation/repartition/ncrna}"
FUNCTIONAL_ROOT="${FUNCTIONAL_ROOT:-${ROOT_DIR}/annotation/repartition/functional}"
SPECIES_SET="${SPECIES_SET:-pzi psp ppa}"
STRICT="${STRICT:-0}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing or empty $2: $1"; }

require_file "${CONDA_SH}" "Conda initialization script"
require_file "${SCRIPT}" "haplotype summary script"
[[ "${STRICT}" == "0" || "${STRICT}" == "1" ]] || die "STRICT must be 0 or 1"

source "${CONDA_SH}"
command -v conda >/dev/null 2>&1 || die "Conda is unavailable after sourcing ${CONDA_SH}"
conda activate "${PYTHON_ENV}"
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || die "Python is unavailable in ${PYTHON_ENV}: ${PYTHON_BIN}"
PYTHON_VERSION="$(${PYTHON_BIN} -c 'import sys; print(".".join(map(str, sys.version_info[:3])))')" \
    || die "Unable to query Python version: ${PYTHON_BIN}"
"${PYTHON_BIN}" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 6) else 1)' \
    || die "Python >= 3.6 is required; found ${PYTHON_VERSION}"

mkdir -p "${OUTPUT_DIR}" "${ROOT_DIR}/logs"
args=(
    "${SCRIPT}"
    --root "${ROOT_DIR}"
    --output-dir "${OUTPUT_DIR}"
    --input-root "${INPUT_ROOT}"
    --gene-root "${GENE_ROOT}"
    --lai-root "${LAI_ROOT}"
    --ncrna-root "${NCRNA_ROOT}"
    --functional-root "${FUNCTIONAL_ROOT}"
)
read -r -a SPECIES_LIST <<< "${SPECIES_SET//,/ }"
[[ "${#SPECIES_LIST[@]}" -gt 0 ]] || die "SPECIES_SET is empty"
args+=(--species "${SPECIES_LIST[@]}")
[[ "${STRICT}" == "1" ]] && args+=(--strict)

echo "=== Haplotype annotation summary ==="
echo "Species       : ${SPECIES_SET}"
echo "Input root    : ${INPUT_ROOT}"
echo "Gene root     : ${GENE_ROOT}"
echo "LAI root      : ${LAI_ROOT}"
echo "ncRNA root    : ${NCRNA_ROOT}"
echo "Functional root: ${FUNCTIONAL_ROOT}"
echo "Output        : ${OUTPUT_DIR}"
echo "Python        : $(command -v "${PYTHON_BIN}") (${PYTHON_VERSION}); Conda env: ${PYTHON_ENV}"
echo "Strict        : ${STRICT}"

"${PYTHON_BIN}" "${args[@]}"
echo "[OK] Haplotype summary completed: ${OUTPUT_DIR}/haplotype_annotation_summary.tsv"
