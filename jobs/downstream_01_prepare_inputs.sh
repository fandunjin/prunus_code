#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/analysis}"
HAP_RENAME_DIR="${HAP_RENAME_DIR:-${ROOT_DIR}/haploid_rename}"
PACKAGE_CHECKSUMS="${HAP_RENAME_DIR}/checksums.md5"
SAVED_PACKAGE_CHECKSUMS="${ANALYSIS_DIR}/.haploid_rename.checksums.md5"
MANIFEST="${MANIFEST:-${ROOT_DIR}/downstream_analysis/config/species_manifest.tsv}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
PY_ENV="${PY_ENV:-haphic}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
HAP_LABELS="pzi_hapA1,pzi_hapA2,pzi_hapB1,pzi_hapB2"
HAP_REFRESH_ARGS=()

mkdir -p "${ANALYSIS_DIR}" "${ANALYSIS_DIR}/logs"

for label in pzi_hapA1 pzi_hapA2 pzi_hapB1 pzi_hapB2; do
    for input in \
        "${HAP_RENAME_DIR}/genome/${label}.fa" \
        "${HAP_RENAME_DIR}/gff/${label}.gff3" \
        "${HAP_RENAME_DIR}/cds/${label}.cds.fa" \
        "${HAP_RENAME_DIR}/pep/${label}.pep.fa"; do
        [[ -s "${input}" ]] || { echo "[ERROR] Missing renamed haplotype input: ${input}" >&2; exit 2; }
    done
done
[[ -s "${PACKAGE_CHECKSUMS}" ]] || { echo "[ERROR] Missing renamed-package checksums: ${PACKAGE_CHECKSUMS}" >&2; exit 2; }
if ! cmp -s "${PACKAGE_CHECKSUMS}" "${SAVED_PACKAGE_CHECKSUMS}"; then
    HAP_REFRESH_ARGS+=(--force)
    echo "[INFO] Renamed haplotype package changed; refreshing canonical inputs."
fi

set +u
source "${CONDA_SH}"
conda activate "${PY_ENV}"
set -u

# Refresh the four renamed haplotypes even when analysis/canonical already
# contains files or symlinks staged by an earlier chr01-style run.
"${PYTHON_BIN}" "${ROOT_DIR}/scripts/downstream_prepare_inputs.py" \
    --root "${ROOT_DIR}" \
    --manifest "${MANIFEST}" \
    --out-dir "${ANALYSIS_DIR}" \
    --labels "${HAP_LABELS}" \
    "${HAP_REFRESH_ARGS[@]}" \
    --skip-orthofinder-sets

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/downstream_prepare_inputs.py" \
    --root "${ROOT_DIR}" \
    --manifest "${MANIFEST}" \
    --out-dir "${ANALYSIS_DIR}"

cp -f "${PACKAGE_CHECKSUMS}" "${SAVED_PACKAGE_CHECKSUMS}"

echo "Renamed haplotype package: ${HAP_RENAME_DIR}"
echo "Canonical manifest:"
cat "${ANALYSIS_DIR}/canonical_manifest.tsv"
echo
echo "Main OrthoFinder proteins:"
cat "${ANALYSIS_DIR}/orthofinder/main/pep_manifest.tsv"
echo
echo "Sensitivity OrthoFinder proteins:"
cat "${ANALYSIS_DIR}/orthofinder/sensitivity/pep_manifest.tsv"
