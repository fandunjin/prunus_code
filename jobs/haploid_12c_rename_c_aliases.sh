#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/haploid_rename}"
OUT_DIR="${OUT_DIR:-${HAP_RENAME_DIR}}"
SCRIPT="${ROOT_DIR}/scripts/rename_pzi_haplotype_package.py"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
ENV_NAME="${ENV_NAME:-circos}"
FORCE_RENAME="${FORCE_RENAME:-0}"

die() {
    echo "[ERROR] $*" >&2
    exit 1
}

[[ -s "${SCRIPT}" ]] || die "Missing rename script: ${SCRIPT}"
[[ -s "${CONDA_SH}" ]] || die "Missing conda setup: ${CONDA_SH}"
mkdir -p "${OUT_DIR}" "${ROOT_DIR}/logs"

# shellcheck disable=SC1090
source "${CONDA_SH}"
conda activate "${ENV_NAME}"
command -v python >/dev/null 2>&1 || die "python is unavailable in ${ENV_NAME}"

echo "=== PZI complete chromosome-renamed haplotype package ==="
echo "Root   : ${ROOT_DIR}"
echo "Output : ${OUT_DIR}"
echo "Env    : ${ENV_NAME}"
cmd=(python "${SCRIPT}" --root "${ROOT_DIR}" --out-dir "${OUT_DIR}")
if [[ "${FORCE_RENAME}" == "1" ]]; then
    cmd+=(--force)
fi
"${cmd[@]}"

echo
echo "=== Outputs ==="
find "${OUT_DIR}" -maxdepth 2 -type f -size +0 -printf '%p\n' | sort
