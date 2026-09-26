#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# End-to-end WSL driver. Run from the repository root or from any directory.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
THREADS="${THREADS:-$(nproc)}"
FORCE="${FORCE:-1}"

cd "${ROOT_DIR}"
bash jobs/haploid_repartition_00_wsl_setup.sh

# Activate the environment in this shell before invoking the child scripts.
ENV_NAME="${ENV_NAME:-prunus-circos}"
if command -v conda >/dev/null 2>&1; then
    CONDA_BASE="$(conda info --base)"
    # shellcheck disable=SC1091
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate "${ENV_NAME}"
elif command -v micromamba >/dev/null 2>&1; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate "${ENV_NAME}"
else
    echo "[ERROR] conda/micromamba is not available after environment setup" >&2
    exit 2
fi

THREADS="${THREADS}" FORCE="${FORCE}" bash jobs/haploid_repartition_06_jcvi_wsl.sh
FORCE="${FORCE}" bash jobs/haploid_repartition_07_circos_wsl.sh
echo "[OK] WSL JCVI + native Circos workflow completed"
