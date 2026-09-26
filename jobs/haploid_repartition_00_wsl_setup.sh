#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Set up a WSL Conda environment containing the tools used by the revised
# JCVI/Circos workflow. This script does not download anything unless the
# environment is missing and conda is available.

set -euo pipefail

ENV_NAME="${ENV_NAME:-prunus-circos}"
CONDA_BIN="${CONDA_BIN:-$(command -v conda || true)}"
MAMBA_BIN="${MAMBA_BIN:-$(command -v micromamba || true)}"

die() { echo "[ERROR] $*" >&2; exit 2; }

if [[ -n "${CONDA_BIN}" ]]; then
    CONDA_BASE="$(conda info --base)"
    # shellcheck disable=SC1091
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    if ! conda env list | awk '{print $1}' | grep -Fxq "${ENV_NAME}"; then
        echo "[INFO] Creating Conda environment ${ENV_NAME}"
        conda create -y -n "${ENV_NAME}" -c conda-forge -c bioconda \
            python=3.11 jcvi diamond circos
    else
        echo "[INFO] Existing Conda environment: ${ENV_NAME}"
    fi
    conda activate "${ENV_NAME}"
elif [[ -n "${MAMBA_BIN}" ]]; then
    eval "$(micromamba shell hook --shell bash)"
    if ! micromamba env list | awk '{print $1}' | grep -Fxq "${ENV_NAME}"; then
        echo "[INFO] Creating micromamba environment ${ENV_NAME}"
        micromamba create -y -n "${ENV_NAME}" -c conda-forge -c bioconda \
            python=3.11 jcvi diamond circos
    else
        echo "[INFO] Existing micromamba environment: ${ENV_NAME}"
    fi
    micromamba activate "${ENV_NAME}"
else
    die "Neither conda nor micromamba is available. Install Miniforge in WSL, then rerun this script."
fi

for exe in python diamond circos perl; do
    command -v "${exe}" >/dev/null 2>&1 || die "${exe} is not available in ${ENV_NAME}"
done
python -c 'import jcvi; print("jcvi", jcvi.__version__ if hasattr(jcvi, "__version__") else "available")'
circos -version >/dev/null 2>&1 || die "Circos is installed but could not start"
printf '[OK] WSL environment ready: %s\n' "${ENV_NAME}"
