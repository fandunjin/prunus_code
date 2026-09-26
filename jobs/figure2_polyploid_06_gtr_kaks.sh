#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}${CONDA_ROOT}/etc/profile.d/conda.sh}"
WGD_ENV="${WGD_ENV:-wgd_kaks}"

set +u
source "${CONDA_SH}"
conda activate "${WGD_ENV}"
set -u

exec bash "${ROOT_DIR}/jobs/figure2_polyploid_06_gtr_kaks_local.sh"
