#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
set -euo pipefail
ROOT_DIR="${PROJECT_ROOT}"
cd "${ROOT_DIR}"
export PATH="${CONDA_ROOT}/bin:${PATH}"
exec bash "${ROOT_DIR}/jobs/repartition_phylo_post_astral_local.sh"
