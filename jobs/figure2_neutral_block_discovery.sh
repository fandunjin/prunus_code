#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# JSUB archive wrapper for Figure 2 stages 2--4.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}${CONDA_ROOT}/etc/profile.d/conda.sh}"
JCVI_ENV="${JCVI_ENV:-jcvi}"
THREADS="${THREADS:-32}"
PYTHON_BIN="${PYTHON_BIN:-python}"

source "${CONDA_SH}"
conda activate "${JCVI_ENV}"
ROOT_DIR="${ROOT_DIR}" THREADS="${THREADS}" PYTHON_BIN="${PYTHON_BIN}" \
    CONDA_SH="" JCVI_ENV="" \
    bash "${ROOT_DIR}/jobs/figure2_neutral_block_discovery_local.sh"
