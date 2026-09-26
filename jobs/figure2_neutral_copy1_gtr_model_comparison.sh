#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# JSUB archive wrapper for the neutral copy1/GTR model comparison.
# The implementation is shared with the local WSL job so that parameters do
# not silently diverge between HPC and local runs.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}${CONDA_ROOT}/etc/profile.d/conda.sh}"
JCVI_ENV="${JCVI_ENV:-jcvi}"
THREADS="${THREADS:-32}"
PYTHON_BIN="${PYTHON_BIN:-python}"

source "${CONDA_SH}"
conda activate "${JCVI_ENV}"
ROOT_DIR="${ROOT_DIR}" THREADS="${THREADS}" PYTHON_BIN="${PYTHON_BIN}" \
    CONDA_SH="" JCVI_ENV="" \
    bash "${ROOT_DIR}/jobs/figure2_neutral_copy1_gtr_model_comparison_local.sh"
