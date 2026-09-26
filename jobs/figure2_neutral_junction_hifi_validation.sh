#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# JSUB archive wrapper for Figure 2 stage-5 HiFi junction review.
# HIFI_BAM and HIFI_BAI must be supplied explicitly in the submitted script or
# shell environment after verifying that the BAM is aligned to copy1.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}${CONDA_ROOT}/etc/profile.d/conda.sh}"
JCVI_ENV="${JCVI_ENV:-jcvi}"
PYTHON_BIN="${PYTHON_BIN:-python}"

source "${CONDA_SH}"
conda activate "${JCVI_ENV}"
ROOT_DIR="${ROOT_DIR}" PYTHON_BIN="${PYTHON_BIN}" \
    bash "${ROOT_DIR}/jobs/figure2_neutral_junction_hifi_validation_local.sh"
