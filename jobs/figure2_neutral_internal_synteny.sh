#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# JSUB archive wrapper for neutral copy1/copy2 and independent haplotype synteny.
# Uses the local workflow unchanged so that HPC and WSL parameters cannot diverge.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}${CONDA_ROOT}/etc/profile.d/conda.sh}"
JCVI_ENV="${JCVI_ENV:-jcvi}"
THREADS="${THREADS:-32}"
PYTHON_BIN="${PYTHON_BIN:-python}"

source "${CONDA_SH}"
conda activate "${JCVI_ENV}"
ROOT_DIR="${ROOT_DIR}" THREADS="${THREADS}" PYTHON_BIN="${PYTHON_BIN}" \
    bash "${ROOT_DIR}/jobs/figure2_neutral_internal_synteny_local.sh"
