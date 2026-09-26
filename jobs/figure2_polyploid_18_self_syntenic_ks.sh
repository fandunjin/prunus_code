#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

# Archived HPC equivalent of figure2_polyploid_18_self_syntenic_ks_local.sh.
# Submit one dataset at a time, e.g.:
#   jsub -J fig2_selfks_pzi_a DATASET=PZI_A jobs/figure2_polyploid_18_self_syntenic_ks.sh

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}${CONDA_ROOT}}"
WGD_ENV="${WGD_ENV:-wgd2-reference-clean}"
TOOLS_ENV="${TOOLS_ENV:-wgd2-reference}"
THREADS="${THREADS:-24}"
DATASET="${DATASET:-}"

mkdir -p "${ROOT_DIR}/downstream_analysis/figure2_polyploid/logs"
ROOT_DIR="${ROOT_DIR}" CONDA_BASE="${CONDA_BASE}" WGD_ENV="${WGD_ENV}" \
TOOLS_ENV="${TOOLS_ENV}" THREADS="${THREADS}" DATASET="${DATASET}" \
  bash "${ROOT_DIR}/jobs/figure2_polyploid_18_self_syntenic_ks_local.sh"
