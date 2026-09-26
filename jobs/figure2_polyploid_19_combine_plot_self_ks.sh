#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

# Archived HPC equivalent of figure2_polyploid_19_combine_plot_self_ks_local.sh.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}${CONDA_ROOT}}"
WGD_ENV="${WGD_ENV:-wgd2-reference-clean}"
mkdir -p "${ROOT_DIR}/downstream_analysis/figure2_polyploid/logs"
ROOT_DIR="${ROOT_DIR}" CONDA_BASE="${CONDA_BASE}" WGD_ENV="${WGD_ENV}" \
  bash "${ROOT_DIR}/jobs/figure2_polyploid_19_combine_plot_self_ks_local.sh"
