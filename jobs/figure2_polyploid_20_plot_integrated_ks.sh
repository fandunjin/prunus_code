#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

# Archived HPC equivalent. The job only redraws existing Ks source tables.
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export FORCE_REPLOT="${FORCE_REPLOT:-1}"
bash "${ROOT_DIR}/jobs/figure2_polyploid_20_plot_integrated_ks_local.sh"
