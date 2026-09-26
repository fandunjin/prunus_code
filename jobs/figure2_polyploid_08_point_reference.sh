#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

# HPC archive of the native POInT likelihood runner.
set -euo pipefail
ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/downstream_analysis/figure2_polyploid/logs"
ROOT_DIR="${ROOT_DIR}" bash "${ROOT_DIR}/jobs/figure2_polyploid_08_point_official_local.sh"
