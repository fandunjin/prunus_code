#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

# HPC archive of the tested local runner. Override executable/environment paths
# only when the cluster installation differs from the local conda layout.
set -euo pipefail
ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/downstream_analysis/figure2_polyploid/logs"
ROOT_DIR="${ROOT_DIR}" bash "${ROOT_DIR}/jobs/figure2_polyploid_07_fractbias_reference_local.sh"
