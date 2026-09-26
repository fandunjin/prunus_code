#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail
ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/downstream_analysis/figure2_polyploid/logs"
ROOT_DIR="${ROOT_DIR}" bash "${ROOT_DIR}/jobs/figure2_polyploid_09_panel_a_reference_preflight_local.sh"
