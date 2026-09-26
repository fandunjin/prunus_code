#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Compatibility entry point for the native POInT likelihood workflow.

set -euo pipefail
ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
exec bash "${ROOT_DIR}/jobs/figure2_polyploid_08_point_official_local.sh" "$@"
