#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail
ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
THREADS="${THREADS:-32}" exec bash "${ROOT_DIR}/jobs/figure2_repartition_04_iadhore_self_pairs_local.sh"
