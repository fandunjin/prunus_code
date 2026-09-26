#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

# Submit one group per job, for example:
# GROUP=external THREADS=12 jsub < jobs/figure2_repartition_07_codeml_run.sh
set -euo pipefail
ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
exec bash "${ROOT_DIR}/jobs/figure2_repartition_07_codeml_run_local.sh"
