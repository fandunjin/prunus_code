#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# HPC archive for the local/Python all-PZI retention plotting workflow.
set -euo pipefail
cd "${JSUB_SUBMISSION_DIR:-$PWD}"
bash jobs/figure2_polyploid_07c_fractbias_all_gtr_all_pzi_local.sh
