#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# jsub archive; submit from the repository root.
set -euo pipefail
cd "${JSUB_SUBMISSION_DIR:-$PWD}"
bash jobs/figure2_corrected_haplotype_duplication_mapping_local.sh
