#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# jsub archive for rendering the native OrthoFinder duplication-mapping panel.
# Submit from the repository root: jsub < jobs/figure2_polyploid_15_plot_orthofinder_duplications.sh

set -euo pipefail
ROOT_DIR="${ROOT_DIR:-$PWD}"
cd "${ROOT_DIR}"
bash jobs/figure2_polyploid_15_plot_orthofinder_duplications_local.sh
