#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Reference-method MAFFT-default + FastTree-default family trees for the
# i-ADHoRe-derived family set.

set -euo pipefail
ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping}"
INPUT="${INPUT:-${BASE}/00_family_selection_reference_iadhore_allastral}"
OUT="${OUT:-${BASE}/02_family_trees_reference_iadhore_allastral}"
export ROOT_DIR BASE INPUT OUT
exec bash "${ROOT_DIR}/jobs/figure2_polyploid_04_wgd_family_trees_local.sh"
