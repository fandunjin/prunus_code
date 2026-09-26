#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Run the dependent analyses only after the active ASTRAL/concordance stage.
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/analysisphylogenetics}"
WAIT_SECONDS="${WAIT_SECONDS:-60}"
export PATH="${CONDA_ROOT}/bin:${PATH}"

while [[ ! -e "${BASE}/06_species_trees/COMPLETE" ]]; do
  sleep "${WAIT_SECONDS}"
done

if [[ ! -e "${BASE}/07_topology_frequencies/COMPLETE" ]]; then
  bash "${ROOT_DIR}/jobs/repartition_phylo_06_topology_local.sh"
fi
if [[ ! -e "${BASE}/08_final_summary/COMPLETE" ]]; then
  PYTHON_BIN="${PYTHON_BIN:-/usr/bin/python3}" \
    bash "${ROOT_DIR}/jobs/repartition_phylo_07_support_plot_local.sh"
fi
if [[ ! -e "${BASE}/MCMCtree_COMPLETE" ]]; then
  bash "${ROOT_DIR}/jobs/repartition_phylo_08_mcmctree_local.sh"
fi

echo "[OK] Repartitioned post-ASTRAL workflow completed: ${BASE}"
