#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Full local workflow for the revised haplotype partition.
set -euo pipefail
ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/analysisphylogenetics}"

if [[ ! -e "${BASE}/05_gene_trees/COMPLETE" ]]; then
  bash "${ROOT_DIR}/jobs/repartition_phylo_run_core_local.sh"
fi
if [[ -z "${ASTRAL_JAR:-}" ]]; then
  echo "[NEXT] Core stages are complete. Set ASTRAL_JAR and rerun this command to continue from stage 05."
  exit 0
fi
if [[ ! -e "${BASE}/06_species_trees/COMPLETE" ]]; then
  BASE="${BASE}" ASTRAL_JAR="${ASTRAL_JAR}" bash "${ROOT_DIR}/jobs/repartition_phylo_05_astral_concordance_local.sh"
fi
if [[ ! -e "${BASE}/07_topology_frequencies/COMPLETE" ]]; then
  bash "${ROOT_DIR}/jobs/repartition_phylo_06_topology_local.sh"
fi
if [[ ! -e "${BASE}/08_final_summary/COMPLETE" ]]; then
  bash "${ROOT_DIR}/jobs/repartition_phylo_07_support_plot_local.sh"
fi
if [[ ! -e "${BASE}/MCMCtree_COMPLETE" ]]; then
  bash "${ROOT_DIR}/jobs/repartition_phylo_08_mcmctree_local.sh"
fi
echo "[OK] Revised haplotype phylogenetic workflow is complete: ${BASE}"
