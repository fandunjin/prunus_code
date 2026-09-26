#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/analysisphylogenetics}"
ASTRAL_JAR="${ASTRAL_JAR:-${ROOT_DIR}/software/astral_5.7.8/Astral/astral.5.7.8.jar}"
THREADS="${THREADS:-20}"

cd "${ROOT_DIR}"
mkdir -p "${BASE}/logs"
THREADS="${THREADS}" bash "${ROOT_DIR}/jobs/repartition_phylo_04_gene_trees_local.sh"
ASTRAL_JAR="${ASTRAL_JAR}" THREADS=4 bash "${ROOT_DIR}/jobs/repartition_phylo_05_astral_concordance_local.sh"
bash "${ROOT_DIR}/jobs/repartition_phylo_06_topology_local.sh"
bash "${ROOT_DIR}/jobs/repartition_phylo_07_support_plot_local.sh"
bash "${ROOT_DIR}/jobs/repartition_phylo_08_mcmctree_local.sh"
printf '%s\n' "Revised repartitioned phylogenetic pipeline completed." > "${BASE}/PIPELINE_COMPLETE"
