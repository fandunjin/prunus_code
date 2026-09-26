#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/analysisphylogenetics}"
TREES="${BASE}/06_species_trees"
OUT="${OUT:-${BASE}/07_topology_frequencies}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
SEED="${SEED:-20260825}"

[[ -e "${TREES}/COMPLETE" ]] || { echo "[ERROR] Run repartition_phylo_05_astral_concordance_local.sh first" >&2; exit 2; }
[[ ! -e "${OUT}" ]] || { echo "[ERROR] Output exists: ${OUT}" >&2; exit 2; }
mkdir -p "${OUT}"
for dataset in main sensitivity; do
  "${PYTHON_BIN}" "${ROOT_DIR}/scripts/fig2_phylo_topology_frequencies.py" \
    --trees "${TREES}/${dataset}/all_gene_trees.bs10.tre" \
    --quartets "${BASE}/config/focal_quartets.tsv" \
    --out-dir "${OUT}/${dataset}" --replicates 1000 --seed "${SEED}"
done
touch "${OUT}/COMPLETE"
echo "[OK] New repartitioned topology frequencies completed: ${OUT}"
