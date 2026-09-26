#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail
ROOT="${PROJECT_ROOT}/downstream_analysis/phylogenetics"
TREES="${BASE}/06_species_trees"
OUT="${BASE}/07_topology_frequencies"
PYTHON_BIN="${PYTHON_BIN:-python3}"
SEED="${SEED:-20260813}"

test -e "${TREES}/COMPLETE" || { echo "[ERROR] Run figure2_phylo_05 first." >&2; exit 3; }
test ! -e "${OUT}" || { echo "[ERROR] Output exists: ${OUT}" >&2; exit 3; }
mkdir -p "${OUT}"

for dataset in main sensitivity; do
    "${PYTHON_BIN}" "${ROOT}/scripts/fig2_phylo_topology_frequencies.py" \
        --trees "${TREES}/${dataset}/all_gene_trees.bs10.tre" \
        --quartets "${BASE}/config/focal_quartets.tsv" \
        --out-dir "${OUT}/${dataset}" \
        --replicates 1000 \
        --seed "${SEED}"
done
touch "${OUT}/COMPLETE"
