#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Run the new haploid-repartition core phylogeny through support plotting.
set -euo pipefail
ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
for stage in \
  repartition_phylo_00_freeze_inputs_local.sh \
  repartition_phylo_01_orthofinder_local.sh \
  repartition_phylo_02_extract_single_copy_local.sh \
  repartition_phylo_03_align_trim_local.sh \
  repartition_phylo_04_gene_trees_local.sh; do
  bash "${ROOT_DIR}/jobs/${stage}"
done
echo "[NEXT] Set ASTRAL_JAR and run repartition_phylo_05_astral_concordance_local.sh"
