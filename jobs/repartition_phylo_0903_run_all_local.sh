#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Run the complete phylogenomic workflow against the 0903 haplotype inputs.
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/0903/analysisphylogenetics}"
MANIFEST="${MANIFEST:-${BASE}/config/input_manifest.tsv}"
ASTRAL_JAR="${ASTRAL_JAR:-${ROOT_DIR}/software/astral_5.7.8/Astral/astral.5.7.8.jar}"
THREADS="${THREADS:-20}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
IQTREE_BIN="${IQTREE_BIN:-}"
MCMCTREE_BIN="${MCMCTREE_BIN:-mcmctree}"

die() { echo "[ERROR] $*" >&2; exit 2; }
[[ -s "${MANIFEST}" ]] || die "Missing 0903 manifest: ${MANIFEST}"
[[ -s "${ASTRAL_JAR}" ]] || die "Missing ASTRAL jar: ${ASTRAL_JAR}"
[[ ! -e "${BASE}/PIPELINE_COMPLETE" ]] || die "0903 workflow is already complete: ${BASE}"

cd "${ROOT_DIR}"
mkdir -p "${BASE}/logs"

if [[ ! -e "${BASE}/01_frozen_inputs/COMPLETE" ]]; then
  BASE="${BASE}" MANIFEST="${MANIFEST}" SOURCE_DIRECTORY="${ROOT_DIR}/0903/haploid_repartition" bash "${ROOT_DIR}/jobs/repartition_phylo_00_freeze_inputs_local.sh"
fi
if [[ ! -e "${BASE}/02_orthofinder/COMPLETE" ]]; then
  BASE="${BASE}" bash "${ROOT_DIR}/jobs/repartition_phylo_01_orthofinder_local.sh"
fi
if [[ ! -e "${BASE}/03_single_copy/COMPLETE" ]]; then
  BASE="${BASE}" bash "${ROOT_DIR}/jobs/repartition_phylo_02_extract_single_copy_local.sh"
fi
if [[ ! -e "${BASE}/04_alignments/COMPLETE" ]]; then
  BASE="${BASE}" THREADS="${THREADS}" bash "${ROOT_DIR}/jobs/repartition_phylo_03_align_trim_local.sh"
fi
if [[ ! -e "${BASE}/05_gene_trees/COMPLETE" ]]; then
  BASE="${BASE}" THREADS="${THREADS}" IQTREE_BIN="${IQTREE_BIN}" bash "${ROOT_DIR}/jobs/repartition_phylo_04_gene_trees_local.sh"
fi
if [[ ! -e "${BASE}/06_species_trees/COMPLETE" ]]; then
  BASE="${BASE}" ASTRAL_JAR="${ASTRAL_JAR}" THREADS=4 bash "${ROOT_DIR}/jobs/repartition_phylo_05_astral_concordance_local.sh"
fi
if [[ ! -e "${BASE}/07_topology_frequencies/COMPLETE" ]]; then
  BASE="${BASE}" bash "${ROOT_DIR}/jobs/repartition_phylo_06_topology_local.sh"
fi
if [[ ! -e "${BASE}/08_final_summary/COMPLETE" ]]; then
  BASE="${BASE}" bash "${ROOT_DIR}/jobs/repartition_phylo_07_support_plot_local.sh"
fi
if [[ ! -e "${BASE}/MCMCtree_COMPLETE" ]]; then
  BASE="${BASE}" IQTREE_BIN="${IQTREE_BIN}" MCMCTREE_BIN="${MCMCTREE_BIN}" THREADS=8 bash "${ROOT_DIR}/jobs/repartition_phylo_08_mcmctree_local.sh"
fi
DATED_OUT="${BASE}/08_final_summary/main/dated_astral_tree"
if [[ ! -e "${DATED_OUT}/COMPLETE" ]]; then
  [[ -s "${BASE}/08_final_summary/main/final_astral_rch_rooted.annotated.tre" ]] || die "Missing new ASTRAL support tree for dated plotting"
  MCMCTREE_TREE="${BASE}/09_mcmctree_subgenome17_outgroup_fossil55_35/03_mcmctree_chain1/FigTree.tre"
  [[ -s "${MCMCTREE_TREE}" ]] || die "Missing new subgenome17 MCMCtree FigTree tree: ${MCMCTREE_TREE}"
  "${PYTHON_BIN}" "${ROOT_DIR}/scripts/plot_dated_astral_support_tree.py" \
    --astral "${BASE}/08_final_summary/main/final_astral_rch_rooted.annotated.tre" \
    --mcmctree "${MCMCTREE_TREE}" \
    --manifest "${MANIFEST}" \
    --outdir "${DATED_OUT}"
fi

printf '%s\n' "0903 phylogenomic workflow completed." > "${BASE}/PIPELINE_COMPLETE"
echo "[OK] 0903 phylogenomic workflow completed: ${BASE}"
