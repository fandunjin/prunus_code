#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Resume ASTRAL/concordance analyses without deleting completed outputs.
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/analysisphylogenetics}"
ALIGN="${BASE}/04_alignments"
GENES="${BASE}/05_gene_trees"
OUT="${OUT:-${BASE}/06_species_trees}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
IQTREE_BIN="${IQTREE_BIN:-}"
ASTRAL_JAR="${ASTRAL_JAR:-}"
ASTRAL_MEM="${ASTRAL_MEM:-16G}"
ASTRAL_BOOTSTRAPS="${ASTRAL_BOOTSTRAPS:-1000}"
THREADS="${THREADS:-4}"
SEED="${SEED:-20260825}"

[[ -e "${GENES}/COMPLETE" ]] || { echo "[ERROR] Run repartition_phylo_04_gene_trees_local.sh first" >&2; exit 2; }
[[ -s "${ASTRAL_JAR}" ]] || { echo "[ERROR] Set ASTRAL_JAR to an ASTRAL jar" >&2; exit 2; }
if [[ -z "${IQTREE_BIN}" ]]; then
  for candidate in iqtree3 iqtree2 iqtree; do
    if command -v "${candidate}" >/dev/null 2>&1; then IQTREE_BIN="$(command -v "${candidate}")"; break; fi
  done
fi
[[ -n "${IQTREE_BIN}" ]] || { echo "[ERROR] IQ-TREE unavailable; set IQTREE_BIN" >&2; exit 2; }
command -v java >/dev/null 2>&1 || { echo "[ERROR] java unavailable" >&2; exit 2; }
mkdir -p "${OUT}"

run_dataset() {
  local dataset="$1"
  local dout="${OUT}/${dataset}"
  mkdir -p "${dout}/astral"
  find "${GENES}/${dataset}/genes" -maxdepth 1 -type f -name '*.treefile' -printf '%f\n' | sort |
    while read -r name; do cat "${GENES}/${dataset}/genes/${name}"; done > "${dout}/all_gene_trees.uncontracted.tre"

  for threshold in 10 30 50; do
    if [[ ! -s "${dout}/all_gene_trees.bs${threshold}.tre" ]]; then
      "${PYTHON_BIN}" "${ROOT_DIR}/scripts/fig2_phylo_contract_trees.py" --tree-dir "${GENES}/${dataset}/genes" --threshold "${threshold}" --out "${dout}/all_gene_trees.bs${threshold}.tre"
    fi
    if [[ ! -s "${dout}/astral/species_tree.bs${threshold}.full_annotation.tre" ]]; then
      java -Xmx"${ASTRAL_MEM}" -jar "${ASTRAL_JAR}" -i "${dout}/all_gene_trees.bs${threshold}.tre" -t 2 -o "${dout}/astral/species_tree.bs${threshold}.full_annotation.tre" 2> "${dout}/astral/species_tree.bs${threshold}.log"
    fi
  done

  if [[ ! -s "${dout}/astral/species_tree.bs10.clean.tre" ]]; then
    java -Xmx"${ASTRAL_MEM}" -jar "${ASTRAL_JAR}" -i "${dout}/all_gene_trees.bs10.tre" -t 0 -o "${dout}/astral/species_tree.bs10.clean.tre" 2> "${dout}/astral/species_tree.bs10.clean.log"
  fi
  if [[ ! -s "${dout}/astral/gene_site_bootstrap.final_tree.tre" ]]; then
    find "${GENES}/${dataset}/genes" -maxdepth 1 -type f -name '*.treefile' -printf '%f\n' | sort | sed 's/\.treefile$/.ufboot/' | awk -v dir="${GENES}/${dataset}/genes" '{print dir "/" $0}' > "${dout}/astral/bootstrap_file_list.txt"
    if [[ ! -s "${dout}/astral/gene_site_bootstrap.all_trees.tre" ]]; then
      java -Xmx"${ASTRAL_MEM}" -jar "${ASTRAL_JAR}" -i "${dout}/all_gene_trees.bs10.tre" -b "${dout}/astral/bootstrap_file_list.txt" -g -r "${ASTRAL_BOOTSTRAPS}" -s "${SEED}" -o "${dout}/astral/gene_site_bootstrap.all_trees.tre" 2> "${dout}/astral/gene_site_bootstrap.log"
    fi
    tail -n 1 "${dout}/astral/gene_site_bootstrap.all_trees.tre" > "${dout}/astral/gene_site_bootstrap.final_tree.tre"
  fi
  cp "${dout}/astral/species_tree.bs10.full_annotation.tre" "${dout}/astral/species_tree.bs10.q123.tre"

  if [[ ! -s "${dout}/concat_partitioned.treefile" ]]; then
    "${IQTREE_BIN}" -s "${ALIGN}/${dataset}/qc/concat_protein.fa" -p "${ALIGN}/${dataset}/qc/concat_partitions.nex" -m MFP+MERGE -B 1000 --alrt 1000 --bnni -T "${THREADS}" -seed "${SEED}" --prefix "${dout}/concat_partitioned"
  fi
  if [[ ! -s "${dout}/concordance.cf.tree" ]]; then
    "${IQTREE_BIN}" -t "${dout}/astral/species_tree.bs10.clean.tre" --gcf "${dout}/all_gene_trees.bs10.tre" -T "${THREADS}" -seed "${SEED}" --prefix "${dout}/concordance"
  fi
  if [[ ! -s "${dout}/concordance_lscf.cf.tree" ]]; then
    "${IQTREE_BIN}" -te "${dout}/astral/species_tree.bs10.clean.tre" -s "${ALIGN}/${dataset}/qc/concat_protein.fa" -m LG+F+G4 --scfl 1000 -T "${THREADS}" -seed "${SEED}" --prefix "${dout}/concordance_lscf"
  fi
}

run_dataset main
run_dataset sensitivity
{
  "${IQTREE_BIN}" --version
  java -version
  echo "ASTRAL_JAR=${ASTRAL_JAR}"
  sha256sum "${ASTRAL_JAR}"
  echo "gCF and likelihood-sCF run as separate IQ-TREE invocations"
} > "${OUT}/software_versions.txt" 2>&1
touch "${OUT}/COMPLETE"
echo "[OK] ASTRAL and concordance analyses completed: ${OUT}"
