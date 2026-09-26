#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
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
[[ ! -e "${OUT}" ]] || { echo "[ERROR] Output exists: ${OUT}" >&2; exit 2; }
[[ -s "${ASTRAL_JAR}" ]] || { echo "[ERROR] Set ASTRAL_JAR to ASTRAL 5.7.1 jar" >&2; exit 2; }
if [[ -z "${IQTREE_BIN}" ]]; then
  for candidate in iqtree3 iqtree2 iqtree; do command -v "${candidate}" >/dev/null 2>&1 && IQTREE_BIN="$(command -v "${candidate}")" && break; done
fi
[[ -n "${IQTREE_BIN}" ]] || { echo "[ERROR] IQ-TREE unavailable; set IQTREE_BIN" >&2; exit 2; }
command -v java >/dev/null 2>&1 || { echo "[ERROR] java unavailable" >&2; exit 2; }
mkdir -p "${OUT}"

for dataset in main sensitivity; do
  DOUT="${OUT}/${dataset}"; mkdir -p "${DOUT}/astral"
  find "${GENES}/${dataset}/genes" -maxdepth 1 -type f -name '*.treefile' -printf '%f\n' | sort | while read -r name; do cat "${GENES}/${dataset}/genes/${name}"; done > "${DOUT}/all_gene_trees.uncontracted.tre"
  for threshold in 10 30 50; do
    "${PYTHON_BIN}" "${ROOT_DIR}/scripts/fig2_phylo_contract_trees.py" --tree-dir "${GENES}/${dataset}/genes" --threshold "${threshold}" --out "${DOUT}/all_gene_trees.bs${threshold}.tre"
    java -Xmx"${ASTRAL_MEM}" -jar "${ASTRAL_JAR}" -i "${DOUT}/all_gene_trees.bs${threshold}.tre" -t 2 -o "${DOUT}/astral/species_tree.bs${threshold}.full_annotation.tre" 2> "${DOUT}/astral/species_tree.bs${threshold}.log"
  done
  java -Xmx"${ASTRAL_MEM}" -jar "${ASTRAL_JAR}" -i "${DOUT}/all_gene_trees.bs10.tre" -t 0 -o "${DOUT}/astral/species_tree.bs10.clean.tre" 2> "${DOUT}/astral/species_tree.bs10.clean.log"
  find "${GENES}/${dataset}/genes" -maxdepth 1 -type f -name '*.treefile' -printf '%f\n' | sort | sed 's/\.treefile$/.ufboot/' | awk -v dir="${GENES}/${dataset}/genes" '{print dir "/" $0}' > "${DOUT}/astral/bootstrap_file_list.txt"
  java -Xmx"${ASTRAL_MEM}" -jar "${ASTRAL_JAR}" -i "${DOUT}/all_gene_trees.bs10.tre" -b "${DOUT}/astral/bootstrap_file_list.txt" -g -r "${ASTRAL_BOOTSTRAPS}" -s "${SEED}" -o "${DOUT}/astral/gene_site_bootstrap.all_trees.tre" 2> "${DOUT}/astral/gene_site_bootstrap.log"
  tail -n 1 "${DOUT}/astral/gene_site_bootstrap.all_trees.tre" > "${DOUT}/astral/gene_site_bootstrap.final_tree.tre"
  cp "${DOUT}/astral/species_tree.bs10.full_annotation.tre" "${DOUT}/astral/species_tree.bs10.q123.tre"
  "${IQTREE_BIN}" -s "${ALIGN}/${dataset}/qc/concat_protein.fa" -p "${ALIGN}/${dataset}/qc/concat_partitions.nex" -m MFP+MERGE -B 1000 --alrt 1000 --bnni -T "${THREADS}" -seed "${SEED}" --prefix "${DOUT}/concat_partitioned"
  # gCF and likelihood-sCF are separate IQ-TREE analyses. IQ-TREE rejects
  # --gcf and --scfl in the same invocation.
  "${IQTREE_BIN}" -t "${DOUT}/astral/species_tree.bs10.clean.tre" --gcf "${DOUT}/all_gene_trees.bs10.tre" -T "${THREADS}" -seed "${SEED}" --prefix "${DOUT}/concordance"
  "${IQTREE_BIN}" -te "${DOUT}/astral/species_tree.bs10.clean.tre" -s "${ALIGN}/${dataset}/qc/concat_protein.fa" -m LG+F+G4 --scfl 1000 -T "${THREADS}" -seed "${SEED}" --prefix "${DOUT}/concordance_lscf"
done
{
  echo "iqtree=${IQTREE_BIN}"; "${IQTREE_BIN}" --version || true
  echo "astral_jar=${ASTRAL_JAR}"; sha256sum "${ASTRAL_JAR}"
  echo "astral_bootstraps=${ASTRAL_BOOTSTRAPS}"; echo "seed=${SEED}"
} > "${OUT}/software_versions.txt"
touch "${OUT}/COMPLETE"
echo "[OK] New repartitioned ASTRAL and concordance analyses completed: ${OUT}"
