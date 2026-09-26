#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Corrected-haplotype Figure 2b: map predefined A/B homeolog LCAs to the
# collapsed biological-species backbone, then render an editable figure.
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
GENE_TREES="${GENE_TREES:-${ROOT_DIR}/analysisphylogenetics/05_gene_trees/main/genes}"
OUT="${OUT:-${ROOT_DIR}/analysisphylogenetics/09_duplication_node_mapping/main}"

die(){ echo "[ERROR] $*" >&2; exit 2; }
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || die "Python unavailable: ${PYTHON_BIN}"
"${PYTHON_BIN}" -c 'import matplotlib' >/dev/null 2>&1 || die "Python matplotlib unavailable: ${PYTHON_BIN}"
[[ -d "${GENE_TREES}" ]] || die "Missing corrected gene-tree directory: ${GENE_TREES}"
[[ -s "${ROOT_DIR}/scripts/fig2_map_corrected_homeolog_nodes.py" ]] || die "Missing mapping script"
[[ -s "${ROOT_DIR}/scripts/plot_fig2b_corrected_homeolog_node_mapping.py" ]] || die "Missing plotting script"

mkdir -p "${OUT}"
"${PYTHON_BIN}" "${ROOT_DIR}/scripts/fig2_map_corrected_homeolog_nodes.py" \
  --gene-trees "${GENE_TREES}" --outdir "${OUT}" --ufboot-cutoff 95
"${PYTHON_BIN}" "${ROOT_DIR}/scripts/plot_fig2b_corrected_homeolog_node_mapping.py" \
  --summary "${OUT}/homeolog_pair_node_summary.tsv" \
  --shared "${OUT}/shared_high_support_homeolog_mappings.tsv" \
  --outdir "${OUT}/figure2b"
echo "[OK] Corrected Figure 2b node mapping and visualization completed: ${OUT}/figure2b"
