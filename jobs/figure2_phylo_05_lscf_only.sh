#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT="${PROJECT_ROOT}/downstream_analysis/phylogenetics"
ALIGN="${BASE}/04_alignments"
OUT="${BASE}/06_species_trees"
CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
TREE_ENV="${TREE_ENV:-phylo}"
THREADS="${THREADS:-48}"
SEED="${SEED:-20260813}"

source "${CONDA_SH}"
conda activate "${TREE_ENV}"

IQTREE_BIN=""
for candidate in iqtree3 iqtree2 iqtree; do
    if command -v "${candidate}" >/dev/null 2>&1; then
        IQTREE_BIN="$(command -v "${candidate}")"
        break
    fi
done
[[ -n "${IQTREE_BIN}" ]] || { echo "[ERROR] IQ-TREE unavailable" >&2; exit 2; }

for dataset in main sensitivity; do
    DOUT="${OUT}/${dataset}"
    TREE="${DOUT}/astral/species_tree.bs10.clean.tre"
    GENE_TREES="${DOUT}/all_gene_trees.bs10.tre"
    ALIGNMENT="${ALIGN}/${dataset}/qc/concat_protein.fa"
    PREFIX="${DOUT}/concordance_lscf"

    for file in "${TREE}" "${GENE_TREES}" "${ALIGNMENT}"; do
        [[ -s "${file}" ]] || { echo "[ERROR] Missing input: ${file}" >&2; exit 3; }
    done

    if [[ ! -s "${PREFIX}.cf.tree" ]]; then
        "${IQTREE_BIN}" \
            -t "${TREE}" \
            -s "${ALIGNMENT}" \
            --scfl 1000 \
            -T "${THREADS}" \
            -seed "${SEED}" \
            --prefix "${PREFIX}"
    else
        echo "[SKIP] Existing likelihood-sCF output: ${PREFIX}.cf.tree"
    fi

    [[ -s "${PREFIX}.cf.tree" ]] || { echo "[ERROR] Missing ${PREFIX}.cf.tree" >&2; exit 4; }
    [[ -s "${PREFIX}.cf.stat" ]] || { echo "[ERROR] Missing ${PREFIX}.cf.stat" >&2; exit 4; }
done

{
    "${IQTREE_BIN}" --version
    echo "TREE_ENV=${TREE_ENV}"
    echo "THREADS=${THREADS}"
    echo "SCFL_REPLICATES=1000"
} > "${OUT}/lscf_software_versions.txt" 2>&1

touch "${OUT}/LSCF_COMPLETE"
echo "[OK] Likelihood-based sCF completed for main and sensitivity."
