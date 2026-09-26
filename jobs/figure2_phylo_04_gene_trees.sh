#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail
ROOT="${PROJECT_ROOT}/downstream_analysis/phylogenetics"
ALIGN="${BASE}/04_alignments"
OUT="${BASE}/05_gene_trees"
CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
TREE_ENV="${TREE_ENV:-phylo}"
THREADS="${THREADS:-48}"
UFB_REPLICATES="${UFB_REPLICATES:-2000}"
ALRT_REPLICATES="${ALRT_REPLICATES:-1000}"
SEED="${SEED:-20260813}"

test -e "${ALIGN}/COMPLETE" || { echo "[ERROR] Run figure2_phylo_03 first." >&2; exit 3; }
test ! -e "${OUT}" || { echo "[ERROR] Output exists: ${OUT}" >&2; exit 3; }
mkdir -p "${OUT}"
set +u
source "${CONDA_SH}"
conda activate "${TREE_ENV}"
set -u

IQTREE_BIN=""
for candidate in iqtree3 iqtree2 iqtree; do
    if command -v "${candidate}" >/dev/null 2>&1; then IQTREE_BIN="$(command -v "${candidate}")"; break; fi
done
test -n "${IQTREE_BIN}" || { echo "[ERROR] IQ-TREE unavailable" >&2; exit 2; }
command -v parallel >/dev/null || { echo "[ERROR] GNU parallel unavailable" >&2; exit 2; }
export IQTREE_BIN UFB_REPLICATES ALRT_REPLICATES SEED

for dataset in main sensitivity; do
    mkdir -p "${OUT}/${dataset}/genes"
    find "${ALIGN}/${dataset}/qc/retained" -maxdepth 1 -type f -name '*.fa' -print0 | sort -z | \
        parallel -0 -j "${THREADS}" --halt now,fail=1 \
        'base=$(basename {} .fa); "${IQTREE_BIN}" -s {} -m MFP -B "${UFB_REPLICATES}" --alrt "${ALRT_REPLICATES}" --bnni -T 1 -seed "${SEED}" --prefix "'"${OUT}/${dataset}/genes"'/${base}"'

    expected="$(wc -l < "${ALIGN}/${dataset}/qc/retained_orthogroups.txt")"
    observed="$(find "${OUT}/${dataset}/genes" -maxdepth 1 -type f -name '*.treefile' -size +0c | wc -l)"
    bootstrapped="$(find "${OUT}/${dataset}/genes" -maxdepth 1 -type f -name '*.ufboot' -size +0c | wc -l)"
    if [ "${observed}" -ne "${expected}" ] || [ "${bootstrapped}" -ne "${expected}" ]; then
        echo "[ERROR] ${dataset}: expected ${expected} trees/ufboot files, observed ${observed}/${bootstrapped}" >&2
        exit 3
    fi
done

"${IQTREE_BIN}" --version > "${OUT}/software_versions.txt" 2>&1
touch "${OUT}/COMPLETE"
