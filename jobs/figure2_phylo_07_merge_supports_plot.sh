#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT="${PROJECT_ROOT}/downstream_analysis/phylogenetics"
OUT="${BASE}/08_final_summary"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
PLOT_ENV="${PLOT_ENV:-akr}"

test -e "${BASE}/07_topology_frequencies/COMPLETE" || {
    echo "[ERROR] Missing 07_topology_frequencies/COMPLETE" >&2
    exit 3
}
test ! -e "${OUT}" || {
    echo "[ERROR] Output exists: ${OUT}" >&2
    exit 3
}

set +u
source "${CONDA_SH}"
conda activate "${PLOT_ENV}"
set -u

for module in matplotlib; do
    python -c "import ${module}" >/dev/null 2>&1 || {
        echo "[ERROR] Python module unavailable in ${PLOT_ENV}: ${module}" >&2
        exit 2
    }
done

python "${ROOT}/scripts/fig2_phylo_merge_supports_plot.py" \
    --base "${BASE}" \
    --manifest "${BASE}/config/input_manifest.tsv" \
    --out-dir "${OUT}" \
    --datasets main sensitivity \
    --outgroup rch \
    --width 7.2 \
    --height 7.0

for dataset in main sensitivity; do
    test -s "${OUT}/${dataset}/final_astral_support_tree.pdf" || {
        echo "[ERROR] Missing PDF figure: ${dataset}" >&2
        exit 4
    }
    test -s "${OUT}/${dataset}/final_node_support.tsv" || exit 4
    test -s "${OUT}/${dataset}/final_astral_rch_rooted.annotated.tre" || exit 4
done

{
    echo "python=$(python --version 2>&1)"
    echo "matplotlib=$(python -c 'import matplotlib; print(matplotlib.__version__)')"
    echo "plot_env=${PLOT_ENV}"
    echo "outgroup=rch"
    echo "node_label_order=bootstrap/PP/gCF/likelihood-sCF"
} > "${OUT}/software_and_parameters.txt"

echo "[OK] Final phylogenomic support figures completed: ${OUT}"
