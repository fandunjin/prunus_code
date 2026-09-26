#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/downstream_analysis/ks}"
PHYLO_ROOT="${PHYLO_ROOT:-${ROOT_DIR}/downstream_analysis/phylogenetics}"
ANCHOR_ROOT="${ANCHOR_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/reference_anchor_states}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/fig3abcd_like}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
PLOT_ENV="${PLOT_ENV:-akr}"
PLOT_SCRIPT="${ROOT_DIR}/scripts/plot_fig2_polyploid_fig3abcd.py"
KS_INPUT="${KS_INPUT:-${KS_ROOT}/kaks/collinear_all_pairs_YN/summary/all_merged_4dtv_kaks.tsv}"
TOPOLOGY_INPUT="${TOPOLOGY_INPUT:-${PHYLO_ROOT}/07_topology_frequencies/main/topology_frequencies.tsv}"
SUPPORT_INPUT="${SUPPORT_INPUT:-${PHYLO_ROOT}/08_final_summary/main/final_node_support.tsv}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }

[[ -e "${ANCHOR_ROOT}/COMPLETE" ]] || die "Missing completion marker: ${ANCHOR_ROOT}/COMPLETE"
[[ ! -e "${OUT_DIR}/COMPLETE" ]] || die "Output is already complete: ${OUT_DIR}"
require_file "${PLOT_SCRIPT}"
require_file "${KS_INPUT}"
require_file "${TOPOLOGY_INPUT}"
require_file "${SUPPORT_INPUT}"
require_file "${ANCHOR_ROOT}/summary/gtr_anchor_state_summary.tsv"
require_file "${ANCHOR_ROOT}/summary/gtr_joint_anchor_states.tsv"

mkdir -p "${OUT_DIR}" "${ROOT_DIR}/downstream_analysis/figure2_polyploid/logs"
set +u
source "${CONDA_SH}"
conda activate "${PLOT_ENV}"
set -u
for module in numpy matplotlib; do
    python -c "import ${module}" >/dev/null 2>&1 || die "Python module unavailable in ${PLOT_ENV}: ${module}"
done

python "${PLOT_SCRIPT}" \
    --ks "${KS_INPUT}" \
    --topologies "${TOPOLOGY_INPUT}" \
    --support "${SUPPORT_INPUT}" \
    --anchor-summary "${ANCHOR_ROOT}/summary/gtr_anchor_state_summary.tsv" \
    --joint-states "${ANCHOR_ROOT}/summary/gtr_joint_anchor_states.tsv" \
    --outdir "${OUT_DIR}"

for extension in svg pdf tiff; do
    require_file "${OUT_DIR}/figure2_polyploid_fig3abcd.${extension}"
done
for table in panel_a_ks_source_data.tsv panel_b_topology_source_data.tsv panel_c_anchor_source_data.tsv panel_d_joint_source_data.tsv figure_manifest.tsv; do
    require_file "${OUT_DIR}/${table}"
done

{
    echo "python=$(python --version 2>&1)"
    echo "numpy=$(python -c 'import numpy; print(numpy.__version__)')"
    echo "matplotlib=$(python -c 'import matplotlib; print(matplotlib.__version__)')"
    echo "ks_max=3.0"
    echo "kde_bandwidth=0.04"
    echo "primary_haplotypes=A1/B1"
    echo "sensitivity_haplotypes=A2/B2_in_anchor_tables"
    echo "anchor_boundary=not_anchored_is_not_a_gene_loss"
} > "${OUT_DIR}/software_and_parameters.txt"

echo "[OK] Figure 2 Fig. 3a-d-like bundle completed: ${OUT_DIR}"
