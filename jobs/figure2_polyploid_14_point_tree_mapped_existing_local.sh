#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Render a reference-style, tree-mapped POInT Fig. 2d from completed outputs.
# No likelihood calculation is rerun by this job.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
POINT_ROOT="${POINT_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/09_point_conditional_fixed_r1}"
PARTITION_ROOT="${PARTITION_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/07_point_partition_models_r1}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/11_point_panel_d_tree_mapped}"
PLOT_PYTHON_BIN="${PLOT_PYTHON_BIN:-}"
PLOT="${ROOT_DIR}/scripts/plot_fig2_panel_d_point_tree_mapped.py"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }

require_file "${POINT_ROOT}/04_summary/point_model_fit_all.tsv"
require_file "${POINT_ROOT}/04_summary/point_branch_transition_summary.tsv"
require_file "${POINT_ROOT}/03_model_runs/reference_WGD_bias_fix/point_species_tree.nex.out"
require_file "${POINT_ROOT}/03_model_runs/reference_WGD_bias_fix/normalized/point_posterior_summary.tsv"
require_file "${PARTITION_ROOT}/03_summary/five_model_comparison.tsv"
require_file "${PLOT}"
PILLAR_COUNT=$(awk 'NR > 1 && NF {n++} END {print n + 0}' "${POINT_ROOT}/03_model_runs/reference_WGD_bias_fix/normalized/point_posterior_summary.tsv")

if [[ -z "${PLOT_PYTHON_BIN}" ]]; then
    for candidate in \
        "${CONDA_ROOT}/envs/wgd2-reference-clean/bin/python" \
        "/opt/plant-wakr-venv/bin/python3" \
        "python3"; do
        if { [[ "${candidate}" == */* && -x "${candidate}" ]] || command -v "${candidate}" >/dev/null 2>&1; } \
            && "${candidate}" -c 'import matplotlib, numpy' >/dev/null 2>&1; then
            PLOT_PYTHON_BIN="${candidate}"
            break
        fi
    done
fi
[[ -n "${PLOT_PYTHON_BIN}" ]] || die "No Python interpreter with matplotlib and numpy found; set PLOT_PYTHON_BIN"
[[ ! -e "${OUT_ROOT}" ]] || die "Output directory already exists: ${OUT_ROOT}; set a new OUT_ROOT for a rerun"

"${PLOT_PYTHON_BIN}" "${PLOT}" \
    --fit "${POINT_ROOT}/04_summary/point_model_fit_all.tsv" \
    --transition "${POINT_ROOT}/04_summary/point_branch_transition_summary.tsv" \
    --tree "${POINT_ROOT}/03_model_runs/reference_WGD_bias_fix/point_species_tree.nex.out" \
    --model "${POINT_ROOT}/03_model_runs/reference_WGD_bias_fix/optimized_model.txt" \
    --partitions "${PARTITION_ROOT}/03_summary/five_model_comparison.tsv" \
    --outdir "${OUT_ROOT}"

cat > "${OUT_ROOT}/software_and_parameters.txt" <<EOF
method=official POInT v1.63; native conditional transition output
visual_layout=reference-adapted two-copy tree map
topology=PPA,(PSP,PZI)
initial_copy_count_per_track=${PILLAR_COUNT} optimized pillars
copy1_loss_transition=Dupl-->Copy2
copy2_loss_transition=Dupl-->Copy1
model_fit=${POINT_ROOT}/04_summary/point_model_fit_all.tsv
conditional_transitions=${POINT_ROOT}/04_summary/point_branch_transition_summary.tsv
optimized_tree=${POINT_ROOT}/03_model_runs/reference_WGD_bias_fix/point_species_tree.nex.out
partition_models=${PARTITION_ROOT}/03_summary/five_model_comparison.tsv
EOF
echo "[OK] Tree-mapped Figure 2d completed: ${OUT_ROOT}"
