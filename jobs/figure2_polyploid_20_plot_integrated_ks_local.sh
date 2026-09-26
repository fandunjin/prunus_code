#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Combine the three completed Figure 2 Ks result classes on one shared Ks axis.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CONDA_BASE="${CONDA_BASE:-${HOME}/miniconda3}"
PLOT_ENV="${PLOT_ENV:-wgd2-reference-clean}"
PYTHON_BIN="${PYTHON_BIN:-${CONDA_BASE}/envs/${PLOT_ENV}/bin/python3}"
FORCE_REPLOT="${FORCE_REPLOT:-0}"

die() { echo "[ERROR] $*" >&2; exit 2; }
[[ -x "${PYTHON_BIN}" ]] || die "Missing plotting Python: ${PYTHON_BIN}"
"${PYTHON_BIN}" -c 'import matplotlib, numpy' || die "Plotting Python requires matplotlib and numpy"

SELF_KS="${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/09_within_genome_ks/06_combined_source_data/within_genome_ks.tsv"
GTR_ORTHOLOGUE_KS="${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/06_panel_a_reference_exact/08_panel_a_reference_data/panel_a_lower_ks.tsv"
OUTDIR="${ROOT_DIR}/downstream_analysis/figure2_polyploid/fig2_integrated_ks"
PLOTTER="${ROOT_DIR}/scripts/plot_fig2_integrated_ks.py"

for file in "${SELF_KS}" "${GTR_ORTHOLOGUE_KS}" "${PLOTTER}"; do
  [[ -s "${file}" ]] || die "Missing input: ${file}"
done
if [[ -e "${OUTDIR}/COMPLETE" && "${FORCE_REPLOT}" != "1" ]]; then
  die "Output is already complete: ${OUTDIR}; rerun with FORCE_REPLOT=1"
fi

mkdir -p "${OUTDIR}"
PLOT_ARGS=(
  --self-ks "${SELF_KS}"
  --gtr-orthologue-ks "${GTR_ORTHOLOGUE_KS}"
  --outdir "${OUTDIR}"
  --bandwidth 0.035
  --orthologue-bandwidth 0.04
  --xmax 3.0
  --grid-points 3001
)
if [[ "${FORCE_REPLOT}" == "1" ]]; then
  PLOT_ARGS+=(--overwrite)
fi
MPLBACKEND=Agg QT_QPA_PLATFORM=offscreen "${PYTHON_BIN}" "${PLOTTER}" "${PLOT_ARGS[@]}"

for extension in svg pdf png tiff; do
  [[ -s "${OUTDIR}/figure2_integrated_ks.${extension}" ]] || die "Missing output: ${OUTDIR}/figure2_integrated_ks.${extension}"
done
echo "[OK] Figure 2 integrated Ks plot: ${OUTDIR}"
