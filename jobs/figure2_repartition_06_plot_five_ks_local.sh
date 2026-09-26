#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Legacy-only MAFFT/YN visualization. It is retained to preserve the earlier
# exploratory output, but it is not part of the formal Figure 2 workflow.
#
# Formal unified CodeML figure: run jobs/figure2_repartition_08_combine_plot_codeml_local.sh

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/downstream_analysis/figure2_repartition_ks}"
OUT="${OUT:-${BASE}/05_five_ks_figure}"
PLOT_PYTHON="${PLOT_PYTHON:-/opt/plant-wakr-venv/bin/python3}"
XMAX="${XMAX:-5.0}"
SELF_PLOT_XMAX="${SELF_PLOT_XMAX:-3.0}"
SELF_AXIS_XMAX="${SELF_AXIS_XMAX:-4.5}"
BANDWIDTH="${BANDWIDTH:-0.035}"
FORCE="${FORCE:-0}"

die() { echo "[ERROR] $*" >&2; exit 2; }
need() { [[ -s "$1" ]] || die "Missing input: $1"; }

need "${BASE}/03_yn_kaks/COMPLETE"
need "${ROOT_DIR}/scripts/plot_figure2_repartition_five_ks.py"
command -v "${PLOT_PYTHON}" >/dev/null 2>&1 || die "Python unavailable: ${PLOT_PYTHON}"
"${PLOT_PYTHON}" -c 'import matplotlib, numpy; print("matplotlib=" + matplotlib.__version__)' \
  || die "PLOT_PYTHON must provide matplotlib and numpy"

args=(--base "${BASE}" --outdir "${OUT}" --xmax "${XMAX}" --self-plot-xmax "${SELF_PLOT_XMAX}" --self-axis-xmax "${SELF_AXIS_XMAX}" --bandwidth "${BANDWIDTH}")
[[ "${FORCE}" == "1" ]] && args+=(--force)
"${PLOT_PYTHON}" "${ROOT_DIR}/scripts/plot_figure2_repartition_five_ks.py" "${args[@]}"

for file in pdf svg png tiff; do need "${OUT}/figures/figure2_repartition_five_ks.${file}"; done
echo "[OK] Corrected five Ks panels: ${OUT}/figures"
