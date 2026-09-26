#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Local audit plot: all nine GTR chromosomes from the formal Fig. 2c anchors.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
if [[ -z "${PYTHON_BIN:-}" ]]; then
    if [[ -x /opt/plant-wakr-venv/bin/python3 ]]; then
        PYTHON_BIN="/opt/plant-wakr-venv/bin/python3"
    else
        PYTHON_BIN="python3"
    fi
fi
FORMAL_ROOT="${FORMAL_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/04_fractbias_8track_reference_exact}"
OUT_ROOT="${OUT_ROOT:-${FORMAL_ROOT}/07_all_gtr_chromosome_audit_v2}"
PLOT_DIR="${PLOT_DIR:-${OUT_ROOT}/03_plot_direct_labels}"
ANCHORS="${FORMAL_ROOT}/04_quota_align/normalized/quota_filtered_anchors.tsv"
CATALOG="${FORMAL_ROOT}/01_catalog"
CALCULATOR="${ROOT_DIR}/scripts/fig2_fractbias_all_gtr_input.py"
PLOTTER="${ROOT_DIR}/scripts/plot_fig2_panel_c_fractbias_all_gtr.py"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }

for file in "$ANCHORS" "$CATALOG/gtr_gene_order.tsv" "$CALCULATOR" "$PLOTTER"; do require_file "$file"; done
command -v "$PYTHON_BIN" >/dev/null 2>&1 || die "Python unavailable: $PYTHON_BIN"

if [[ ! -s "$OUT_ROOT/01_windows/COMPLETE" ]]; then
    [[ ! -e "$OUT_ROOT/01_windows" ]] || die "Incomplete window directory already exists: $OUT_ROOT/01_windows"
    mkdir -p "$OUT_ROOT"
    "$PYTHON_BIN" "$CALCULATOR" --anchors "$ANCHORS" --catalog "$CATALOG" --outdir "$OUT_ROOT/01_windows"
else
    echo "[OK] Reusing completed all-GTR windows: $OUT_ROOT/01_windows"
fi
[[ ! -e "$PLOT_DIR" ]] || die "Plot directory already exists: $PLOT_DIR"
"$PYTHON_BIN" "$PLOTTER" \
    --fractbias "$OUT_ROOT/01_windows/all_gtr_fractbias_windows.tsv" \
    --selection "$OUT_ROOT/01_windows/all_gtr_selected_tracks.tsv" \
    --outdir "$PLOT_DIR"

printf '[OK] All-nine-GTR audit completed: %s\n' "$OUT_ROOT"
