#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/downstream_analysis/karyotype/run_20260803_205037}"
PAIR_DIR="${RUN_DIR}/pairs/psp_vs_ppa"
FIG_DIR="${RUN_DIR}/figures"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
JCVI_ENV="${JCVI_ENV:-jcvi}"
THREADS="${THREADS:-32}"
CSCORE="${CSCORE:-0.70}"
MIN_ANCHORS="${MIN_ANCHORS:-4}"
MINSPAN="${MINSPAN:-30}"
FORCE="${FORCE:-0}"

die() {
    echo "[ERROR] $*" >&2
    exit 1
}

require_file() {
    [[ -s "$1" ]] || die "Missing or empty file: $1"
}

mkdir -p "$PAIR_DIR" "$FIG_DIR"
export MPLBACKEND="${MPLBACKEND:-Agg}"
export MPLCONFIGDIR="${RUN_DIR}/mplconfig_psp_ppa"
mkdir -p "$MPLCONFIGDIR"

[[ -s "$CONDA_SH" ]] || die "Conda initialization script not found: $CONDA_SH"
source "$CONDA_SH"
conda activate "$JCVI_ENV"

command -v python >/dev/null 2>&1 || die "python is unavailable in $JCVI_ENV"
command -v diamond >/dev/null 2>&1 || die "diamond is unavailable in $JCVI_ENV"
python -c 'import jcvi' >/dev/null 2>&1 || die "jcvi is not importable in $JCVI_ENV"

INPUT_ROOT="${RUN_DIR}/inputs"
for f in \
    "$INPUT_ROOT/psp/psp.bed" \
    "$INPUT_ROOT/psp/psp.pep" \
    "$INPUT_ROOT/ppa/ppa.bed" \
    "$INPUT_ROOT/ppa/ppa.pep"; do
    require_file "$f"
done

ln -sfn "$INPUT_ROOT/psp/psp.bed" "$PAIR_DIR/psp.bed"
ln -sfn "$INPUT_ROOT/psp/psp.pep" "$PAIR_DIR/psp.pep"
ln -sfn "$INPUT_ROOT/ppa/ppa.bed" "$PAIR_DIR/ppa.bed"
ln -sfn "$INPUT_ROOT/ppa/ppa.pep" "$PAIR_DIR/ppa.pep"

cd "$PAIR_DIR"
RAW_ANCHORS="psp.ppa.anchors"
SCREENED_ANCHORS="psp.ppa.anchors.new"
SIMPLE_ANCHORS="psp.ppa.anchors.simple"
DOTPLOT="$FIG_DIR/psp_vs_ppa.dotplot.png"

echo "=== JCVI psp vs ppa dotplot ==="
echo "Run directory : $RUN_DIR"
echo "Pair directory: $PAIR_DIR"
echo "JCVI env      : $JCVI_ENV"
echo "DIAMOND       : $(command -v diamond)"
echo "Force         : $FORCE"

if [[ "$FORCE" == "1" || ! -s "$RAW_ANCHORS" ]]; then
    python -m jcvi.compara.catalog ortholog psp ppa \
        --dbtype=prot \
        --align_soft=diamond_blastp \
        --cscore="$CSCORE" \
        --min_size="$MIN_ANCHORS" \
        --cpus="$THREADS" \
        --no_strip_names \
        --no_dotplot
else
    echo "[INFO] Existing anchors found; skip catalog: $RAW_ANCHORS"
fi
require_file "$RAW_ANCHORS"

if [[ "$FORCE" == "1" || ! -s "$SCREENED_ANCHORS" || ! -s "$SIMPLE_ANCHORS" ]]; then
    python -m jcvi.compara.synteny screen \
        --qbed=psp.bed \
        --sbed=ppa.bed \
        --minspan="$MINSPAN" \
        --simple \
        "$RAW_ANCHORS" "$SCREENED_ANCHORS"
else
    echo "[INFO] Existing screened anchors found; skip screen"
fi
require_file "$SCREENED_ANCHORS"

if [[ ! -s "$SIMPLE_ANCHORS" ]]; then
    if [[ -s "$SCREENED_ANCHORS.simple" ]]; then
        ln -sfn "$SCREENED_ANCHORS.simple" "$SIMPLE_ANCHORS"
    else
        die "Six-column simple anchors not found: $SIMPLE_ANCHORS"
    fi
fi

if [[ "$FORCE" == "1" || ! -s "$DOTPLOT" ]]; then
    python -m jcvi.graphics.dotplot "$RAW_ANCHORS" \
        --qbed=psp.bed \
        --sbed=ppa.bed \
        --figsize=12x12 \
        --dpi=300 \
        --format=png \
        --notex \
        --outfile="$DOTPLOT"
else
    echo "[INFO] Existing dotplot found; skip: $DOTPLOT"
fi
require_file "$DOTPLOT"

echo "[OK] psp vs ppa dotplot completed: $DOTPLOT"
