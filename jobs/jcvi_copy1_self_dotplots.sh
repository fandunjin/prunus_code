#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

# Self-synteny analyses only:
#   pzi_copy1 vs pzi_copy1
#   psp_copy1 vs psp_copy1
# The PSP input used in the existing JCVI run is relabelled psp_copy1 here.

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/downstream_analysis/karyotype/run_20260803_205037}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/downstream_analysis/karyotype/copy1_self_jcvi}"
PAIR_ROOT="${OUT_DIR}/pairs"
FIG_DIR="${OUT_DIR}/figures"
SUMMARY="${OUT_DIR}/self_synteny_summary.tsv"

CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
JCVI_ENV="${JCVI_ENV:-jcvi}"
THREADS="${THREADS:-32}"
CSCORE="${CSCORE:-0.70}"
MIN_ANCHORS="${MIN_ANCHORS:-4}"
MINSPAN="${MINSPAN:-30}"

die() {
    echo "[ERROR] $*" >&2
    exit 1
}

require_file() {
    [[ -s "$1" ]] || die "Missing or empty file: $1"
}

anchor_count() {
    awk 'NF >= 2 && $1 !~ /^#/ {n++} END {print n+0}' "$1"
}

plot_dotplot() {
    local anchors="$1"
    local bed="$2"
    local output="$3"
    local format="$4"

    python -m jcvi.graphics.dotplot "$anchors" \
        --qbed="$bed" \
        --sbed="$bed" \
        --figsize=12x12 \
        --dpi=300 \
        --format="$format" \
        --notex \
        --outfile="$output"
}

run_self_synteny() {
    local target_label="$1"
    local source_label="$2"
    local source_dir="${SOURCE_RUN}/inputs/${source_label}"
    local pair_dir="${PAIR_ROOT}/${target_label}_self"
    local source_bed="${source_dir}/${source_label}.bed"
    local source_pep="${source_dir}/${source_label}.pep"
    local bed="${target_label}.bed"
    local pep="${target_label}.pep"
    local raw="${target_label}.${target_label}.anchors"
    local screened="${target_label}.${target_label}.anchors.new"
    local simple="${target_label}.${target_label}.anchors.simple"
    local raw_png="${FIG_DIR}/${target_label}_self.raw.dotplot.png"
    local raw_pdf="${FIG_DIR}/${target_label}_self.raw.dotplot.pdf"
    local screened_png="${FIG_DIR}/${target_label}_self.screened.dotplot.png"
    local screened_pdf="${FIG_DIR}/${target_label}_self.screened.dotplot.pdf"

    require_file "$source_bed"
    require_file "$source_pep"
    mkdir -p "$pair_dir"
    ln -sfn "$source_bed" "${pair_dir}/${bed}"
    ln -sfn "$source_pep" "${pair_dir}/${pep}"

    cd "$pair_dir"
    echo
    echo "=== ${target_label} self-synteny ==="
    echo "BED       : $source_bed"
    echo "Protein   : $source_pep"
    echo "Pair dir  : $pair_dir"

    if [[ ! -s "$raw" ]]; then
        echo "[RUN] JCVI catalog ortholog ${target_label} ${target_label}"
        python -m jcvi.compara.catalog ortholog \
            "$target_label" "$target_label" \
            --dbtype=prot \
            --align_soft=diamond_blastp \
            --cscore="$CSCORE" \
            --min_size="$MIN_ANCHORS" \
            --cpus="$THREADS" \
            --no_strip_names \
            --no_dotplot
    else
        echo "[SKIP] Existing raw anchors: ${pair_dir}/${raw}"
    fi
    require_file "$raw"

    if [[ ! -s "$screened" ]]; then
        echo "[RUN] Screen self-synteny anchors"
        python -m jcvi.compara.synteny screen \
            --qbed="$bed" \
            --sbed="$bed" \
            --minspan="$MINSPAN" \
            --simple \
            "$raw" "$screened"
    else
        echo "[SKIP] Existing screened anchors: ${pair_dir}/${screened}"
    fi
    require_file "$screened"

    # Depending on JCVI version, --simple may name this file after either the
    # raw or screened anchors. Keep the conventional name for downstream use.
    if [[ ! -s "$simple" && -s "${screened}.simple" ]]; then
        ln -sfn "${screened}.simple" "$simple"
    fi

    if [[ ! -s "$raw_png" ]]; then
        plot_dotplot "$raw" "$bed" "$raw_png" png
    else
        echo "[SKIP] Existing figure: $raw_png"
    fi
    if [[ ! -s "$raw_pdf" ]]; then
        plot_dotplot "$raw" "$bed" "$raw_pdf" pdf
    else
        echo "[SKIP] Existing figure: $raw_pdf"
    fi
    if [[ ! -s "$screened_png" ]]; then
        plot_dotplot "$screened" "$bed" "$screened_png" png
    else
        echo "[SKIP] Existing figure: $screened_png"
    fi
    if [[ ! -s "$screened_pdf" ]]; then
        plot_dotplot "$screened" "$bed" "$screened_pdf" pdf
    else
        echo "[SKIP] Existing figure: $screened_pdf"
    fi

    require_file "$raw_png"
    require_file "$raw_pdf"
    require_file "$screened_png"
    require_file "$screened_pdf"

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$target_label" "$source_label" "$(anchor_count "$raw")" \
        "$(anchor_count "$screened")" "${pair_dir}/${raw}" \
        "${pair_dir}/${screened}" "$raw_png" "$screened_png" >> "$SUMMARY"
}

mkdir -p "$OUT_DIR" "$PAIR_ROOT" "$FIG_DIR"
require_file "$CONDA_SH"

# shellcheck disable=SC1090
source "$CONDA_SH"
conda activate "$JCVI_ENV"

export MPLBACKEND="${MPLBACKEND:-Agg}"
export MPLCONFIGDIR="${OUT_DIR}/mplconfig"
mkdir -p "$MPLCONFIGDIR"

command -v python >/dev/null 2>&1 || die "python is unavailable in ${JCVI_ENV}"
command -v diamond >/dev/null 2>&1 || die "diamond is unavailable in ${JCVI_ENV}"
python -c 'import jcvi' >/dev/null 2>&1 || die "jcvi is not importable in ${JCVI_ENV}"

echo "=== PZI/PSP copy1 self-synteny with JCVI ==="
echo "Source run : $SOURCE_RUN"
echo "Output     : $OUT_DIR"
echo "JCVI env   : $JCVI_ENV"
echo "DIAMOND    : $(command -v diamond)"
echo "Threads    : $THREADS"
echo "C-score    : $CSCORE"
echo "Min anchors: $MIN_ANCHORS"
echo "Min span   : $MINSPAN"

printf 'label\tsource_label\traw_anchor_pairs\tscreened_anchor_pairs\traw_anchors\tscreened_anchors\traw_dotplot_png\tscreened_dotplot_png\n' > "$SUMMARY"

run_self_synteny pzi_copy1 pzi_copy1
run_self_synteny psp_copy1 psp

echo
echo "[OK] Both self-synteny analyses completed."
echo "Summary : $SUMMARY"
echo "Figures : $FIG_DIR"

