#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/analysis}"
SCRIPTS_DIR="${SCRIPTS_DIR:-${ROOT_DIR}/scripts}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
SYNTENY_ENV="${SYNTENY_ENV:-synteny}"
JCVI_DIR="${JCVI_DIR:-${ANALYSIS_DIR}/karyotype/jcvi}"
FIG_DIR="${FIG_DIR:-${JCVI_DIR}/figures}"
CUSTOM_FIG_DIR="${CUSTOM_FIG_DIR:-${FIG_DIR}/custom_svg}"
NATIVE_FIG_DIR="${NATIVE_FIG_DIR:-${FIG_DIR}/jcvi_native}"
MIN_BLOCK_PAIRS="${MIN_BLOCK_PAIRS:-10}"
MAX_DOTPLOT_POINTS="${MAX_DOTPLOT_POINTS:-120000}"
RUN_JCVI_NATIVE="${RUN_JCVI_NATIVE:-0}"
# On the current HPC synteny environment, JCVI/Matplotlib vector backends can
# fail with "object __array__ method not producing an array" and leave broken
# SVG/PDF files. PNG uses the Agg raster backend and is the reliable default.
NATIVE_FORMATS="${NATIVE_FORMATS:-png}"

export MPLBACKEND="${MPLBACKEND:-Agg}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${FIG_DIR}/mplconfig}"

mkdir -p "${FIG_DIR}" "${CUSTOM_FIG_DIR}" "${NATIVE_FIG_DIR}" "${MPLCONFIGDIR}" "${ANALYSIS_DIR}/logs"

set +u
source "${CONDA_SH}"
conda activate "${SYNTENY_ENV}"
set -u

if ! command -v python >/dev/null 2>&1; then
    echo "[ERROR] python not found in SYNTENY_ENV=${SYNTENY_ENV}" >&2
    exit 2
fi

for script in downstream_plot_jcvi_synteny.py downstream_make_jcvi_native_plot_inputs.py; do
    if [ ! -s "${SCRIPTS_DIR}/${script}" ]; then
        echo "[ERROR] Missing plotting script: ${SCRIPTS_DIR}/${script}" >&2
        exit 2
    fi
done

for table in jcvi_pair_status.tsv jcvi_syntenic_pairs.tsv jcvi_block_summary.tsv; do
    if [ ! -s "${JCVI_DIR}/summary/${table}" ]; then
        echo "[ERROR] Missing 09 result table: ${JCVI_DIR}/summary/${table}" >&2
        echo "[ERROR] Run jobs/downstream_09_karyotype_jcvi.sh successfully before 09b." >&2
        exit 2
    fi
done

STATUS_TABLE="${JCVI_DIR}/summary/jcvi_pair_status.tsv"
NORMALIZED_STATUS="${JCVI_DIR}/summary/jcvi_pair_status.for_09b.tsv"
status_header=""
if [ -s "${STATUS_TABLE}" ]; then
    IFS= read -r status_header < "${STATUS_TABLE}"
fi
if [ "${status_header}" = "comparison	query_label	ref_label	comparison_type	status	n_blocks	n_syntenic_pairs	screened_anchors	jcvi_toaligns" ]; then
    STATUS_FOR_PLOTS="${STATUS_TABLE}"
else
    echo "[WARN] Normalizing old 09 status table without comparison_type; rerun 09 for haplotype-internal plots." >&2
    awk 'BEGIN{FS=OFS="\t"} NR==1{print "comparison","query_label","ref_label","comparison_type","status","n_blocks","n_syntenic_pairs","screened_anchors","jcvi_toaligns"; next} NF>=8{print $1,$2,$3,"external_reference",$4,$5,$6,$7,$8}' "${STATUS_TABLE}" > "${NORMALIZED_STATUS}"
    STATUS_FOR_PLOTS="${NORMALIZED_STATUS}"
fi

validate_pdf() {
    pdf_file="$1"
    if [ ! -s "${pdf_file}" ]; then
        echo "[ERROR] Expected PDF was not created: ${pdf_file}" >&2
        return 1
    fi
    if ! tail -c 2048 "${pdf_file}" | grep -a "%%EOF" >/dev/null 2>&1; then
        echo "[ERROR] Broken PDF without %%EOF marker: ${pdf_file}" >&2
        echo "[ERROR] JCVI/Matplotlib probably failed while saving. Check this job stderr." >&2
        return 1
    fi
}

validate_png() {
    png_file="$1"
    if [ ! -s "${png_file}" ]; then
        echo "[ERROR] Expected PNG was not created: ${png_file}" >&2
        return 1
    fi
    if ! head -c 8 "${png_file}" | od -An -tx1 | grep -qi "89 50 4e 47 0d 0a 1a 0a"; then
        echo "[ERROR] Broken PNG header: ${png_file}" >&2
        return 1
    fi
}

run_jcvi_dotplot() {
    anchor_file="$1"
    qbed_file="$2"
    sbed_file="$3"
    outfile="$4"
    format="$5"
    python -m jcvi.graphics.dotplot \
        "${anchor_file}" \
        --qbed="${qbed_file}" \
        --sbed="${sbed_file}" \
        --figsize=12x12 \
        --dpi=300 \
        --format="${format}" \
        --notex \
        --outfile="${outfile}"
    if [ "${format}" = "pdf" ]; then
        validate_pdf "${outfile}"
    elif [ "${format}" = "png" ]; then
        validate_png "${outfile}"
    fi
}

run_jcvi_karyotype() {
    seqids_file="$1"
    layout_file="$2"
    outfile="$3"
    format="$4"
    python -m jcvi.graphics.karyotype \
        "${seqids_file}" \
        "${layout_file}" \
        --format="${format}" \
        --notex \
        --outfile="${outfile}"
    if [ "${format}" = "pdf" ]; then
        validate_pdf "${outfile}"
    elif [ "${format}" = "png" ]; then
        validate_png "${outfile}"
    fi
}

python "${SCRIPTS_DIR}/downstream_plot_jcvi_synteny.py" \
    --jcvi-dir "${JCVI_DIR}" \
    --outdir "${CUSTOM_FIG_DIR}" \
    --min-block-pairs "${MIN_BLOCK_PAIRS}" \
    --max-dotplot-points "${MAX_DOTPLOT_POINTS}"

if [ "${RUN_JCVI_NATIVE}" != "1" ]; then
    echo "[INFO] Skipping JCVI native graphics because RUN_JCVI_NATIVE=${RUN_JCVI_NATIVE}."
    echo "[INFO] The current synteny environment fails inside jcvi.graphics savefig with: object __array__ method not producing an array."
    echo "[INFO] Use the JCVI-anchor SVG figures in ${CUSTOM_FIG_DIR}, or rerun with RUN_JCVI_NATIVE=1 after fixing the plotting environment."
else
while IFS=$'\t' read -r comparison query ref comparison_type status n_blocks n_pairs screened_anchors jcvi_toaligns; do
    if [ "${comparison}" = "comparison" ]; then
        continue
    fi
    if [ "${status}" != "ok" ]; then
        echo "[WARN] Skip native JCVI plot for ${comparison}: status=${status}" >&2
        continue
    fi

    pair_work="${JCVI_DIR}/pairs/${comparison}"
    if [ "${comparison_type}" = "internal_haplotype" ]; then
        native_work="${NATIVE_FIG_DIR}/haplotype_internal/${comparison}"
    else
        native_work="${NATIVE_FIG_DIR}/reference_external/${comparison}"
    fi
    mkdir -p "${native_work}"

    for required in "${pair_work}/${query}.bed" "${pair_work}/${query}.pep" "${pair_work}/${ref}.bed" "${pair_work}/${ref}.pep" "${pair_work}/${query}.${ref}.anchors" "${pair_work}/${query}.${ref}.anchors.simple"; do
        if [ ! -s "${required}" ]; then
            echo "[WARN] Skip native JCVI plot for ${comparison}: missing ${required}" >&2
            continue 2
        fi
    done

    ln -sf "${pair_work}/${query}.bed" "${native_work}/${query}.bed"
    ln -sf "${pair_work}/${query}.pep" "${native_work}/${query}.pep"
    ln -sf "${pair_work}/${ref}.bed" "${native_work}/${ref}.bed"
    ln -sf "${pair_work}/${ref}.pep" "${native_work}/${ref}.pep"
    ln -sf "${pair_work}/${query}.${ref}.anchors" "${native_work}/${query}.${ref}.anchors"
    ln -sf "${pair_work}/${query}.${ref}.anchors.simple" "${native_work}/${query}.${ref}.anchors.simple"

    python "${SCRIPTS_DIR}/downstream_make_jcvi_native_plot_inputs.py" \
        --query "${query}" \
        --ref "${ref}" \
        --query-seqids "${JCVI_DIR}/inputs/${query}/${query}.seqids" \
        --ref-seqids "${JCVI_DIR}/inputs/${ref}/${ref}.seqids" \
        --query-bed "${native_work}/${query}.bed" \
        --ref-bed "${native_work}/${ref}.bed" \
        --anchors-simple "${native_work}/${query}.${ref}.anchors.simple" \
        --out-seqids "${native_work}/${comparison}.seqids" \
        --out-layout "${native_work}/${comparison}.layout"

    cd "${native_work}"
    IFS=',' read -r -a native_formats <<< "${NATIVE_FORMATS}"
    for native_format in "${native_formats[@]}"; do
        run_jcvi_dotplot \
            "${query}.${ref}.anchors" \
            "${query}.bed" \
            "${ref}.bed" \
            "${comparison}.jcvi_dotplot.${native_format}" \
            "${native_format}"

        run_jcvi_karyotype \
            "${comparison}.seqids" \
            "${comparison}.layout" \
            "${comparison}.jcvi_karyotype.${native_format}" \
            "${native_format}"
    done

    cd "${ROOT_DIR}"
done < "${STATUS_FOR_PLOTS}"
fi

cat > "${FIG_DIR}/README.txt" <<'EOF'
Generated by downstream_09b_jcvi_plots.sh from downstream_09 JCVI outputs.

Figure sets:
  custom_svg/dotplots/*.svg
    Script-generated SVG dotplots for every successful comparison. X-axis is
    reference, Y-axis is PZI query. Points are screened syntenic gene pairs.

  custom_svg/ribbons/*.svg
    Script-generated macro-synteny ribbon plots for every successful comparison.
    Only blocks with at least MIN_BLOCK_PAIRS screened syntenic gene pairs are
    drawn by default.

  custom_svg/summary_heatmaps/*.svg
    Global heatmaps for syntenic gene-pair counts, block counts, and mean pairs
    per block across all PZI query/reference combinations.

Optional native JCVI graphics:
  Set RUN_JCVI_NATIVE=1 to run python -m jcvi.graphics.dotplot/karyotype. On the
  current synteny environment this can fail during Matplotlib savefig with:
  object __array__ method not producing an array. The default RUN_JCVI_NATIVE=0
  therefore keeps the job successful and relies on the JCVI-anchor SVG figures
  above.

  jcvi_native/reference_external/<comparison>/*.jcvi_dotplot.png
  jcvi_native/haplotype_internal/<comparison>/*.jcvi_dotplot.png
    Native JCVI dotplots produced by python -m jcvi.graphics.dotplot from the
    raw JCVI anchors.

  jcvi_native/reference_external/<comparison>/*.jcvi_karyotype.png
  jcvi_native/haplotype_internal/<comparison>/*.jcvi_karyotype.png
    Native JCVI karyotype/ribbon-style plots produced by
    python -m jcvi.graphics.karyotype from .anchors.simple, seqids, and layout.

  PNG is the default native JCVI format because it is reliable on this HPC
  Matplotlib/Numpy environment. SVG/PDF can be requested with NATIVE_FORMATS,
  but vector backends may fail in the current synteny environment. PDF files
  are checked for the final %%EOF marker so failed saves do not silently leave
  unreadable partial files.

Recommended manuscript triage:
  Main karyotype figures: prioritize pzi_A2B2 vs ppe/pdu/pmu_refseq/pcam/pspe.
  Subgenome comparison: prioritize pzi_hapA1/A2/B1/B2 vs ppe/pcam/pspe.
  Outgroup controls: use rch/fve mainly as supplementary checks.
EOF

echo "[INFO] JCVI figures finished: ${FIG_DIR}"
