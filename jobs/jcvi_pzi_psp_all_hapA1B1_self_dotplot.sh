#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

# Full self-synteny for the 16 chromosomes in PZI/PSP copy1, while retaining
# the haplotype chromosome names. The JCVI inputs contain hapA1 + hapB1:
# C01_A1--C08_A1 and C01_B1--C08_B1. No chr01--chr16 relabelling is done.

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/haploid_rename}"
PSP_GENES_DIR="${PSP_GENES_DIR:-${ROOT_DIR}/data/psp/haploid/annotation/genes}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/downstream_analysis/karyotype/pzi_psp_all_self_jcvi}"
STAGE_DIR="${OUT_DIR}/stage"
INPUT_ROOT="${OUT_DIR}/inputs"
PAIR_ROOT="${OUT_DIR}/pairs"
FIG_DIR="${OUT_DIR}/figures"
SUMMARY_DIR="${OUT_DIR}/summary"

CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
JCVI_ENV="${JCVI_ENV:-jcvi}"
THREADS="${THREADS:-32}"
CSCORE="${CSCORE:-0.70}"
MIN_ANCHORS="${MIN_ANCHORS:-4}"
MINSPAN="${MINSPAN:-30}"
FORCE="${FORCE:-0}"

PREPARE_SCRIPT="${ROOT_DIR}/scripts/downstream_prepare_jcvi_inputs.py"

die() { echo "[ERROR] $*" >&2; exit 1; }
require_file() { [[ -s "$1" ]] || die "Missing or empty file: $1"; }
anchor_count() { awk 'NF >= 2 && $1 !~ /^#/ {n++} END {print n+0}' "$1"; }

prepare_species() {
    local species="$1" gff_a="$2" pep_a="$3" gff_b="$4" pep_b="$5"
    local canonical="${STAGE_DIR}/canonical/${species}"
    mkdir -p "$canonical"
    require_file "$gff_a"; require_file "$pep_a"; require_file "$gff_b"; require_file "$pep_b"

    if [[ "$FORCE" == "1" || ! -s "${canonical}/primary.gff3" ]]; then
        {
            printf '##gff-version 3\n'
            awk '!/^#/ && NF >= 9 {print}' "$gff_a"
            awk '!/^#/ && NF >= 9 {print}' "$gff_b"
        } > "${canonical}/primary.gff3.tmp"
        mv "${canonical}/primary.gff3.tmp" "${canonical}/primary.gff3"
    fi
    if [[ "$FORCE" == "1" || ! -s "${canonical}/primary.pep.fa" ]]; then
        cat "$pep_a" "$pep_b" > "${canonical}/primary.pep.fa.tmp"
        mv "${canonical}/primary.pep.fa.tmp" "${canonical}/primary.pep.fa"
    fi
    require_file "${canonical}/primary.gff3"
    require_file "${canonical}/primary.pep.fa"
}

run_self() {
    local species="$1"
    local label="${species}_hapA1B1"
    local input_dir="${INPUT_ROOT}/${species}"
    local pair_dir="${PAIR_ROOT}/${label}_self"
    local bed="${label}.bed" pep="${label}.pep"
    local raw="${label}.${label}.anchors"
    local screened="${label}.${label}.anchors.new"
    local simple="${label}.${label}.anchors.simple"
    local raw_png="${FIG_DIR}/${label}_self.raw.dotplot.png"
    local raw_pdf="${FIG_DIR}/${label}_self.raw.dotplot.pdf"
    local screened_png="${FIG_DIR}/${label}_self.screened.dotplot.png"
    local screened_pdf="${FIG_DIR}/${label}_self.screened.dotplot.pdf"

    mkdir -p "$input_dir" "$pair_dir"
    python "$PREPARE_SCRIPT" \
        --analysis-dir "$STAGE_DIR" --outdir "$INPUT_ROOT" \
        --summary "${SUMMARY_DIR}/${species}_input_summary.tsv" --labels "$species"
    require_file "${input_dir}/${species}.bed"
    require_file "${input_dir}/${species}.pep"
    ln -sfn "${input_dir}/${species}.bed" "${pair_dir}/${bed}"
    ln -sfn "${input_dir}/${species}.pep" "${pair_dir}/${pep}"

    if [[ "$FORCE" == "1" ]]; then
        rm -f "${pair_dir}/${raw}" "${pair_dir}/${raw}.simple" \
            "${pair_dir}/${screened}" "${pair_dir}/${screened}.simple" \
            "${pair_dir}/${simple}" "$raw_png" "$raw_pdf" "$screened_png" "$screened_pdf"
    fi

    cd "$pair_dir"
    if [[ ! -s "$raw" ]]; then
        echo "[RUN] ${species}: JCVI catalog ortholog ${label} ${label}"
        python -m jcvi.compara.catalog ortholog "$label" "$label" \
            --dbtype=prot --align_soft=diamond_blastp \
            --cscore="$CSCORE" --min_size="$MIN_ANCHORS" --cpus="$THREADS" \
            --no_strip_names --no_dotplot
    else
        echo "[SKIP] ${species}: existing raw anchors"
    fi
    require_file "$raw"
    if [[ ! -s "$screened" ]]; then
        echo "[RUN] ${species}: screen anchors"
        python -m jcvi.compara.synteny screen \
            --qbed="$bed" --sbed="$bed" --minspan="$MINSPAN" --simple "$raw" "$screened"
    else
        echo "[SKIP] ${species}: existing screened anchors"
    fi
    require_file "$screened"
    if [[ ! -s "$simple" ]]; then
        for candidate in "${screened}.simple" "${raw}.simple"; do
            if [[ -s "$candidate" ]]; then ln -sfn "$candidate" "$simple"; break; fi
        done
    fi

    for spec in "${raw}|${raw_png}|png" "${raw}|${raw_pdf}|pdf" \
                "${screened}|${screened_png}|png" "${screened}|${screened_pdf}|pdf"; do
        IFS='|' read -r anchors output format <<< "$spec"
        if [[ ! -s "$output" ]]; then
            python -m jcvi.graphics.dotplot "$anchors" \
                --qbed="$bed" --sbed="$bed" --figsize=12x12 --dpi=300 \
                --format="$format" --nochpf --notex --outfile="$output"
        else
            echo "[SKIP] ${species}: existing figure $output"
        fi
        require_file "$output"
    done

    local nchr ngenes
    nchr="$(awk 'NF >= 4 {seen[$1]=1} END {for (x in seen) n++; print n+0}' "${input_dir}/${species}.bed")"
    ngenes="$(awk 'NF >= 4 {n++} END {print n+0}' "${input_dir}/${species}.bed")"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$species" "$label" "$nchr" "$ngenes" "$(anchor_count "$raw")" \
        "$(anchor_count "$screened")" "${pair_dir}/${raw}" "${pair_dir}/${screened}" "$raw_png" \
        >> "${SUMMARY_DIR}/self_synteny_summary.tsv"
}

mkdir -p "${STAGE_DIR}/canonical" "$INPUT_ROOT" "$PAIR_ROOT" "$FIG_DIR" "$SUMMARY_DIR"
require_file "$CONDA_SH"; require_file "$PREPARE_SCRIPT"

prepare_species pzi \
    "${PZI_DIR}/gff/pzi_hapA1.gff3" "${PZI_DIR}/pep/pzi_hapA1.pep.fa" \
    "${PZI_DIR}/gff/pzi_hapB1.gff3" "${PZI_DIR}/pep/pzi_hapB1.pep.fa"
prepare_species psp \
    "${PSP_GENES_DIR}/hapA1_helixer_edta_masked/postprocess/psp_hapA1.edta_masked.helixer.genes.gff3" \
    "${PSP_GENES_DIR}/hapA1_helixer_edta_masked/postprocess/psp_hapA1.edta_masked.helixer.protein.fa" \
    "${PSP_GENES_DIR}/hapB1_helixer_edta_masked/postprocess/psp_hapB1.edta_masked.helixer.genes.gff3" \
    "${PSP_GENES_DIR}/hapB1_helixer_edta_masked/postprocess/psp_hapB1.edta_masked.helixer.protein.fa"

source "$CONDA_SH"
conda activate "$JCVI_ENV"
export MPLBACKEND="${MPLBACKEND:-Agg}"
export MPLCONFIGDIR="${OUT_DIR}/mplconfig"
mkdir -p "$MPLCONFIGDIR"
command -v python >/dev/null 2>&1 || die "python is unavailable in ${JCVI_ENV}"
command -v diamond >/dev/null 2>&1 || die "diamond is unavailable in ${JCVI_ENV}"
python -c 'import jcvi' >/dev/null 2>&1 || die "jcvi is not importable in ${JCVI_ENV}"

echo "=== PZI and PSP full hapA1 + hapB1 self-synteny with JCVI ==="
echo "Output: ${OUT_DIR}"
echo "JCVI: ${JCVI_ENV}; DIAMOND: $(command -v diamond)"
echo "C-score=${CSCORE}; min_size=${MIN_ANCHORS}; minspan=${MINSPAN}"
printf 'species\tlabel\tchromosomes\tgenes\traw_anchor_pairs\tscreened_anchor_pairs\traw_anchors\tscreened_anchors\traw_dotplot_png\n' > "${SUMMARY_DIR}/self_synteny_summary.tsv"
run_self pzi
run_self psp

cat > "${OUT_DIR}/README.txt" <<EOF
Full copy1 self-synteny for PZI and PSP with JCVI.
Each input contains hapA1 + hapB1 and retains C01_A1--C08_A1/C01_B1--C08_B1 IDs.
JCVI parameters: cscore=${CSCORE}, min_size=${MIN_ANCHORS}, minspan=${MINSPAN}.
Raw and screened PNG/PDF dotplots are in ${FIG_DIR}.
EOF
echo "[OK] PZI and PSP full hapA1/hapB1 self-synteny completed."
echo "Summary: ${SUMMARY_DIR}/self_synteny_summary.tsv"
echo "Figures: ${FIG_DIR}"
