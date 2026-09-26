#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

# Pairwise JCVI synteny among peach (PPE), almond (PCAM) and GTR.
# The three comparisons are: ppe-pcam, ppe-gtr and pcam-gtr.

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/downstream_analysis/karyotype/ppe_pcam_gtr_pairwise}"
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

stage_species() {
    local label="$1" genome="$2" gff="$3" pep="$4" cds="$5"
    local dest="${STAGE_DIR}/canonical/${label}"
    for f in "$genome" "$gff" "$pep" "$cds"; do require_file "$f"; done
    mkdir -p "$dest"
    ln -sfn "$genome" "${dest}/genome.fa"
    ln -sfn "$gff" "${dest}/primary.gff3"
    ln -sfn "$pep" "${dest}/primary.pep.fa"
    ln -sfn "$cds" "${dest}/primary.cds.fa"
}

run_pair() {
    local query="$1" ref="$2"
    local comparison="${query}_vs_${ref}"
    local pair_dir="${PAIR_ROOT}/${comparison}"
    local qbed="${query}.bed" rbed="${ref}.bed"
    local qpep="${query}.pep" rpep="${ref}.pep"
    local raw="${query}.${ref}.anchors"
    local screened="${query}.${ref}.anchors.new"
    local simple="${query}.${ref}.anchors.simple"
    local raw_png="${FIG_DIR}/${comparison}.raw.dotplot.png"
    local raw_pdf="${FIG_DIR}/${comparison}.raw.dotplot.pdf"
    local screened_png="${FIG_DIR}/${comparison}.screened.dotplot.png"
    local screened_pdf="${FIG_DIR}/${comparison}.screened.dotplot.pdf"

    mkdir -p "$pair_dir"
    require_file "${INPUT_ROOT}/${query}/${qbed}"
    require_file "${INPUT_ROOT}/${ref}/${rbed}"
    require_file "${INPUT_ROOT}/${query}/${qpep}"
    require_file "${INPUT_ROOT}/${ref}/${rpep}"
    ln -sfn "${INPUT_ROOT}/${query}/${qbed}" "${pair_dir}/${qbed}"
    ln -sfn "${INPUT_ROOT}/${ref}/${rbed}" "${pair_dir}/${rbed}"
    ln -sfn "${INPUT_ROOT}/${query}/${qpep}" "${pair_dir}/${qpep}"
    ln -sfn "${INPUT_ROOT}/${ref}/${rpep}" "${pair_dir}/${rpep}"

    if [[ "$FORCE" == "1" ]]; then
        rm -f "${pair_dir}/${raw}" "${pair_dir}/${raw}.simple" \
            "${pair_dir}/${screened}" "${pair_dir}/${screened}.simple" \
            "${pair_dir}/${simple}" "$raw_png" "$raw_pdf" "$screened_png" "$screened_pdf"
    fi

    (
        cd "$pair_dir"
        if [[ ! -s "$raw" ]]; then
            echo "[RUN] ${comparison}: JCVI catalog ortholog"
            python -m jcvi.compara.catalog ortholog "$query" "$ref" \
                --dbtype=prot --align_soft=diamond_blastp \
                --cscore="$CSCORE" --min_size="$MIN_ANCHORS" \
                --cpus="$THREADS" --no_strip_names --no_dotplot
        else
            echo "[SKIP] ${comparison}: existing raw anchors"
        fi
        require_file "$raw"

        if [[ ! -s "$screened" ]]; then
            python -m jcvi.compara.synteny screen \
                --qbed="$qbed" --sbed="$rbed" --minspan="$MINSPAN" \
                --simple "$raw" "$screened"
        else
            echo "[SKIP] ${comparison}: existing screened anchors"
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
                    --qbed="$qbed" --sbed="$rbed" --figsize=12x12 --dpi=300 \
                    --format="$format" --nochpf --notex --outfile="$output"
            else
                echo "[SKIP] ${comparison}: existing figure $output"
            fi
            require_file "$output"
        done
    )

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$comparison" "$query" "$ref" "$(anchor_count "${pair_dir}/${raw}")" \
        "$(anchor_count "${pair_dir}/${screened}")" "${pair_dir}/${raw}" \
        "${pair_dir}/${screened}" "$raw_png" "$raw_pdf" "$screened_png" \
        >> "${SUMMARY_DIR}/pairwise_summary.tsv"
}

mkdir -p "${STAGE_DIR}/canonical" "$INPUT_ROOT" "$PAIR_ROOT" "$FIG_DIR" "$SUMMARY_DIR"
require_file "$CONDA_SH"
require_file "$PREPARE_SCRIPT"

stage_species pcam \
    "${ROOT_DIR}/data/pcam/GCA_050084535.1_PCAv2.0_genomic.fna.gz" \
    "${ROOT_DIR}/data/pcam/GCA_050084535.1_PCAv2.0_genomic.gff.gz" \
    "${ROOT_DIR}/data/pcam/GCA_050084535.1_PCAv2.0_protein.faa.gz" \
    "${ROOT_DIR}/data/pcam/GCA_050084535.1_PCAv2.0_cds_from_genomic.fna.gz"
stage_species ppe \
    "${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fna.gz" \
    "${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_genomic.gff.gz" \
    "${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_protein.faa.gz" \
    "${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_cds_from_genomic.fna.gz"
stage_species gtr \
    "${ROOT_DIR}/data/gtr/GCA_018257905.1_V1_genomic.fna" \
    "${ROOT_DIR}/data/gtr/Gtr.NCBI.gene.gff3" \
    "${ROOT_DIR}/data/gtr/Gtr.NCBI.prot.faa" \
    "${ROOT_DIR}/data/gtr/Gtr.NCBI.gene.modId.fna"

source "$CONDA_SH"
conda activate "$JCVI_ENV"
export MPLBACKEND="${MPLBACKEND:-Agg}"
export MPLCONFIGDIR="${OUT_DIR}/mplconfig"
mkdir -p "$MPLCONFIGDIR"
command -v python >/dev/null 2>&1 || die "python is unavailable in ${JCVI_ENV}"
command -v diamond >/dev/null 2>&1 || die "diamond is unavailable in ${JCVI_ENV}"
python -c 'import jcvi' >/dev/null 2>&1 || die "jcvi is not importable in ${JCVI_ENV}"

python "${PREPARE_SCRIPT}" \
    --analysis-dir "${STAGE_DIR}" \
    --outdir "${INPUT_ROOT}" \
    --summary "${SUMMARY_DIR}/input_summary.tsv" \
    --labels pcam ppe gtr

echo "=== JCVI pairwise synteny: PPE, PCAM and GTR ==="
echo "Output: ${OUT_DIR}"
echo "Parameters: cscore=${CSCORE}, min_size=${MIN_ANCHORS}, minspan=${MINSPAN}"
printf 'comparison\tquery\tref\traw_anchor_pairs\tscreened_anchor_pairs\traw_anchors\tscreened_anchors\traw_png\traw_pdf\tscreened_png\n' \
    > "${SUMMARY_DIR}/pairwise_summary.tsv"

run_pair ppe pcam
run_pair ppe gtr
run_pair pcam gtr

cat > "${OUT_DIR}/README.txt" <<EOF
Pairwise JCVI synteny for ppe, pcam and gtr.
Comparisons: ppe-pcam, ppe-gtr and pcam-gtr.
Parameters: cscore=${CSCORE}, min_size=${MIN_ANCHORS}, minspan=${MINSPAN}.
EOF
echo "[OK] Three pairwise JCVI dotplots completed."
echo "Summary: ${SUMMARY_DIR}/pairwise_summary.tsv"
