#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

# Compare only C02_A and C06_A within each species:
#   PZI: C02_A1 vs C06_A1
#   PSP: C02_A1 vs C06_A1
#   PPA: C02_A  vs C06_A
# Each comparison is independent and produces raw and screened JCVI dotplots.

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/haploid_rename}"
PPA_DIR="${PPA_DIR:-${ROOT_DIR}/data/ppa/haploid_rename}"
PSP_GENES_DIR="${PSP_GENES_DIR:-${ROOT_DIR}/data/psp/haploid/annotation/genes}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/downstream_analysis/karyotype/C02A_C06A_self_jcvi}"
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

extract_fasta_records() {
    local ids="$1"
    local fasta="$2"
    local output="$3"
    awk -v ids="$ids" '
        BEGIN {
            while ((getline line < ids) > 0) wanted[line] = 1
            close(ids)
        }
        /^>/ {
            if (keep && seq != "") print seq
            header = substr($0, 2)
            split(header, a, /[[:space:]]+/)
            keep = (a[1] in wanted)
            seq = keep ? $0 "\n" : ""
            next
        }
        {if (keep) seq = seq $0 "\n"}
        END {if (keep && seq != "") print seq}
    ' "$fasta" > "$output"
}

prepare_species() {
    local species="$1"
    local gff="$2"
    local pep="$3"
    local chr_a="$4"
    local chr_b="$5"
    local canonical="${STAGE_DIR}/canonical/${species}"
    local ids="${canonical}/protein_ids.txt"

    require_file "$gff"
    require_file "$pep"
    mkdir -p "$canonical"

    if [[ "$FORCE" == "1" || ! -s "${canonical}/primary.gff3" ]]; then
        {
            printf '##gff-version 3\n'
            awk -F '\t' -v a="$chr_a" -v b="$chr_b" \
                '$0 ~ /^#/ {next} NF >= 9 && ($1 == a || $1 == b) {print}' "$gff"
        } > "${canonical}/primary.gff3.tmp"
        mv "${canonical}/primary.gff3.tmp" "${canonical}/primary.gff3"
    fi

    # Keep only mRNA/protein IDs represented by the two requested chromosomes.
    awk -F '\t' -v a="$chr_a" -v b="$chr_b" '
        $0 !~ /^#/ && NF >= 9 && ($1 == a || $1 == b) && $3 == "mRNA" {
            if (match($9, /(^|;)ID=[^;]+/)) {
                x = substr($9, RSTART, RLENGTH)
                sub(/^.*ID=/, "", x)
                print x
            }
        }
    ' "$gff" | sort -u > "${ids}.tmp"
    mv "${ids}.tmp" "$ids"
    [[ -s "$ids" ]] || die "No mRNA IDs found for ${species} ${chr_a}/${chr_b}"

    if [[ "$FORCE" == "1" || ! -s "${canonical}/primary.pep.fa" ]]; then
        extract_fasta_records "$ids" "$pep" "${canonical}/primary.pep.fa.tmp"
        mv "${canonical}/primary.pep.fa.tmp" "${canonical}/primary.pep.fa"
    fi
    require_file "${canonical}/primary.gff3"
    require_file "${canonical}/primary.pep.fa"
}

plot_dotplot() {
    local pair_dir="$1"
    local anchors="$2"
    local bed="$3"
    local output="$4"
    local format="$5"
    (
        cd "$pair_dir"
        python -m jcvi.graphics.dotplot "$anchors" \
            --qbed="$bed" \
            --sbed="$bed" \
            --figsize=12x12 \
            --dpi=300 \
            --format="$format" \
            --nochpf \
            --notex \
            --outfile="$output"
    )
}

run_self() {
    local species="$1"
    local canonical="${STAGE_DIR}/canonical/${species}"
    local input_dir="${INPUT_ROOT}/${species}"
    local pair_dir="${PAIR_ROOT}/${species}_C02A_C06A"
    local label="${species}_C02A_C06A"
    local raw="${label}.${label}.anchors"
    local screened="${label}.${label}.anchors.new"
    local simple="${label}.${label}.anchors.simple"
    local bed="${label}.bed"
    local pep="${label}.pep"
    local raw_png="${FIG_DIR}/${species}_C02A_C06A.raw.dotplot.png"
    local raw_pdf="${FIG_DIR}/${species}_C02A_C06A.raw.dotplot.pdf"
    local screened_png="${FIG_DIR}/${species}_C02A_C06A.screened.dotplot.png"
    local screened_pdf="${FIG_DIR}/${species}_C02A_C06A.screened.dotplot.pdf"

    mkdir -p "$input_dir" "$pair_dir"
    python "$PREPARE_SCRIPT" \
        --analysis-dir "$STAGE_DIR" \
        --outdir "$INPUT_ROOT" \
        --summary "${SUMMARY_DIR}/${species}_input_summary.tsv" \
        --labels "$species"

    require_file "${input_dir}/${species}.bed"
    require_file "${input_dir}/${species}.pep"
    ln -sfn "${input_dir}/${species}.bed" "${pair_dir}/${bed}"
    ln -sfn "${input_dir}/${species}.pep" "${pair_dir}/${pep}"

    if [[ "$FORCE" == "1" ]]; then
        rm -f "${pair_dir}/${raw}" "${pair_dir}/${raw}.simple" \
            "${pair_dir}/${screened}" "${pair_dir}/${screened}.simple" \
            "${pair_dir}/${simple}" "$raw_png" "$raw_pdf" "$screened_png" "$screened_pdf"
    fi

    (
        cd "$pair_dir"
        if [[ ! -s "$raw" ]]; then
            echo "[RUN] ${species}: JCVI catalog ortholog ${label} ${label}"
            python -m jcvi.compara.catalog ortholog "$label" "$label" \
                --dbtype=prot \
                --align_soft=diamond_blastp \
                --cscore="$CSCORE" \
                --min_size="$MIN_ANCHORS" \
                --cpus="$THREADS" \
                --no_strip_names \
                --no_dotplot
        else
            echo "[SKIP] ${species}: existing raw anchors"
        fi
        require_file "$raw"
        if [[ ! -s "$screened" ]]; then
            python -m jcvi.compara.synteny screen \
                --qbed="$bed" --sbed="$bed" \
                --minspan="$MINSPAN" --simple "$raw" "$screened"
        else
            echo "[SKIP] ${species}: existing screened anchors"
        fi
        require_file "$screened"
        if [[ ! -s "$simple" ]]; then
            for candidate in "${screened}.simple" "${raw}.simple"; do
                if [[ -s "$candidate" ]]; then ln -sfn "$candidate" "$simple"; break; fi
            done
        fi
    )

    [[ -s "$raw_png" ]] || plot_dotplot "$pair_dir" "$raw" "$bed" "$raw_png" png
    [[ -s "$raw_pdf" ]] || plot_dotplot "$pair_dir" "$raw" "$bed" "$raw_pdf" pdf
    [[ -s "$screened_png" ]] || plot_dotplot "$pair_dir" "$screened" "$bed" "$screened_png" png
    [[ -s "$screened_pdf" ]] || plot_dotplot "$pair_dir" "$screened" "$bed" "$screened_pdf" pdf

    for f in "$raw_png" "$raw_pdf" "$screened_png" "$screened_pdf"; do require_file "$f"; done
    local nchr nbed
    nchr="$(awk 'NF >= 4 {seen[$1]=1} END {for (x in seen) n++; print n+0}' "${input_dir}/${species}.bed")"
    nbed="$(awk 'NF >= 4 {n++} END {print n+0}' "${input_dir}/${species}.bed")"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$species" "$nchr" "$nbed" "$(anchor_count "${pair_dir}/${raw}")" \
        "$(anchor_count "${pair_dir}/${screened}")" \
        "${pair_dir}/${raw}" "${pair_dir}/${screened}" "${raw_png}" \
        >> "${SUMMARY_DIR}/self_synteny_summary.tsv"
}

mkdir -p "$STAGE_DIR/canonical" "$INPUT_ROOT" "$PAIR_ROOT" "$FIG_DIR" "$SUMMARY_DIR"
require_file "$CONDA_SH"
require_file "$PREPARE_SCRIPT"

# PZI and PSP use A1 chromosome IDs; PPA's renamed A haplotype uses C##_A.
prepare_species pzi \
    "${PZI_DIR}/gff/pzi_hapA1.gff3" \
    "${PZI_DIR}/pep/pzi_hapA1.pep.fa" C02_A1 C06_A1
prepare_species psp \
    "${PSP_GENES_DIR}/hapA1_helixer_edta_masked/postprocess/psp_hapA1.edta_masked.helixer.genes.gff3" \
    "${PSP_GENES_DIR}/hapA1_helixer_edta_masked/postprocess/psp_hapA1.edta_masked.helixer.protein.fa" C02_A1 C06_A1
prepare_species ppa \
    "${PPA_DIR}/gff/ppa_hapA.gff3" \
    "${PPA_DIR}/pep/ppa_hapA.pep.fa" C02_A C06_A

source "$CONDA_SH"
conda activate "$JCVI_ENV"
export MPLBACKEND="${MPLBACKEND:-Agg}"
export MPLCONFIGDIR="${OUT_DIR}/mplconfig"
mkdir -p "$MPLCONFIGDIR"
command -v python >/dev/null 2>&1 || die "python is unavailable in ${JCVI_ENV}"
command -v diamond >/dev/null 2>&1 || die "diamond is unavailable in ${JCVI_ENV}"
python -c 'import jcvi' >/dev/null 2>&1 || die "jcvi is not importable in ${JCVI_ENV}"

echo "=== C02_A vs C06_A self-synteny: PZI, PSP, PPA ==="
echo "Output: $OUT_DIR"
echo "JCVI: ${JCVI_ENV}; DIAMOND: $(command -v diamond)"
echo "C-score=${CSCORE}; min_size=${MIN_ANCHORS}; minspan=${MINSPAN}"
printf 'species\tchromosomes\tgenes\traw_anchor_pairs\tscreened_anchor_pairs\traw_anchors\tscreened_anchors\traw_dotplot_png\n' > "${SUMMARY_DIR}/self_synteny_summary.tsv"

run_self pzi
run_self psp
run_self ppa

cat > "${OUT_DIR}/README.txt" <<EOF
C02_A versus C06_A JCVI self-synteny for PZI, PSP and PPA.
PZI/PSP input chromosomes: C02_A1 and C06_A1.
PPA input chromosomes: C02_A and C06_A.
JCVI parameters: cscore=${CSCORE}, min_size=${MIN_ANCHORS}, minspan=${MINSPAN}.
Raw and screened PNG/PDF dotplots are in ${FIG_DIR}.
EOF

echo "[OK] Three C02_A/C06_A self-comparisons completed."
echo "Summary: ${SUMMARY_DIR}/self_synteny_summary.tsv"
echo "Figures: ${FIG_DIR}"
