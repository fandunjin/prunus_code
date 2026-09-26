#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

# JCVI self-synteny and dotplot for the two PPA haplotypes (ppa_hapA +
# ppa_hapB).  The two renamed GFF3 files are combined into one annotation;
# JCVI then searches for duplicated/homeologous blocks within PPA.

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/data/ppa/haploid_rename}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/downstream_analysis/karyotype/ppa_self_jcvi}"
STAGE_DIR="${OUT_DIR}/stage"
CANONICAL_DIR="${STAGE_DIR}/canonical/ppa"
INPUT_ROOT="${OUT_DIR}/inputs"
PAIR_DIR="${OUT_DIR}/pairs/ppa_self"
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
LABEL="ppa"
RAW="${LABEL}.${LABEL}.anchors"
SCREENED="${LABEL}.${LABEL}.anchors.new"
SIMPLE="${LABEL}.${LABEL}.anchors.simple"

RAW_PNG="${FIG_DIR}/ppa_self.raw.dotplot.png"
RAW_PDF="${FIG_DIR}/ppa_self.raw.dotplot.pdf"
SCREENED_PNG="${FIG_DIR}/ppa_self.screened.dotplot.png"
SCREENED_PDF="${FIG_DIR}/ppa_self.screened.dotplot.pdf"

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
    local output="$2"
    local format="$3"
    python -m jcvi.graphics.dotplot "$anchors" \
        --qbed="ppa.bed" \
        --sbed="ppa.bed" \
        --figsize=12x12 \
        --dpi=300 \
        --format="${format}" \
        --nochpf \
        --notex \
        --outfile="${output}"
}

mkdir -p "${CANONICAL_DIR}" "${INPUT_ROOT}" "${PAIR_DIR}" "${FIG_DIR}" "${SUMMARY_DIR}"
require_file "${CONDA_SH}"
require_file "${PREPARE_SCRIPT}"
for hap in hapA hapB; do
    require_file "${PPA_DIR}/gff/ppa_${hap}.gff3"
    require_file "${PPA_DIR}/pep/ppa_${hap}.pep.fa"
done

# Build canonical files only when missing (set FORCE=1 to rebuild them).
if [[ "${FORCE}" == "1" || ! -s "${CANONICAL_DIR}/primary.gff3" ]]; then
    {
        printf '##gff-version 3\n'
        for hap in hapA hapB; do
            awk '!/^#/ && NF >= 9 {print}' "${PPA_DIR}/gff/ppa_${hap}.gff3"
        done
    } > "${CANONICAL_DIR}/primary.gff3.tmp"
    mv "${CANONICAL_DIR}/primary.gff3.tmp" "${CANONICAL_DIR}/primary.gff3"
fi
if [[ "${FORCE}" == "1" || ! -s "${CANONICAL_DIR}/primary.pep.fa" ]]; then
    cat "${PPA_DIR}/pep/ppa_hapA.pep.fa" "${PPA_DIR}/pep/ppa_hapB.pep.fa" \
        > "${CANONICAL_DIR}/primary.pep.fa.tmp"
    mv "${CANONICAL_DIR}/primary.pep.fa.tmp" "${CANONICAL_DIR}/primary.pep.fa"
fi

# shellcheck disable=SC1090
source "${CONDA_SH}"
conda activate "${JCVI_ENV}"
export MPLBACKEND="${MPLBACKEND:-Agg}"
export MPLCONFIGDIR="${OUT_DIR}/mplconfig"
mkdir -p "${MPLCONFIGDIR}"

command -v python >/dev/null 2>&1 || die "python is unavailable in ${JCVI_ENV}"
command -v diamond >/dev/null 2>&1 || die "diamond is unavailable in ${JCVI_ENV}"
python -c 'import jcvi' >/dev/null 2>&1 || die "jcvi is not importable in ${JCVI_ENV}"

echo "=== PPA self-synteny with JCVI ==="
echo "PPA input   : ${PPA_DIR}"
echo "Output      : ${OUT_DIR}"
echo "JCVI env    : ${JCVI_ENV}"
echo "DIAMOND     : $(command -v diamond)"
echo "C-score     : ${CSCORE}"
echo "Min anchors : ${MIN_ANCHORS}"
echo "Min span    : ${MINSPAN}"

python "${PREPARE_SCRIPT}" \
    --analysis-dir "${STAGE_DIR}" \
    --outdir "${INPUT_ROOT}" \
    --summary "${SUMMARY_DIR}/input_summary.tsv" \
    --labels "${LABEL}"

BED="${INPUT_ROOT}/${LABEL}/${LABEL}.bed"
PEP="${INPUT_ROOT}/${LABEL}/${LABEL}.pep"
require_file "${BED}"
require_file "${PEP}"
ln -sfn "${BED}" "${PAIR_DIR}/ppa.bed"
ln -sfn "${PEP}" "${PAIR_DIR}/ppa.pep"

if [[ "${FORCE}" == "1" ]]; then
    rm -f "${PAIR_DIR}/${RAW}" "${PAIR_DIR}/${RAW}.simple" \
        "${PAIR_DIR}/${SCREENED}" "${PAIR_DIR}/${SCREENED}.simple" \
        "${PAIR_DIR}/${SIMPLE}" "${RAW_PNG}" "${RAW_PDF}" \
        "${SCREENED_PNG}" "${SCREENED_PDF}"
fi

cd "${PAIR_DIR}"
if [[ ! -s "${RAW}" ]]; then
    echo "[RUN] JCVI catalog ortholog ppa ppa"
    python -m jcvi.compara.catalog ortholog ppa ppa \
        --dbtype=prot \
        --align_soft=diamond_blastp \
        --cscore="${CSCORE}" \
        --min_size="${MIN_ANCHORS}" \
        --cpus="${THREADS}" \
        --no_strip_names \
        --no_dotplot
else
    echo "[SKIP] Existing raw anchors: ${PAIR_DIR}/${RAW}"
fi
require_file "${RAW}"

if [[ ! -s "${SCREENED}" ]]; then
    echo "[RUN] Screen PPA self-synteny anchors"
    python -m jcvi.compara.synteny screen \
        --qbed=ppa.bed \
        --sbed=ppa.bed \
        --minspan="${MINSPAN}" \
        --simple \
        "${RAW}" "${SCREENED}"
else
    echo "[SKIP] Existing screened anchors: ${PAIR_DIR}/${SCREENED}"
fi
require_file "${SCREENED}"

if [[ ! -s "${SIMPLE}" ]]; then
    for candidate in "${SCREENED}.simple" "${RAW}.simple"; do
        if [[ -s "${candidate}" ]]; then
            ln -sfn "${candidate}" "${SIMPLE}"
            break
        fi
    done
fi

[[ -s "${RAW_PNG}" ]] || plot_dotplot "${RAW}" "${RAW_PNG}" png
[[ -s "${RAW_PDF}" ]] || plot_dotplot "${RAW}" "${RAW_PDF}" pdf
[[ -s "${SCREENED_PNG}" ]] || plot_dotplot "${SCREENED}" "${SCREENED_PNG}" png
[[ -s "${SCREENED_PDF}" ]] || plot_dotplot "${SCREENED}" "${SCREENED_PDF}" pdf

for f in "${RAW_PNG}" "${RAW_PDF}" "${SCREENED_PNG}" "${SCREENED_PDF}"; do
    require_file "${f}"
done

N_CHROMOSOMES="$(awk 'NF >= 4 {seen[$1]=1} END {for (x in seen) n++; print n+0}' "${BED}")"
N_GENES="$(awk 'NF >= 4 {n++} END {print n+0}' "${BED}")"
printf 'label\tcomponents\tchromosomes\tgenes\traw_anchor_pairs\tscreened_anchor_pairs\traw_anchors\tscreened_anchors\traw_dotplot\tscreened_dotplot\n' > "${SUMMARY_DIR}/self_synteny_summary.tsv"
printf 'ppa\thapA+hapB\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${N_CHROMOSOMES}" "${N_GENES}" "$(anchor_count "${RAW}")" \
    "$(anchor_count "${SCREENED}")" "${PAIR_DIR}/${RAW}" "${PAIR_DIR}/${SCREENED}" \
    "${RAW_PNG}" "${SCREENED_PNG}" >> "${SUMMARY_DIR}/self_synteny_summary.tsv"

cat > "${OUT_DIR}/README.txt" <<EOF
PPA self-synteny generated with JCVI.
Input: ${PPA_DIR}/gff/ppa_hapA.gff3 + ppa_hapB.gff3 and matching protein FASTA files.
JCVI parameters: cscore=${CSCORE}, min_size=${MIN_ANCHORS}, minspan=${MINSPAN}.
Raw and screened dotplots are in ${FIG_DIR}.
EOF

echo "[OK] PPA self-synteny completed."
echo "Summary: ${SUMMARY_DIR}/self_synteny_summary.tsv"
echo "Figures: ${FIG_DIR}"
