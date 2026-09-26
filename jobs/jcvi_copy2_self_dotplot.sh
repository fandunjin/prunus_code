#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

# Gene-level self-synteny of PZI copy2. The copy2 annotation is assembled from
# the chromosome-renamed primary hapA2 and hapB2 GFF/CDS/protein sets.

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/haploid_rename}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/downstream_analysis/karyotype/copy2_self_jcvi}"
STAGE_DIR="${OUT_DIR}/stage"
CANONICAL_DIR="${STAGE_DIR}/canonical/pzi_copy2"
INPUT_ROOT="${OUT_DIR}/inputs"
INPUT_DIR="${INPUT_ROOT}/pzi_copy2"
PAIR_DIR="${OUT_DIR}/pairs/pzi_copy2_self"
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
PACKAGE_CHECKSUMS="${HAP_RENAME_DIR}/checksums.md5"
CHROMOSOME_MAP="${HAP_RENAME_DIR}/chromosome_rename_map.tsv"
COPY2_SEQID_MAP="${SUMMARY_DIR}/copy2_seqid_map.tsv"
SIGNATURE_FILE="${SUMMARY_DIR}/copy2_self_jcvi.signature"
PLOT_SIGNATURE_FILE="${SUMMARY_DIR}/copy2_dotplot.signature"
INPUT_SUMMARY="${SUMMARY_DIR}/input_summary.tsv"
RUN_SUMMARY="${SUMMARY_DIR}/self_synteny_summary.tsv"

LABEL="pzi_copy2"
RAW="${LABEL}.${LABEL}.anchors"
SCREENED="${LABEL}.${LABEL}.anchors.new"
SIMPLE="${LABEL}.${LABEL}.anchors.simple"

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

    python -m jcvi.graphics.dotplot "${anchors}" \
        --qbed="${LABEL}.bed" \
        --sbed="${LABEL}.bed" \
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
require_file "${PACKAGE_CHECKSUMS}"
require_file "${CHROMOSOME_MAP}"

for hap in hapA2 hapB2; do
    require_file "${HAP_RENAME_DIR}/gff/pzi_${hap}.gff3"
    require_file "${HAP_RENAME_DIR}/cds/pzi_${hap}.cds.fa"
    require_file "${HAP_RENAME_DIR}/pep/pzi_${hap}.pep.fa"
done

awk -F '\t' '
    BEGIN {OFS = "\t"}
    NR == 1 {next}
    $6 != "" {
        if ($3 in seen_source) {
            print "Duplicate copy2 source chromosome: " $3 > "/dev/stderr"
            exit 1
        }
        if ($6 in seen_target) {
            print "Duplicate copy2 target chromosome: " $6 > "/dev/stderr"
            exit 1
        }
        seen_source[$3] = 1
        seen_target[$6] = 1
        print $3, $6
        n++
    }
    END {
        if (n != 16) {
            print "Expected 16 copy2 chromosome mappings, observed " n > "/dev/stderr"
            exit 1
        }
    }
' "${CHROMOSOME_MAP}" > "${COPY2_SEQID_MAP}.tmp"
mv "${COPY2_SEQID_MAP}.tmp" "${COPY2_SEQID_MAP}"

PACKAGE_ID="$(cksum "${PACKAGE_CHECKSUMS}" | awk '{print $1 ":" $2}')"
MAP_ID="$(cksum "${COPY2_SEQID_MAP}" | awk '{print $1 ":" $2}')"
SIGNATURE="package=${PACKAGE_ID};cscore=${CSCORE};min_anchors=${MIN_ANCHORS};minspan=${MINSPAN}"
PLOT_SIGNATURE="copy2_chr01_chr16_labels=1;map=${MAP_ID};figsize=12x12;dpi=300"
REBUILD=0
REPLOT=0
if [[ "${FORCE}" == "1" ]]; then
    REBUILD=1
elif [[ -s "${SIGNATURE_FILE}" ]] && ! grep -Fqx -- "${SIGNATURE}" "${SIGNATURE_FILE}"; then
    REBUILD=1
fi
if [[ "${FORCE}" == "1" || ! -s "${PLOT_SIGNATURE_FILE}" ]] || \
   ! grep -Fqx -- "${PLOT_SIGNATURE}" "${PLOT_SIGNATURE_FILE}"; then
    REPLOT=1
fi

echo "=== PZI copy2 self-synteny with JCVI ==="
echo "Components : hapA2 + hapB2"
echo "Hap package: ${HAP_RENAME_DIR}"
echo "Output     : ${OUT_DIR}"
echo "JCVI env   : ${JCVI_ENV}"
echo "Threads    : ${THREADS}"
echo "C-score    : ${CSCORE}"
echo "Min anchors: ${MIN_ANCHORS}"
echo "Min span   : ${MINSPAN}"
echo "Package ID : ${PACKAGE_ID}"
echo "Display IDs: chr01--chr16 (copy2 convention)"
echo "Rebuild    : ${REBUILD}"
echo "Redraw     : ${REPLOT}"

# Build one canonical copy2 annotation. Only GFF column 1 is renamed to the
# copy2 chr01--chr16 convention; identifiers, coordinates and strands remain
# unchanged.
{
    printf '##gff-version 3\n'
    for hap in hapA2 hapB2; do
        awk -F '\t' '
            BEGIN {OFS = "\t"}
            NR == FNR {map[$1] = $2; next}
            $0 ~ /^#/ || NF == 0 {next}
            !($1 in map) {
                print "Unmapped copy2 GFF sequence ID: " $1 > "/dev/stderr"
                exit 1
            }
            {$1 = map[$1]; print}
        ' "${COPY2_SEQID_MAP}" "${HAP_RENAME_DIR}/gff/pzi_${hap}.gff3"
    done
} > "${CANONICAL_DIR}/primary.gff3.tmp"
mv "${CANONICAL_DIR}/primary.gff3.tmp" "${CANONICAL_DIR}/primary.gff3"

cat \
    "${HAP_RENAME_DIR}/pep/pzi_hapA2.pep.fa" \
    "${HAP_RENAME_DIR}/pep/pzi_hapB2.pep.fa" \
    > "${CANONICAL_DIR}/primary.pep.fa.tmp"
mv "${CANONICAL_DIR}/primary.pep.fa.tmp" "${CANONICAL_DIR}/primary.pep.fa"

cat \
    "${HAP_RENAME_DIR}/cds/pzi_hapA2.cds.fa" \
    "${HAP_RENAME_DIR}/cds/pzi_hapB2.cds.fa" \
    > "${CANONICAL_DIR}/primary.cds.fa.tmp"
mv "${CANONICAL_DIR}/primary.cds.fa.tmp" "${CANONICAL_DIR}/primary.cds.fa"

# shellcheck disable=SC1090
source "${CONDA_SH}"
conda activate "${JCVI_ENV}"
export MPLBACKEND="${MPLBACKEND:-Agg}"
export MPLCONFIGDIR="${OUT_DIR}/mplconfig"
mkdir -p "${MPLCONFIGDIR}"

command -v python >/dev/null 2>&1 || die "python is unavailable in ${JCVI_ENV}"
command -v diamond >/dev/null 2>&1 || die "diamond is unavailable in ${JCVI_ENV}"
python -c 'import jcvi' >/dev/null 2>&1 || die "jcvi is not importable in ${JCVI_ENV}"

python "${PREPARE_SCRIPT}" \
    --analysis-dir "${STAGE_DIR}" \
    --outdir "${INPUT_ROOT}" \
    --summary "${INPUT_SUMMARY}" \
    --labels "${LABEL}"

BED="${INPUT_DIR}/${LABEL}.bed"
PEP="${INPUT_DIR}/${LABEL}.pep"
require_file "${BED}"
require_file "${PEP}"

N_CHROMOSOMES="$(awk 'NF >= 4 {seen[$1]=1} END {print length(seen)}' "${BED}")"
N_GENES="$(awk 'NF >= 4 {n++} END {print n+0}' "${BED}")"
DUPLICATE_IDS="$(awk 'NF >= 4 {if (seen[$4]++) n++} END {print n+0}' "${BED}")"
EXPECTED_SEQIDS="chr01,chr02,chr03,chr04,chr05,chr06,chr07,chr08,chr09,chr10,chr11,chr12,chr13,chr14,chr15,chr16"
OBSERVED_SEQIDS="$(awk 'NF >= 4 && !seen[$1]++ {print $1}' "${BED}" | paste -sd, -)"
[[ "${N_CHROMOSOMES}" == "16" ]] || die "Expected 16 copy2 chromosomes, observed ${N_CHROMOSOMES}"
[[ "${DUPLICATE_IDS}" == "0" ]] || die "Observed ${DUPLICATE_IDS} duplicate JCVI gene IDs"
[[ "${OBSERVED_SEQIDS}" == "${EXPECTED_SEQIDS}" ]] || \
    die "Unexpected copy2 chromosome IDs/order: ${OBSERVED_SEQIDS}"

ln -sfn "${BED}" "${PAIR_DIR}/${LABEL}.bed"
ln -sfn "${PEP}" "${PAIR_DIR}/${LABEL}.pep"

RAW_PNG="${FIG_DIR}/${LABEL}_self.raw.dotplot.png"
RAW_PDF="${FIG_DIR}/${LABEL}_self.raw.dotplot.pdf"
SCREENED_PNG="${FIG_DIR}/${LABEL}_self.screened.dotplot.png"
SCREENED_PDF="${FIG_DIR}/${LABEL}_self.screened.dotplot.pdf"

if [[ "${REBUILD}" == "1" ]]; then
    rm -f \
        "${PAIR_DIR}/${LABEL}.${LABEL}.last" \
        "${PAIR_DIR}/${LABEL}.${LABEL}.last.filtered" \
        "${PAIR_DIR}/${LABEL}.dmnd" \
        "${PAIR_DIR}/${LABEL}.${LABEL}.lifted.anchors" \
        "${PAIR_DIR}/${LABEL}.${LABEL}.lifted.anchors.simple" \
        "${PAIR_DIR}/${RAW}" \
        "${PAIR_DIR}/${RAW}.simple" \
        "${PAIR_DIR}/${SCREENED}" \
        "${PAIR_DIR}/${SCREENED}.simple" \
        "${PAIR_DIR}/${SIMPLE}" \
        "${RAW_PNG}" "${RAW_PDF}" "${SCREENED_PNG}" "${SCREENED_PDF}" \
        "${SIGNATURE_FILE}"
fi
if [[ "${REPLOT}" == "1" ]]; then
    rm -f "${RAW_PNG}" "${RAW_PDF}" "${SCREENED_PNG}" "${SCREENED_PDF}" \
        "${PLOT_SIGNATURE_FILE}"
fi

# Record the input/parameter version before the expensive search. If the job is
# interrupted later, rerunning with the same signature resumes existing files.
printf '%s\n' "${SIGNATURE}" > "${SIGNATURE_FILE}"

cd "${PAIR_DIR}"
if [[ ! -s "${RAW}" ]]; then
    echo "[RUN] JCVI catalog ortholog ${LABEL} ${LABEL}"
    python -m jcvi.compara.catalog ortholog \
        "${LABEL}" "${LABEL}" \
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
    echo "[RUN] Screen copy2 self-synteny anchors"
    python -m jcvi.compara.synteny screen \
        --qbed="${LABEL}.bed" \
        --sbed="${LABEL}.bed" \
        --minspan="${MINSPAN}" \
        --simple \
        "${RAW}" "${SCREENED}"
else
    echo "[SKIP] Existing screened anchors: ${PAIR_DIR}/${SCREENED}"
fi
require_file "${SCREENED}"

if [[ ! -s "${SIMPLE}" ]]; then
    for candidate in "${RAW}.simple" "${SCREENED}.simple" "${LABEL}.${LABEL}.simple"; do
        if [[ -s "${candidate}" ]]; then
            ln -sfn "${candidate}" "${SIMPLE}"
            break
        fi
    done
fi

if [[ ! -s "${RAW_PNG}" ]]; then plot_dotplot "${RAW}" "${RAW_PNG}" png; fi
if [[ ! -s "${RAW_PDF}" ]]; then plot_dotplot "${RAW}" "${RAW_PDF}" pdf; fi
if [[ ! -s "${SCREENED_PNG}" ]]; then plot_dotplot "${SCREENED}" "${SCREENED_PNG}" png; fi
if [[ ! -s "${SCREENED_PDF}" ]]; then plot_dotplot "${SCREENED}" "${SCREENED_PDF}" pdf; fi

require_file "${RAW_PNG}"
require_file "${RAW_PDF}"
require_file "${SCREENED_PNG}"
require_file "${SCREENED_PDF}"
printf '%s\n' "${PLOT_SIGNATURE}" > "${PLOT_SIGNATURE_FILE}"

printf 'label\tcomponents\tchromosomes\tgenes\traw_anchor_pairs\tscreened_anchor_pairs\traw_anchors\tscreened_anchors\traw_dotplot\tscreened_dotplot\n' > "${RUN_SUMMARY}"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${LABEL}" "hapA2+hapB2" "${N_CHROMOSOMES}" "${N_GENES}" \
    "$(anchor_count "${RAW}")" "$(anchor_count "${SCREENED}")" \
    "${PAIR_DIR}/${RAW}" "${PAIR_DIR}/${SCREENED}" \
    "${RAW_PNG}" "${SCREENED_PNG}" >> "${RUN_SUMMARY}"

printf '%s\n' "${SIGNATURE}" > "${SIGNATURE_FILE}"

cat > "${OUT_DIR}/README.txt" <<EOF
PZI copy2 self-synteny generated with JCVI.

copy2 composition: chromosome-renamed hapA2 + hapB2 primary annotations
Display chromosome IDs: chr01--chr16 (copy2 convention)
Chromosome map: ${COPY2_SEQID_MAP}
JCVI parameters: cscore=${CSCORE}, min_size=${MIN_ANCHORS}, minspan=${MINSPAN}
Chromosomes: ${N_CHROMOSOMES}
Genes: ${N_GENES}

Raw and screened dotplots:
  ${RAW_PNG}
  ${RAW_PDF}
  ${SCREENED_PNG}
  ${SCREENED_PDF}
EOF

echo
echo "[OK] PZI copy2 self-synteny completed."
echo "Summary : ${RUN_SUMMARY}"
echo "Figures : ${FIG_DIR}"
