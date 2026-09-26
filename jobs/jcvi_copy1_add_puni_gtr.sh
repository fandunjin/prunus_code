#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/downstream_analysis/karyotype/run_20260803_205037}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
JCVI_ENV="${JCVI_ENV:-jcvi}"
THREADS="${THREADS:-32}"
CSCORE="${CSCORE:-0.70}"
MIN_ANCHORS="${MIN_ANCHORS:-4}"
MINSPAN="${MINSPAN:-30}"

SCRIPTS_DIR="${ROOT_DIR}/scripts"
PREPARE_SCRIPT="${SCRIPTS_DIR}/downstream_prepare_jcvi_inputs.py"
NATIVE_INPUT_SCRIPT="${SCRIPTS_DIR}/downstream_make_jcvi_native_plot_inputs.py"

# Override these variables if the server-side filenames differ.
PUNI_GENOME="${PUNI_GENOME:-${ROOT_DIR}/data/puni/GCA_032719365.1_ASM3271936v1_genomic.fna}"
PUNI_GFF="${PUNI_GFF:-${ROOT_DIR}/data/puni/Pun.protein.best.gff}"
PUNI_PEP="${PUNI_PEP:-${ROOT_DIR}/data/puni/Pun.protein.best.fa}"
PUNI_CDS="${PUNI_CDS:-${ROOT_DIR}/data/puni/Pun.cds.best.fa}"

GTR_GENOME="${GTR_GENOME:-${ROOT_DIR}/data/gtr/GCA_018257905.1_V1_genomic.fna}"
GTR_GFF="${GTR_GFF:-${ROOT_DIR}/data/gtr/Gtr.NCBI.gene.gff3}"
GTR_PEP="${GTR_PEP:-${ROOT_DIR}/data/gtr/Gtr.NCBI.prot.faa}"
GTR_CDS="${GTR_CDS:-${ROOT_DIR}/data/gtr/Gtr.NCBI.gene.modId.fna}"

INPUT_DIR="${RUN_DIR}/inputs_extra"
STAGE_DIR="${RUN_DIR}/stage_extra"
PAIR_ROOT="${RUN_DIR}/pairs"
FIG_DIR="${RUN_DIR}/figures"
SUMMARY_DIR="${RUN_DIR}/summary"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }

for f in "${PUNI_GENOME}" "${PUNI_GFF}" "${PUNI_PEP}" "${GTR_GENOME}" "${GTR_GFF}" "${GTR_PEP}"; do
    require_file "${f}"
done

mkdir -p "${INPUT_DIR}" "${STAGE_DIR}/canonical" "${PAIR_ROOT}" "${FIG_DIR}" "${SUMMARY_DIR}"

stage_label() {
    local label="$1" genome="$2" gff="$3" pep="$4" cds="$5"
    local d="${STAGE_DIR}/canonical/${label}"
    mkdir -p "${d}"
    ln -sfn "${genome}" "${d}/genome.fa"
    ln -sfn "${gff}" "${d}/primary.gff3"
    ln -sfn "${pep}" "${d}/primary.pep.fa"
    if [[ -s "${cds}" ]]; then
        ln -sfn "${cds}" "${d}/primary.cds.fa"
    else
        : > "${d}/primary.cds.fa"
        echo "[WARN] No CDS supplied for ${label}; protein synteny will continue." >&2
    fi
}

stage_label puni "${PUNI_GENOME}" "${PUNI_GFF}" "${PUNI_PEP}" "${PUNI_CDS}"
stage_label gtr  "${GTR_GENOME}"  "${GTR_GFF}"  "${GTR_PEP}"  "${GTR_CDS}"

source "${CONDA_SH}"
conda activate "${JCVI_ENV}"
export MPLBACKEND=Agg
export MPLCONFIGDIR="${RUN_DIR}/mplconfig_extra"
mkdir -p "${MPLCONFIGDIR}"

command -v python >/dev/null 2>&1 || die "python unavailable in ${JCVI_ENV}"
command -v diamond >/dev/null 2>&1 || die "diamond unavailable in ${JCVI_ENV}"
python -c 'import jcvi' >/dev/null 2>&1 || die "jcvi is not importable in ${JCVI_ENV}"
require_file "${PREPARE_SCRIPT}"
require_file "${NATIVE_INPUT_SCRIPT}"

python "${PREPARE_SCRIPT}" \
    --analysis-dir "${STAGE_DIR}" \
    --outdir "${INPUT_DIR}" \
    --summary "${SUMMARY_DIR}/input_summary_extra.tsv" \
    --labels puni gtr

awk -F '\t' 'NR > 1 && $4 == 0 {bad=1; print "zero usable BED/PEP records for " $1 > "/dev/stderr"} END {exit bad}' \
    "${SUMMARY_DIR}/input_summary_extra.tsv" || die "Input normalization produced no usable genes"

QUERY="pzi_copy1"
for REF in puni gtr; do
    WORK="${PAIR_ROOT}/${QUERY}_vs_${REF}"
    mkdir -p "${WORK}"
    for label in "${QUERY}" "${REF}"; do
        if [[ "${label}" == "${QUERY}" ]]; then
            src="${RUN_DIR}/inputs/${label}"
        else
            src="${INPUT_DIR}/${label}"
        fi
        for ext in bed pep; do
            require_file "${src}/${label}.${ext}"
            ln -sfn "${src}/${label}.${ext}" "${WORK}/${label}.${ext}"
        done
        if [[ -s "${src}/${label}.cds" ]]; then
            ln -sfn "${src}/${label}.cds" "${WORK}/${label}.cds"
        fi
    done
    cd "${WORK}"

    RAW="${QUERY}.${REF}.anchors"
    SCREENED="${QUERY}.${REF}.anchors.new"
    SIMPLE="${QUERY}.${REF}.anchors.simple"
    if [[ ! -s "${RAW}" ]]; then
        python -m jcvi.compara.catalog ortholog "${QUERY}" "${REF}" \
            --dbtype=prot --align_soft=diamond_blastp \
            --cscore="${CSCORE}" --min_size="${MIN_ANCHORS}" \
            --cpus="${THREADS}" --no_strip_names --no_dotplot
    fi
    require_file "${RAW}"
    if [[ ! -s "${SCREENED}" ]]; then
        python -m jcvi.compara.synteny screen \
            --qbed="${QUERY}.bed" --sbed="${REF}.bed" \
            --minspan="${MINSPAN}" --simple "${RAW}" "${SCREENED}"
    fi
    require_file "${SCREENED}"
    if [[ ! -s "${SIMPLE}" ]]; then
        [[ -s "${SCREENED}.simple" ]] || die "Missing simple anchors for ${QUERY} vs ${REF}"
        ln -sfn "${SCREENED}.simple" "${SIMPLE}"
    fi

    python -m jcvi.graphics.dotplot "${RAW}" \
        --qbed="${QUERY}.bed" --sbed="${REF}.bed" \
        --figsize=12x12 --dpi=300 --format=png --notex \
        --outfile="${FIG_DIR}/${QUERY}_vs_${REF}.dotplot.png"

    python "${NATIVE_INPUT_SCRIPT}" \
        --query "${QUERY}" --ref "${REF}" \
        --query-seqids "${RUN_DIR}/inputs/${QUERY}/${QUERY}.seqids" \
        --ref-seqids "${INPUT_DIR}/${REF}/${REF}.seqids" \
        --query-bed "${QUERY}.bed" --ref-bed "${REF}.bed" \
        --anchors-simple "${SIMPLE}" \
        --out-seqids "${FIG_DIR}/${QUERY}_vs_${REF}.seqids" \
        --out-layout "${FIG_DIR}/${QUERY}_vs_${REF}.layout"
    python -m jcvi.graphics.karyotype \
        "${FIG_DIR}/${QUERY}_vs_${REF}.seqids" \
        "${FIG_DIR}/${QUERY}_vs_${REF}.layout" \
        --format=png --notex \
        --outfile="${FIG_DIR}/${QUERY}_vs_${REF}.karyotype.png"
    cd "${ROOT_DIR}"
done

printf 'comparison\tstatus\tdotplot\tkaryotype\traw_anchors\tscreened_anchors\n' > "${SUMMARY_DIR}/pzi_copy1_additional_pair_status.tsv"
for REF in puni gtr; do
    WORK="${PAIR_ROOT}/${QUERY}_vs_${REF}"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${QUERY}_vs_${REF}" ok \
        "${FIG_DIR}/${QUERY}_vs_${REF}.dotplot.png" \
        "${FIG_DIR}/${QUERY}_vs_${REF}.karyotype.png" \
        "${WORK}/${QUERY}.${REF}.anchors" \
        "${WORK}/${QUERY}.${REF}.anchors.new" \
        >> "${SUMMARY_DIR}/pzi_copy1_additional_pair_status.tsv"
done

touch "${RUN_DIR}/PZI_COPY1_PUNI_GTR_COMPLETE"
echo "[OK] Added pzi_copy1 vs puni and pzi_copy1 vs gtr JCVI dotplots/karyotypes."
