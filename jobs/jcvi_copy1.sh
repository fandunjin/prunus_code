#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

# Gene-level JCVI synteny and native JCVI figures for pzi_copy1 and phased
# haplotypes. Each run has a unique directory and never overwrites old output.
ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/downstream_analysis/karyotype}"
# Four phased PZI inputs are taken from the chromosome-renamed package.  Other
# species and pzi_copy1 continue to use their canonical/source locations.
HAP_RENAME_DIR="${HAP_RENAME_DIR:-${ROOT_DIR}/haploid_rename}"
HAP_GENOME_DIR="${HAP_RENAME_DIR}/genome"
HAP_GFF_DIR="${HAP_RENAME_DIR}/gff"
HAP_CDS_DIR="${HAP_RENAME_DIR}/cds"
HAP_PEP_DIR="${HAP_RENAME_DIR}/pep"
HAP_PACKAGE_CHECKSUMS="${HAP_RENAME_DIR}/checksums.md5"
RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
if [[ -n "${RUN_DIR:-}" ]]; then
    :
elif [[ "${RESUME:-0}" == "1" ]]; then
    latest_resume="$(ls -dt "${OUT_ROOT}"/run_* 2>/dev/null | head -n 1 || true)"
    [[ -n "${latest_resume}" ]] || { echo "[ERROR] RESUME=1 but no run_* directory exists under ${OUT_ROOT}" >&2; exit 1; }
    RUN_DIR="${latest_resume}"
else
    RUN_DIR="${OUT_ROOT}/run_${RUN_ID}"
fi
STAGE_ANALYSIS="${RUN_DIR}/stage"
INPUT_DIR="${RUN_DIR}/inputs"
PAIR_DIR="${RUN_DIR}/pairs"
FIG_DIR="${RUN_DIR}/figures"
SUMMARY_DIR="${RUN_DIR}/summary"
LOG_DIR="${ROOT_DIR}/downstream_analysis/logs"

CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
JCVI_ENV="${JCVI_ENV:-jcvi}"
THREADS="${THREADS:-32}"
CSCORE="${CSCORE:-0.70}"
MIN_ANCHORS="${MIN_ANCHORS:-4}"
MINSPAN="${MINSPAN:-30}"
FORCE="${FORCE:-0}"
HAP_INPUTS_CHANGED=0
MPLBACKEND="${MPLBACKEND:-Agg}"
MULTI_KARYO_FIGSIZE="${MULTI_KARYO_FIGSIZE:-14x24}"

PREPARE_SCRIPT="${ROOT_DIR}/scripts/downstream_prepare_jcvi_inputs.py"
NATIVE_INPUT_SCRIPT="${ROOT_DIR}/scripts/downstream_make_jcvi_native_plot_inputs.py"

# Current 17-species tree. pzi_copy1 and the four phased haplotypes are queries.
QUERY_LABELS=(pzi_copy1 pzi_hapA1 pzi_hapA2 pzi_hapB1 pzi_hapB2)
REF_LABELS=(fve rch pcom msyl mdom pcam pdu pmu_refseq ppe pspe psp ppa pmal pav psal par)

die() { echo "[ERROR] $*" >&2; exit 1; }
warn() { echo "[WARN] $*" >&2; }
require_file() { [[ -s "$1" ]] || die "$2: $1"; }

is_complete_triplet() {
    local d="$1"
    [[ -s "${d}/primary.gff3" && -s "${d}/primary.pep.fa" && -s "${d}/primary.cds.fa" ]]
}

# Set SRC_GENOME/SRC_GFF/SRC_PEP/SRC_CDS for one label.
resolve_sources() {
    local label="$1"
    local canonical="${ROOT_DIR}/downstream_analysis/work/canonical/${label}"
    SRC_GENOME=""; SRC_GFF=""; SRC_PEP=""; SRC_CDS=""

    # The renamed haplotype package is authoritative.  Check it before the
    # reusable canonical cache so an older staged chr01-style input cannot
    # silently override the current C01--C08 chromosome names.
    case "${label}" in
        pzi_hapA1|pzi_hapA2|pzi_hapB1|pzi_hapB2)
            SRC_GENOME="${HAP_GENOME_DIR}/${label}.fa"
            SRC_GFF="${HAP_GFF_DIR}/${label}.gff3"
            SRC_PEP="${HAP_PEP_DIR}/${label}.pep.fa"
            SRC_CDS="${HAP_CDS_DIR}/${label}.cds.fa"
            return 0
            ;;
    esac

    if is_complete_triplet "${canonical}"; then
        SRC_GENOME="${canonical}/genome.fa"
        SRC_GFF="${canonical}/primary.gff3"
        SRC_PEP="${canonical}/primary.pep.fa"
        SRC_CDS="${canonical}/primary.cds.fa"
        return 0
    fi

    case "${label}" in
        pzi_copy1)
            SRC_GENOME="${ROOT_DIR}/final_genome/pzi_copy1.fa"
            SRC_GFF="${ROOT_DIR}/final_genome/pzi_copy1.gff"
            SRC_PEP="${ROOT_DIR}/final_genome/pzi_copy1.pep"
            SRC_CDS="${ROOT_DIR}/final_genome/pzi_copy1.cds"
            ;;
        *)
            local data_label="${label}"
            [[ "${label}" == "pmu_refseq" ]] && data_label="pmu"
            local data="${ROOT_DIR}/data/${data_label}"
            SRC_GENOME="$(find "${data}/annotation/eviann/postprocess" -type f \( -name "${label}.genomic.fa" -o -name "${label}.fa" \) -print -quit 2>/dev/null || true)"
            [[ -n "${SRC_GENOME}" ]] || SRC_GENOME="$(find "${data}/annotation/helixer/work" -type f \( -name "${label}.genomic.fa" -o -name "${label}.fa" \) -print -quit 2>/dev/null || true)"
            [[ -n "${SRC_GENOME}" ]] || SRC_GENOME="$(find "${data}" -type f \( -name '*genomic.fna' -o -name '*genomic.fa' -o -name '*genomic.fna.gz' -o -name '*genomic.fa.gz' \) ! -name '*cds_from_genomic*' ! -name '*rna_from_genomic*' -print -quit 2>/dev/null || true)"

            SRC_GFF="$(find "${data}/annotation/eviann/postprocess" -type f \( -name "${label}.eviann.genes.gff3" -o -name '*.genes.gff3' -o -name '*.gff3' \) -print -quit 2>/dev/null || true)"
            SRC_PEP="$(find "${data}/annotation/eviann/postprocess" -type f \( -name "${label}.eviann.protein.fa" -o -name '*.protein.fa' -o -name '*.protein.faa' \) -print -quit 2>/dev/null || true)"
            SRC_CDS="$(find "${data}/annotation/eviann/postprocess" -type f \( -name "${label}.eviann.cds.fa" -o -name '*.cds.fa' \) -print -quit 2>/dev/null || true)"

            [[ -n "${SRC_GFF}" ]] || SRC_GFF="$(find "${data}" -type f \( -name '*genomic.gff3' -o -name '*genomic.gff3.gz' -o -name '*genomic.gff' -o -name '*genomic.gff.gz' -o -name '*.gff3' -o -name '*.gff3.gz' -o -name '*.gff' -o -name '*.gff.gz' \) -print -quit 2>/dev/null || true)"
            [[ -n "${SRC_PEP}" ]] || SRC_PEP="$(find "${data}" -type f \( -name '*protein.faa' -o -name '*protein.faa.gz' -o -name '*protein.fa' -o -name '*protein.fa.gz' \) -print -quit 2>/dev/null || true)"
            [[ -n "${SRC_CDS}" ]] || SRC_CDS="$(find "${data}" -type f \( -name '*cds_from_genomic.fna' -o -name '*cds_from_genomic.fna.gz' -o -name '*.cds.fa' -o -name '*.cds.fa.gz' \) -print -quit 2>/dev/null || true)"
            ;;
    esac
}

materialize_or_link() {
    local src="$1" dst="$2"
    [[ -n "${src}" && -e "${src}" ]] || return 0
    mkdir -p "$(dirname "${dst}")"
    rm -f "${dst}" "${dst}.tmp"
    if [[ "${src}" == *.gz ]]; then
        gzip -cd "${src}" > "${dst}.tmp"
        mv "${dst}.tmp" "${dst}"
    else
        ln -s "${src}" "${dst}"
    fi
}

stage_label() {
    local label="$1" d
    resolve_sources "${label}"
    d="${STAGE_ANALYSIS}/canonical/${label}"
    mkdir -p "${d}"

    if [[ ! -s "${SRC_GFF}" || ! -s "${SRC_PEP}" ]]; then
        die "Missing required JCVI inputs for ${label}; GFF=${SRC_GFF:-NA}; PEP=${SRC_PEP:-NA}"
    fi

    materialize_or_link "${SRC_GENOME}" "${d}/genome.fa"
    materialize_or_link "${SRC_GFF}" "${d}/primary.gff3"
    materialize_or_link "${SRC_PEP}" "${d}/primary.pep.fa"
    if [[ -n "${SRC_CDS}" && -s "${SRC_CDS}" ]]; then
        materialize_or_link "${SRC_CDS}" "${d}/primary.cds.fa"
    else
        : > "${d}/primary.cds.fa"
        warn "No CDS found for ${label}; JCVI protein synteny will continue without CDS."
    fi
    printf '%s\t%s\t%s\t%s\t%s\n' "${label}" "${SRC_GENOME:-NA}" "${SRC_GFF}" "${SRC_PEP}" "${SRC_CDS:-NA}" >> "${SUMMARY_DIR}/source_manifest.tsv"
    return 0
}

run_dotplot() {
    local q="$1" r="$2" out="$3"
    python -m jcvi.graphics.dotplot "${q}.${r}.anchors" \
        --qbed="${q}.bed" --sbed="${r}.bed" \
        --figsize=12x12 --dpi=300 --format=png --notex --outfile="${out}"
}

run_karyotype() {
    local seqids="$1" layout="$2" out="$3"
    python -m jcvi.graphics.karyotype "${seqids}" "${layout}" \
        --format=png --notex --outfile="${out}"
}

valid_simple_file() {
    local simple="$1"
    [[ -s "${simple}" ]] || return 1
    awk 'NF && $1 !~ /^#/ {if (NF != 6) bad=1; else seen=1} END {exit (!seen || bad)}' "${simple}"
}

mkdir -p "${RUN_DIR}" "${SUMMARY_DIR}" "${PAIR_DIR}" "${FIG_DIR}" "${LOG_DIR}"
for label in pzi_hapA1 pzi_hapA2 pzi_hapB1 pzi_hapB2; do
    require_file "${HAP_GENOME_DIR}/${label}.fa" "renamed ${label} genome"
    require_file "${HAP_GFF_DIR}/${label}.gff3" "renamed ${label} GFF3"
    require_file "${HAP_CDS_DIR}/${label}.cds.fa" "renamed ${label} CDS"
    require_file "${HAP_PEP_DIR}/${label}.pep.fa" "renamed ${label} protein"
done
require_file "${HAP_PACKAGE_CHECKSUMS}" "renamed haplotype package checksums"
SAVED_HAP_PACKAGE_CHECKSUMS="${SUMMARY_DIR}/haploid_rename.checksums.md5"
if ! cmp -s "${HAP_PACKAGE_CHECKSUMS}" "${SAVED_HAP_PACKAGE_CHECKSUMS}"; then
    HAP_INPUTS_CHANGED=1
    echo "[INFO] Renamed haplotype package differs from this run's staged inputs."
fi
if [[ "${FORCE}" == "1" || "${HAP_INPUTS_CHANGED}" == "1" || ! -s "${SUMMARY_DIR}/source_manifest.tsv" ]]; then
    printf 'label\tgenome\tgff\tpep\tcds\n' > "${SUMMARY_DIR}/source_manifest.tsv"
fi
printf 'comparison\tquery\tref\tstatus\traw_anchors\tscreened_anchors\tdotplot\tkaryotype\n' > "${SUMMARY_DIR}/pair_status.tsv"

[[ -s "${CONDA_SH}" ]] || die "Conda initialization script not found: ${CONDA_SH}"
source "${CONDA_SH}"
conda activate "${JCVI_ENV}"
export MPLBACKEND MPLCONFIGDIR="${RUN_DIR}/mplconfig"
mkdir -p "${MPLCONFIGDIR}"

command -v python >/dev/null 2>&1 || die "python is unavailable in ${JCVI_ENV}"
command -v diamond >/dev/null 2>&1 || die "diamond is unavailable in the activated ${JCVI_ENV} environment"
echo "[INFO] Using DIAMOND from activated ${JCVI_ENV}: $(command -v diamond)"
diamond version >/dev/null 2>&1 || die "DIAMOND was found but cannot be executed: $(command -v diamond)"
python -c 'import jcvi' >/dev/null 2>&1 || die "jcvi is not importable in ${JCVI_ENV}"
python -c 'from ete4 import Tree' >/dev/null 2>&1 || die "ete4 is not importable in ${JCVI_ENV}; install ete4 before running JCVI karyotype"
require_file "${PREPARE_SCRIPT}" "JCVI input preparation script"
require_file "${NATIVE_INPUT_SCRIPT}" "JCVI native plot input script"

echo "=== JCVI pzi_copy1/haplotype synteny ==="
echo "Output: ${RUN_DIR}"
echo "Queries: ${QUERY_LABELS[*]}"
echo "References: ${REF_LABELS[*]}"
echo "Resume mode: ${RESUME:-0}"
if [[ -d "${PAIR_DIR}" ]]; then
    echo "Existing raw anchors: $(find "${PAIR_DIR}" -type f -name '*.anchors' -size +0c | wc -l)"
    echo "Existing simple anchors: $(find "${PAIR_DIR}" -type f -name '*.anchors.simple' -size +0c | wc -l)"
    echo "Existing dotplots: $(find "${FIG_DIR}" -type f -name '*.dotplot.png' -size +0c | wc -l)"
    echo "Existing karyotypes: $(find "${FIG_DIR}" -type f -name '*.karyotype.png' -size +0c | wc -l)"
fi

inputs_ready=1
for label in "${QUERY_LABELS[@]}" "${REF_LABELS[@]}"; do
    for required in "${INPUT_DIR}/${label}/${label}.bed" "${INPUT_DIR}/${label}/${label}.pep" "${INPUT_DIR}/${label}/${label}.seqids"; do
        [[ -s "${required}" ]] || inputs_ready=0
    done
done

if [[ "${FORCE}" != "1" && "${HAP_INPUTS_CHANGED}" != "1" && "${inputs_ready}" == "1" && -s "${SUMMARY_DIR}/input_summary.tsv" ]]; then
    echo "[INFO] Existing normalized JCVI inputs found; skip staging and input preparation."
else
    for label in "${QUERY_LABELS[@]}" "${REF_LABELS[@]}"; do
        stage_label "${label}"
    done

    python "${PREPARE_SCRIPT}" \
        --analysis-dir "${STAGE_ANALYSIS}" \
        --outdir "${INPUT_DIR}" \
        --summary "${SUMMARY_DIR}/input_summary.tsv" \
        --labels "${QUERY_LABELS[@]}" "${REF_LABELS[@]}"
fi

# Fail early if a GFF/FASTA identifier convention produced no usable genes.
awk -F '\t' 'NR > 1 && $4 == 0 {print "zero usable BED/PEP records for " $1 > "/dev/stderr"; bad=1} END {exit bad}' \
    "${SUMMARY_DIR}/input_summary.tsv" || die "JCVI input normalization produced an empty species input; inspect input_summary.tsv"

for query in "${QUERY_LABELS[@]}"; do
    for ref in "${REF_LABELS[@]}"; do
        comparison="${query}_vs_${ref}"
        work="${PAIR_DIR}/${comparison}"
        pair_inputs_changed=0
        if [[ "${HAP_INPUTS_CHANGED}" == "1" && "${query}" == pzi_hap* ]]; then
            pair_inputs_changed=1
        fi
        mkdir -p "${work}"
        for label in "${query}" "${ref}"; do
            ln -sfn "${INPUT_DIR}/${label}/${label}.bed" "${work}/${label}.bed"
            ln -sfn "${INPUT_DIR}/${label}/${label}.pep" "${work}/${label}.pep"
            ln -sfn "${INPUT_DIR}/${label}/${label}.cds" "${work}/${label}.cds"
        done

        cd "${work}"
        if [[ "${FORCE}" == "1" || ! -s "${query}.${ref}.anchors" ]]; then
            echo "[INFO] JCVI catalog ortholog: ${comparison}"
            python -m jcvi.compara.catalog ortholog "${query}" "${ref}" \
                --dbtype=prot --align_soft=diamond_blastp \
                --cscore="${CSCORE}" --min_size="${MIN_ANCHORS}" \
                --cpus="${THREADS}" --no_strip_names --no_dotplot
        fi

        raw="${work}/${query}.${ref}.anchors"
        screened="${work}/${query}.${ref}.anchors.new"
        simple="${work}/${query}.${ref}.anchors.simple"
        if [[ ! -s "${raw}" ]]; then
            printf '%s\t%s\t%s\tno_raw_anchors\t%s\t%s\tNA\tNA\n' "${comparison}" "${query}" "${ref}" "${raw}" "${screened}" >> "${SUMMARY_DIR}/pair_status.tsv"
            cd "${ROOT_DIR}"
            continue
        fi

        if [[ "${pair_inputs_changed}" == "1" ]]; then
            rm -f "${screened}" "${simple}" "${screened}.simple"
        fi
        if [[ "${FORCE}" == "1" || "${pair_inputs_changed}" == "1" || ! -s "${screened}" || ! -s "${simple}" ]]; then
            python -m jcvi.compara.synteny screen \
                --qbed="${query}.bed" --sbed="${ref}.bed" \
                --minspan="${MINSPAN}" --simple "${raw}" "${screened}"
        fi

        if [[ ! -s "${screened}" ]]; then
            printf '%s\t%s\t%s\tno_screened_anchors\t%s\t%s\tNA\tNA\n' "${comparison}" "${query}" "${ref}" "${raw}" "${screened}" >> "${SUMMARY_DIR}/pair_status.tsv"
            cd "${ROOT_DIR}"
            continue
        fi

        # JCVI writes the six-column macro-synteny blocks required by
        # graphics.karyotype as <query>.<ref>.anchors.simple.  The screened
        # .anchors.new file contains three-column gene pairs and is only for
        # dotplots/statistics.
        if ! valid_simple_file "${simple}"; then
            alt_simple="${screened}.simple"
            if valid_simple_file "${alt_simple}"; then
                ln -sfn "${alt_simple}" "${simple}"
            else
                die "Missing/invalid six-column JCVI simple file for ${comparison}: ${simple}"
            fi
        fi

        seqids="${FIG_DIR}/${comparison}.seqids"
        layout="${FIG_DIR}/${comparison}.layout"
        dotplot="${FIG_DIR}/${comparison}.dotplot.png"
        karyotype="${FIG_DIR}/${comparison}.karyotype.png"
        python "${NATIVE_INPUT_SCRIPT}" \
            --query "${query}" --ref "${ref}" \
            --query-seqids "${INPUT_DIR}/${query}/${query}.seqids" \
            --ref-seqids "${INPUT_DIR}/${ref}/${ref}.seqids" \
            --query-bed "${query}.bed" --ref-bed "${ref}.bed" \
            --anchors-simple "${simple}" \
            --out-seqids "${seqids}" --out-layout "${layout}"
        if [[ "${FORCE}" == "1" || "${pair_inputs_changed}" == "1" || ! -s "${dotplot}" ]]; then
            run_dotplot "${query}" "${ref}" "${dotplot}"
        else
            echo "[INFO] Existing dotplot found; skip: ${dotplot}"
        fi
        if [[ "${FORCE}" == "1" || "${pair_inputs_changed}" == "1" || ! -s "${karyotype}" ]]; then
            run_karyotype "${seqids}" "${layout}" "${karyotype}"
        else
            echo "[INFO] Existing karyotype found; skip: ${karyotype}"
        fi
        printf '%s\t%s\t%s\tok\t%s\t%s\t%s\t%s\n' "${comparison}" "${query}" "${ref}" "${raw}" "${screened}" "${dotplot}" "${karyotype}" >> "${SUMMARY_DIR}/pair_status.tsv"
        cd "${ROOT_DIR}"
    done
done

# Build one multi-species ribbon plot.  JCVI karyotype uses one seqids line
# and one layout row per species; each edge below links pzi_copy1 to a
# reference species using its screened JCVI anchors.
MULTI_SEQIDS="${FIG_DIR}/ribbon_17_species.seqids"
MULTI_LAYOUT="${FIG_DIR}/ribbon_17_species.layout"
MULTI_KARYOTYPE="${FIG_DIR}/ribbon_17_species.karyotype.png"
MULTI_SPECIES=(pzi_copy1 "${REF_LABELS[@]}")
MULTI_Y=(0.9600 0.9050 0.8500 0.7950 0.7400 0.6850 0.6300 0.5750 0.5200 0.4650 0.4100 0.3550 0.3000 0.2450 0.1900 0.1350 0.0800)

: > "${MULTI_SEQIDS}"
for label in "${MULTI_SPECIES[@]}"; do
    cat "${INPUT_DIR}/${label}/${label}.seqids" >> "${MULTI_SEQIDS}"
    ln -sfn "${INPUT_DIR}/${label}/${label}.bed" "${FIG_DIR}/${label}.bed"
done

{
    echo "# y, xstart, xend, rotation, color, label, va, bed"
    total="${#MULTI_SPECIES[@]}"
    for ((i=0; i<total; i++)); do
        label="${MULTI_SPECIES[$i]}"
        y="${MULTI_Y[$i]}"
        if [[ "${i}" -eq 0 ]]; then
            color="#4477AA"
        else
            color="#CC6677"
        fi
        printf ' %s, .08, .92, 0, %s, %s, center, %s.bed\n' "${y}" "${color}" "${label}" "${label}"
    done
    echo "# edges"
    edge_count=0
    for ((i=1; i<total; i++)); do
        ref="${MULTI_SPECIES[$i]}"
        anchor="${PAIR_DIR}/pzi_copy1_vs_${ref}/pzi_copy1.${ref}.anchors.simple"
        if [[ -s "${anchor}" ]]; then
            ln -sfn "${anchor}" "${FIG_DIR}/pzi_copy1.${ref}.anchors.simple"
            printf ' e, 0, %d, pzi_copy1.%s.anchors.simple\n' "${i}" "${ref}"
            edge_count=$((edge_count + 1))
        else
            warn "No screened anchors for pzi_copy1 vs ${ref}; ${ref} will remain unconnected in the 17-species ribbon."
        fi
    done
    (( edge_count > 0 )) || die "No pzi_copy1 reference anchors are available for the multi-species ribbon."
} > "${MULTI_LAYOUT}"

cd "${FIG_DIR}"
if [[ "${FORCE}" == "1" || ! -s "${MULTI_KARYOTYPE}" ]]; then
    python -m jcvi.graphics.karyotype \
        "${MULTI_SEQIDS}" "${MULTI_LAYOUT}" \
        --figsize="${MULTI_KARYO_FIGSIZE}" \
        --format=png --notex --outfile="${MULTI_KARYOTYPE}"
else
    echo "[INFO] Existing 17-species ribbon found; skip: ${MULTI_KARYOTYPE}"
fi
cd "${ROOT_DIR}"

cp -f "${HAP_PACKAGE_CHECKSUMS}" "${SAVED_HAP_PACKAGE_CHECKSUMS}"

printf '%s\n' "${RUN_DIR}" > "${OUT_ROOT}/latest.txt"
ln -sfn "${RUN_DIR}" "${OUT_ROOT}/latest"
cat > "${RUN_DIR}/README.txt" <<EOF
JCVI synteny run for pzi_copy1 and pzi hapA1/A2/B1/B2.

Queries: ${QUERY_LABELS[*]}
References: ${REF_LABELS[*]}
JCVI parameters: cscore=${CSCORE}, min_size=${MIN_ANCHORS}, minspan=${MINSPAN}, cpus=${THREADS}

Native JCVI figures:
  figures/*.dotplot.png
  figures/*.karyotype.png
  figures/ribbon_17_species.karyotype.png

Input/anchor tables:
  inputs/<label>/*.bed,*.pep,*.cds,*.seqids
  pairs/<query>_vs_<ref>/*.anchors*
  summary/source_manifest.tsv
  summary/input_summary.tsv
  summary/pair_status.tsv
EOF

echo "[OK] JCVI synteny and figures completed: ${RUN_DIR}"
