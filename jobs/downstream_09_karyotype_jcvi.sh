#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/analysis}"
SCRIPTS_DIR="${SCRIPTS_DIR:-${ROOT_DIR}/scripts}"
HAP_RENAME_DIR="${HAP_RENAME_DIR:-${ROOT_DIR}/haploid_rename}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
SYNTENY_ENV="${SYNTENY_ENV:-synteny}"
THREADS="${THREADS:-32}"

KARYO_DIR="${ANALYSIS_DIR}/karyotype"
JCVI_DIR="${KARYO_DIR}/jcvi"
INPUT_DIR="${JCVI_DIR}/inputs"
PAIR_DIR="${JCVI_DIR}/pairs"
SUMMARY_DIR="${JCVI_DIR}/summary"

# JCVI defaults are intentionally not too strict here. PZI A2B2 is expected to
# show duplicated 2:1 synteny against diploid Prunus references, so cscore=.99
# would bias the result toward one reciprocal-best copy.
CSCORE="${CSCORE:-0.70}"
MIN_ANCHORS="${MIN_ANCHORS:-4}"
MINSPAN="${MINSPAN:-30}"
RUN_GENOME_SCREEN="${RUN_GENOME_SCREEN:-0}"
PAF_PRESET="${PAF_PRESET:-asm5}"

QUERY_LABELS=(pzi_A2B2 pzi_hapA1 pzi_hapA2 pzi_hapB1 pzi_hapB2)
REF_LABELS=(ppe pdu pmu_refseq pcam pspe rch fve)
INTERNAL_HAP_PAIRS=(
    pzi_hapA1:pzi_hapA2
    pzi_hapB1:pzi_hapB2
    pzi_hapA1:pzi_hapB1
    pzi_hapA1:pzi_hapB2
    pzi_hapA2:pzi_hapB1
    pzi_hapA2:pzi_hapB2
)
ALL_JCVI_LABELS=(pzi_A2B2 pzi_hapA1 pzi_hapA2 pzi_hapB1 pzi_hapB2 ppe pdu pmu_refseq pcam pspe rch fve)
GENOME_SCREEN_LABELS=(ppa pav par psal pdav pmir pseu pserru)

stage_renamed_haplotypes() {
    local label canonical
    for label in pzi_hapA1 pzi_hapA2 pzi_hapB1 pzi_hapB2; do
        canonical="${ANALYSIS_DIR}/canonical/${label}"
        mkdir -p "${canonical}"
        for source in \
            "${HAP_RENAME_DIR}/genome/${label}.fa" \
            "${HAP_RENAME_DIR}/gff/${label}.gff3" \
            "${HAP_RENAME_DIR}/cds/${label}.cds.fa" \
            "${HAP_RENAME_DIR}/pep/${label}.pep.fa"; do
            if [[ ! -s "${source}" ]]; then
                echo "[ERROR] Missing renamed haplotype input: ${source}" >&2
                exit 2
            fi
        done
        ln -sfn "${HAP_RENAME_DIR}/genome/${label}.fa" "${canonical}/genome.fa"
        ln -sfn "${HAP_RENAME_DIR}/gff/${label}.gff3" "${canonical}/primary.gff3"
        ln -sfn "${HAP_RENAME_DIR}/cds/${label}.cds.fa" "${canonical}/primary.cds.fa"
        ln -sfn "${HAP_RENAME_DIR}/pep/${label}.pep.fa" "${canonical}/primary.pep.fa"
    done
}

mkdir -p "${INPUT_DIR}" "${PAIR_DIR}" "${SUMMARY_DIR}" "${ANALYSIS_DIR}/logs"
stage_renamed_haplotypes

set +u
source "${CONDA_SH}"
conda activate "${SYNTENY_ENV}"
set -u

for exe in python diamond minimap2; do
    if ! command -v "${exe}" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found in SYNTENY_ENV=${SYNTENY_ENV}: ${exe}" >&2
        exit 2
    fi
done

python -c "import jcvi" >/dev/null 2>&1 || {
    echo "[ERROR] Python package jcvi is not importable in SYNTENY_ENV=${SYNTENY_ENV}" >&2
    exit 2
}

for script in downstream_prepare_jcvi_inputs.py downstream_summarize_jcvi_pairs.py; do
    if [ ! -s "${SCRIPTS_DIR}/${script}" ]; then
        echo "[ERROR] Missing helper script: ${SCRIPTS_DIR}/${script}" >&2
        exit 2
    fi
done

echo "[INFO] Preparing JCVI BED/PEP/CDS inputs"
echo "[INFO] Four-haplotype package: ${HAP_RENAME_DIR}"
python "${SCRIPTS_DIR}/downstream_prepare_jcvi_inputs.py" \
    --analysis-dir "${ANALYSIS_DIR}" \
    --outdir "${INPUT_DIR}" \
    --summary "${SUMMARY_DIR}/jcvi_input_summary.tsv" \
    --labels "${ALL_JCVI_LABELS[@]}"

printf "query_label\tref_label\tcomparison_type\n" > "${SUMMARY_DIR}/jcvi_pair_plan.tsv"
for query in "${QUERY_LABELS[@]}"; do
    for ref in "${REF_LABELS[@]}"; do
        printf "%s\t%s\texternal_reference\n" "${query}" "${ref}" >> "${SUMMARY_DIR}/jcvi_pair_plan.tsv"
    done
done
for pair in "${INTERNAL_HAP_PAIRS[@]}"; do
    query="${pair%%:*}"
    ref="${pair#*:}"
    printf "%s\t%s\tinternal_haplotype\n" "${query}" "${ref}" >> "${SUMMARY_DIR}/jcvi_pair_plan.tsv"
done

printf "comparison\tquery_label\tref_label\tcomparison_type\tstatus\tn_blocks\tn_syntenic_pairs\tscreened_anchors\tjcvi_toaligns\n" > "${SUMMARY_DIR}/jcvi_pair_status.tsv"
printf "comparison\tcomparison_type\tblock_id\tquery_label\tref_label\tquery_gene\tref_gene\tquery_seqid\tquery_start0\tquery_end\tquery_strand\tref_seqid\tref_start0\tref_end\tref_strand\n" > "${SUMMARY_DIR}/jcvi_syntenic_pairs.tsv"
printf "comparison\tcomparison_type\tblock_id\tquery_label\tref_label\tn_pairs\tquery_seqids\tquery_start0_min\tquery_end_max\tref_seqids\tref_start0_min\tref_end_max\n" > "${SUMMARY_DIR}/jcvi_block_summary.tsv"

while IFS=$'\t' read -r query ref comparison_type; do
    if [ "${query}" = "query_label" ]; then
        continue
    fi
        comparison="${query}_vs_${ref}"
        pair_path="${PAIR_DIR}/${comparison}"
        mkdir -p "${pair_path}"

        ln -sf "${INPUT_DIR}/${query}/${query}.bed" "${pair_path}/${query}.bed"
        ln -sf "${INPUT_DIR}/${query}/${query}.pep" "${pair_path}/${query}.pep"
        ln -sf "${INPUT_DIR}/${query}/${query}.cds" "${pair_path}/${query}.cds"
        ln -sf "${INPUT_DIR}/${ref}/${ref}.bed" "${pair_path}/${ref}.bed"
        ln -sf "${INPUT_DIR}/${ref}/${ref}.pep" "${pair_path}/${ref}.pep"
        ln -sf "${INPUT_DIR}/${ref}/${ref}.cds" "${pair_path}/${ref}.cds"

        echo "[INFO] JCVI ortholog and synteny: ${comparison}"
        cd "${pair_path}"

        python -m jcvi.compara.catalog ortholog "${query}" "${ref}" \
            --dbtype=prot \
            --align_soft=diamond_blastp \
            --cscore="${CSCORE}" \
            --min_size="${MIN_ANCHORS}" \
            --cpus="${THREADS}" \
            --no_strip_names \
            --no_dotplot

        anchor_file="${query}.${ref}.anchors"
        screened_anchor="${query}.${ref}.anchors.new"
        lifted_anchor="${query}.${ref}.lifted.anchors"
        lifted_screened_anchor="${query}.${ref}.lifted.anchors.new"
        toaligns_file="${comparison}.jcvi_toaligns.tsv"

        if [ ! -s "${anchor_file}" ]; then
            printf "comparison\tquery_label\tref_label\tcomparison_type\tstatus\tn_blocks\tn_syntenic_pairs\tscreened_anchors\tjcvi_toaligns\n" > "${comparison}.status.tsv"
            printf "%s\t%s\t%s\t%s\t%s\t0\t0\t%s\t%s\n" "${comparison}" "${query}" "${ref}" "${comparison_type}" "no_raw_anchors" "${screened_anchor}" "${toaligns_file}" >> "${comparison}.status.tsv"
            tail -n +2 "${comparison}.status.tsv" >> "${SUMMARY_DIR}/jcvi_pair_status.tsv"
            cd "${ROOT_DIR}"
            continue
        fi

        python -m jcvi.compara.synteny screen \
            --qbed="${query}.bed" \
            --sbed="${ref}.bed" \
            --minspan="${MINSPAN}" \
            --simple \
            "${anchor_file}" \
            "${screened_anchor}"

        if [ -s "${lifted_anchor}" ]; then
            python -m jcvi.compara.synteny screen \
                --qbed="${query}.bed" \
                --sbed="${ref}.bed" \
                --minspan="${MINSPAN}" \
                --simple \
                "${lifted_anchor}" \
                "${lifted_screened_anchor}"
        fi

        if [ -s "${screened_anchor}" ]; then
            python -m jcvi.compara.synteny toaligns \
                --prefix="${comparison}_block" \
                --outfile="${toaligns_file}" \
                "${screened_anchor}"

            python "${SCRIPTS_DIR}/downstream_summarize_jcvi_pairs.py" \
                --comparison "${comparison}" \
                --query "${query}" \
                --ref "${ref}" \
                --comparison-type "${comparison_type}" \
                --aligns "${toaligns_file}" \
                --anchors "${screened_anchor}" \
                --qbed "${query}.bed" \
                --sbed "${ref}.bed" \
                --out-pairs "${comparison}.syntenic_pairs.tsv" \
                --out-blocks "${comparison}.block_summary.tsv" \
                --out-status "${comparison}.status.tsv"
        else
            printf "comparison\tquery_label\tref_label\tcomparison_type\tstatus\tn_blocks\tn_syntenic_pairs\tscreened_anchors\tjcvi_toaligns\n" > "${comparison}.status.tsv"
            printf "%s\t%s\t%s\t%s\t%s\t0\t0\t%s\t%s\n" "${comparison}" "${query}" "${ref}" "${comparison_type}" "no_screened_anchors" "${screened_anchor}" "${toaligns_file}" >> "${comparison}.status.tsv"
            printf "comparison\tcomparison_type\tblock_id\tquery_label\tref_label\tquery_gene\tref_gene\tquery_seqid\tquery_start0\tquery_end\tquery_strand\tref_seqid\tref_start0\tref_end\tref_strand\n" > "${comparison}.syntenic_pairs.tsv"
            printf "comparison\tcomparison_type\tblock_id\tquery_label\tref_label\tn_pairs\tquery_seqids\tquery_start0_min\tquery_end_max\tref_seqids\tref_start0_min\tref_end_max\n" > "${comparison}.block_summary.tsv"
        fi

        tail -n +2 "${comparison}.status.tsv" >> "${SUMMARY_DIR}/jcvi_pair_status.tsv"
        tail -n +2 "${comparison}.syntenic_pairs.tsv" >> "${SUMMARY_DIR}/jcvi_syntenic_pairs.tsv"
        tail -n +2 "${comparison}.block_summary.tsv" >> "${SUMMARY_DIR}/jcvi_block_summary.tsv"

        cd "${ROOT_DIR}"
done < "${SUMMARY_DIR}/jcvi_pair_plan.tsv"

if [ "${RUN_GENOME_SCREEN}" = "1" ]; then
    mkdir -p "${KARYO_DIR}/genome_screen"
    for label in "${GENOME_SCREEN_LABELS[@]}"; do
        ref_genome="${ANALYSIS_DIR}/canonical/pzi_A1A2B1B2/genome.fa"
        qry_genome="${ANALYSIS_DIR}/canonical/${label}/genome.fa"
        if [ -s "${ref_genome}" ] && [ -s "${qry_genome}" ]; then
            minimap2 -t "${THREADS}" -cx "${PAF_PRESET}" "${ref_genome}" "${qry_genome}" > "${KARYO_DIR}/genome_screen/pzi32_vs_${label}.${PAF_PRESET}.paf"
        else
            echo "[WARN] Skip genome-level screen for ${label}: genome.fa missing" >&2
        fi
    done
fi

cat > "${JCVI_DIR}/README.txt" <<'EOF'
downstream_09_karyotype_jcvi.sh generates formal gene-level JCVI synteny for:

Queries:
  pzi_A2B2, pzi_hapA1, pzi_hapA2, pzi_hapB1, pzi_hapB2

References:
  ppe, pdu, pmu_refseq, pcam, pspe, rch, fve

Internal haplotype comparisons:
  pzi_hapA1 vs pzi_hapA2
  pzi_hapB1 vs pzi_hapB2
  pzi_hapA1 vs pzi_hapB1
  pzi_hapA1 vs pzi_hapB2
  pzi_hapA2 vs pzi_hapB1
  pzi_hapA2 vs pzi_hapB2

Main outputs:
  inputs/<label>/<label>.bed/.pep/.cds
  inputs/<label>/<label>.id_map.tsv
  pairs/<query>_vs_<ref>/<query>.<ref>.anchors
  pairs/<query>_vs_<ref>/<query>.<ref>.anchors.new
  pairs/<query>_vs_<ref>/<query>.<ref>.anchors.simple
  pairs/<query>_vs_<ref>/<query>_vs_<ref>.jcvi_toaligns.tsv
  pairs/<query>_vs_<ref>/<query>_vs_<ref>.syntenic_pairs.tsv
  pairs/<query>_vs_<ref>/<query>_vs_<ref>.block_summary.tsv
  summary/jcvi_input_summary.tsv
  summary/jcvi_pair_plan.tsv
  summary/jcvi_pair_status.tsv
  summary/jcvi_syntenic_pairs.tsv
  summary/jcvi_block_summary.tsv

Default parameters:
  JCVI catalog ortholog: --dbtype=prot --align_soft=diamond_blastp --cscore=0.70 --min_size=4 --no_strip_names --no_dotplot
  JCVI synteny screen:  --minspan=30 --simple

The default cscore is intentionally 0.70, not 0.99, because pzi_A2B2 versus
diploid references can contain real 2:1 subgenome synteny. The stricter 0.99
RBH setting is better for one-to-one ortholog blocks, but can drop the second
homeologous copy in a polyploid query.
EOF

echo "[INFO] JCVI gene-level synteny finished: ${JCVI_DIR}"
