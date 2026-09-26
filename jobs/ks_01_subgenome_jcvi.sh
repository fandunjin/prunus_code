#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/downstream_analysis/ks}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
JCVI_ENV="${JCVI_ENV:-jcvi}"
THREADS="${THREADS:-32}"
CSCORE="${CSCORE:-0.70}"
MIN_ANCHORS="${MIN_ANCHORS:-4}"
MINSPAN="${MINSPAN:-30}"
MIN_PEP_MATCH_FRAC="${MIN_PEP_MATCH_FRAC:-0.80}"
MIN_CDS_MATCH_FRAC="${MIN_CDS_MATCH_FRAC:-0.80}"

STAGE_DIR="${KS_ROOT}/stage"
JCVI_DIR="${KS_ROOT}/jcvi"
INPUT_DIR="${JCVI_DIR}/inputs"
PAIR_DIR="${JCVI_DIR}/pairs"
SUMMARY_DIR="${JCVI_DIR}/summary"
FIG_DIR="${JCVI_DIR}/figures"
PREPARE_SCRIPT="${ROOT_DIR}/scripts/downstream_prepare_jcvi_inputs.py"
PAIR_SCRIPT="${ROOT_DIR}/scripts/ks_jcvi_anchors_to_pairs.py"
NATIVE_INPUT_SCRIPT="${ROOT_DIR}/scripts/downstream_make_jcvi_native_plot_inputs.py"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }

if [[ -z "${PZI_RENAME_DIR:-}" ]]; then
    for candidate in \
        "${ROOT_DIR}/data/pzi/haploid_rename" \
        "${ROOT_DIR}/haploid_rename"; do
        if [[ -d "${candidate}" ]]; then
            PZI_RENAME_DIR="${candidate}"
            break
        fi
    done
fi
[[ -n "${PZI_RENAME_DIR:-}" ]] || die "Cannot locate the PZI haploid_rename directory"

PSP_DIR="${PSP_DIR:-${ROOT_DIR}/data/psp/haploid/genome}"
PPA_RENAME_DIR="${PPA_RENAME_DIR:-${ROOT_DIR}/data/ppa/haploid_rename}"

mkdir -p "${KS_ROOT}/logs" "${STAGE_DIR}/canonical" "${INPUT_DIR}" \
    "${PAIR_DIR}" "${SUMMARY_DIR}" "${FIG_DIR}"
rm -f "${JCVI_DIR}/COMPLETE"

stage_label() {
    local label="$1" genome="$2" gff="$3" pep="$4" cds="$5"
    local destination="${STAGE_DIR}/canonical/${label}"
    for source in "${genome}" "${gff}" "${pep}" "${cds}"; do
        require_file "${source}"
    done
    mkdir -p "${destination}"
    ln -sfn "${genome}" "${destination}/genome.fa"
    ln -sfn "${gff}" "${destination}/primary.gff3"
    ln -sfn "${pep}" "${destination}/primary.pep.fa"
    ln -sfn "${cds}" "${destination}/primary.cds.fa"
}

for hap in A1 A2 B1 B2; do
    label="pzi_hap${hap}"
    stage_label "${label}" \
        "${PZI_RENAME_DIR}/genome/${label}.fa" \
        "${PZI_RENAME_DIR}/gff/${label}.gff3" \
        "${PZI_RENAME_DIR}/pep/${label}.pep.fa" \
        "${PZI_RENAME_DIR}/cds/${label}.cds.fa"
done

for hap in A1 A2 B1 B2; do
    label="psp_hap${hap}"
    prefix="${PSP_DIR}/${label}.edta_masked.helixer"
    stage_label "${label}" \
        "${PSP_DIR}/${label}.fa" \
        "${prefix}.genes.gff3" \
        "${prefix}.protein.fa" \
        "${prefix}.cds.fa"
done

for hap in A B; do
    label="ppa_hap${hap}"
    stage_label "${label}" \
        "${PPA_RENAME_DIR}/genome/${label}.fa" \
        "${PPA_RENAME_DIR}/gff/${label}.gff3" \
        "${PPA_RENAME_DIR}/pep/${label}.pep.fa" \
        "${PPA_RENAME_DIR}/cds/${label}.cds.fa"
done

stage_label "pcam" \
    "${ROOT_DIR}/data/pcam/GCA_050084535.1_PCAv2.0_genomic.fna.gz" \
    "${ROOT_DIR}/data/pcam/GCA_050084535.1_PCAv2.0_genomic.gff.gz" \
    "${ROOT_DIR}/data/pcam/GCA_050084535.1_PCAv2.0_protein.faa.gz" \
    "${ROOT_DIR}/data/pcam/GCA_050084535.1_PCAv2.0_cds_from_genomic.fna.gz"

stage_label "ppe" \
    "${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fna.gz" \
    "${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_genomic.gff.gz" \
    "${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_protein.faa.gz" \
    "${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_cds_from_genomic.fna.gz"

stage_label "gtr" \
    "${ROOT_DIR}/data/gtr/GCA_018257905.1_V1_genomic.fna" \
    "${ROOT_DIR}/data/gtr/Gtr.NCBI.gene.gff3" \
    "${ROOT_DIR}/data/gtr/Gtr.NCBI.prot.faa" \
    "${ROOT_DIR}/data/gtr/Gtr.NCBI.gene.modId.fna"

LABELS=(
    pzi_hapA1 pzi_hapA2 pzi_hapB1 pzi_hapB2
    psp_hapA1 psp_hapA2 psp_hapB1 psp_hapB2
    ppa_hapA ppa_hapB
    pcam ppe gtr
)

cat > "${SUMMARY_DIR}/pair_plan.tsv" <<'EOF'
comparison	query_label	ref_label	comparison_type	analysis_role
pzi_hapA1_vs_pzi_hapB1	pzi_hapA1	pzi_hapB1	within_species_homeolog	primary
psp_hapA1_vs_psp_hapB1	psp_hapA1	psp_hapB1	within_species_homeolog	primary
ppa_hapA_vs_ppa_hapB	ppa_hapA	ppa_hapB	within_species_homeolog	primary
pzi_hapA1_vs_psp_hapA1	pzi_hapA1	psp_hapA1	cross_species_A	primary
pzi_hapB1_vs_psp_hapB1	pzi_hapB1	psp_hapB1	cross_species_B	primary
pzi_hapA1_vs_ppa_hapA	pzi_hapA1	ppa_hapA	cross_species_A	primary
pzi_hapB1_vs_ppa_hapB	pzi_hapB1	ppa_hapB	cross_species_B	primary
pzi_hapA1_vs_pcam	pzi_hapA1	pcam	cross_species_reference	reference
pzi_hapB1_vs_pcam	pzi_hapB1	pcam	cross_species_reference	reference
pzi_hapA1_vs_ppe	pzi_hapA1	ppe	cross_species_reference	reference
pzi_hapB1_vs_ppe	pzi_hapB1	ppe	cross_species_reference	reference
pzi_hapA1_vs_gtr	pzi_hapA1	gtr	cross_species_reference	reference
pzi_hapB1_vs_gtr	pzi_hapB1	gtr	cross_species_reference	reference
EOF

cat > "${SUMMARY_DIR}/sensitivity_pair_plan.tsv" <<'EOF'
comparison	query_label	ref_label	comparison_type	analysis_role
pzi_hapA2_vs_pzi_hapB2	pzi_hapA2	pzi_hapB2	within_species_homeolog	sensitivity
psp_hapA2_vs_psp_hapB2	psp_hapA2	psp_hapB2	within_species_homeolog	sensitivity
pzi_hapA2_vs_psp_hapA2	pzi_hapA2	psp_hapA2	cross_species_A	sensitivity
pzi_hapB2_vs_psp_hapB2	pzi_hapB2	psp_hapB2	cross_species_B	sensitivity
EOF

source "${CONDA_SH}"
conda activate "${JCVI_ENV}"
export MPLBACKEND=Agg
export MPLCONFIGDIR="${KS_ROOT}/mplconfig"
mkdir -p "${MPLCONFIGDIR}"

command -v python >/dev/null 2>&1 || die "python unavailable in ${JCVI_ENV}"
command -v diamond >/dev/null 2>&1 || die "diamond unavailable in ${JCVI_ENV}"
python -c 'import jcvi' >/dev/null 2>&1 || die "jcvi is not importable in ${JCVI_ENV}"
require_file "${PREPARE_SCRIPT}"
require_file "${PAIR_SCRIPT}"
require_file "${NATIVE_INPUT_SCRIPT}"

python "${PREPARE_SCRIPT}" \
    --analysis-dir "${STAGE_DIR}" \
    --outdir "${INPUT_DIR}" \
    --summary "${SUMMARY_DIR}/input_summary.tsv" \
    --labels "${LABELS[@]}"

awk -F '\t' -v min_pep="${MIN_PEP_MATCH_FRAC}" -v min_cds="${MIN_CDS_MATCH_FRAC}" '
BEGIN {
    OFS="\t";
    print "label", "pep_records", "usable_bed_pep_records", "cds_matched_records", "pep_match_fraction", "cds_match_fraction", "status"
}
NR > 1 {
    pep=$3+0; usable=$4+0; cds=$5+0;
    pep_frac=(pep > 0 ? usable/pep : 0);
    cds_frac=(usable > 0 ? cds/usable : 0);
    status=(usable > 0 && pep_frac >= min_pep && cds_frac >= min_cds ? "PASS" : "FAIL");
    printf "%s\t%d\t%d\t%d\t%.6f\t%.6f\t%s\n", $1, pep, usable, cds, pep_frac, cds_frac, status;
}' "${SUMMARY_DIR}/input_summary.tsv" > "${SUMMARY_DIR}/id_consistency_qc.tsv"

if awk -F '\t' 'NR > 1 && $7 != "PASS" {bad=1} END {exit bad}' "${SUMMARY_DIR}/id_consistency_qc.tsv"; then
    :
else
    cat "${SUMMARY_DIR}/id_consistency_qc.tsv" >&2
    die "GFF/PEP/CDS ID consistency failed; do not proceed to JCVI"
fi

run_comparisons() {
    local plan="$1"
    while IFS=$'\t' read -r comparison query ref comparison_type analysis_role; do
        [[ "${comparison}" == "comparison" ]] && continue
        work="${PAIR_DIR}/${comparison}"
        mkdir -p "${work}"
        for label in "${query}" "${ref}"; do
            for ext in bed pep cds; do
                require_file "${INPUT_DIR}/${label}/${label}.${ext}"
                ln -sfn "${INPUT_DIR}/${label}/${label}.${ext}" "${work}/${label}.${ext}"
            done
        done
        cd "${work}"
        raw="${query}.${ref}.anchors"
        screened="${query}.${ref}.anchors.new"
        simple="${query}.${ref}.anchors.simple"
        if [[ ! -s "${raw}" ]]; then
            python -m jcvi.compara.catalog ortholog "${query}" "${ref}" \
                --dbtype=prot --align_soft=diamond_blastp \
                --cscore="${CSCORE}" --min_size="${MIN_ANCHORS}" \
                --cpus="${THREADS}" --no_strip_names --no_dotplot
        fi
        require_file "${raw}"
        if [[ ! -s "${screened}" ]]; then
            python -m jcvi.compara.synteny screen \
                --qbed="${query}.bed" --sbed="${ref}.bed" \
                --minspan="${MINSPAN}" --simple "${raw}" "${screened}"
        fi
        require_file "${screened}"
        if [[ ! -s "${simple}" ]]; then
            [[ -s "${screened}.simple" ]] || die "Missing simple anchors: ${comparison}"
            ln -sfn "${screened}.simple" "${simple}"
        fi

        if [[ "${comparison_type}" != "cross_species_reference" ]]; then
            python -m jcvi.graphics.dotplot "${raw}" \
                --qbed="${query}.bed" --sbed="${ref}.bed" \
                --figsize=12x12 --dpi=300 --format=png --notex \
                --outfile="${FIG_DIR}/${comparison}.dotplot.png"

            python "${NATIVE_INPUT_SCRIPT}" \
                --query "${query}" --ref "${ref}" \
                --query-seqids "${INPUT_DIR}/${query}/${query}.seqids" \
                --ref-seqids "${INPUT_DIR}/${ref}/${ref}.seqids" \
                --query-bed "${query}.bed" --ref-bed "${ref}.bed" \
                --anchors-simple "${simple}" \
                --out-seqids "${FIG_DIR}/${comparison}.seqids" \
                --out-layout "${FIG_DIR}/${comparison}.layout"
            python -m jcvi.graphics.karyotype \
                "${FIG_DIR}/${comparison}.seqids" "${FIG_DIR}/${comparison}.layout" \
                --format=png --notex --outfile="${FIG_DIR}/${comparison}.karyotype.png"
        fi
        cd "${ROOT_DIR}"
    done < "${plan}"
}

run_comparisons "${SUMMARY_DIR}/pair_plan.tsv"
run_comparisons "${SUMMARY_DIR}/sensitivity_pair_plan.tsv"

{
    head -n 1 "${SUMMARY_DIR}/pair_plan.tsv"
    tail -n +2 "${SUMMARY_DIR}/pair_plan.tsv"
    tail -n +2 "${SUMMARY_DIR}/sensitivity_pair_plan.tsv"
} > "${SUMMARY_DIR}/all_pair_plan.tsv"

python "${PAIR_SCRIPT}" \
    --pair-plan "${SUMMARY_DIR}/all_pair_plan.tsv" \
    --pair-root "${PAIR_DIR}" \
    --out-pairs "${SUMMARY_DIR}/jcvi_syntenic_pairs.tsv" \
    --out-blocks "${SUMMARY_DIR}/jcvi_block_summary.tsv" \
    --out-status "${SUMMARY_DIR}/jcvi_pair_status.tsv"

touch "${JCVI_DIR}/COMPLETE"
echo "[OK] Subgenome JCVI comparisons completed: ${JCVI_DIR}"
