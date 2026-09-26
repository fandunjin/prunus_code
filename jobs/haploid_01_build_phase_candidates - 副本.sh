#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
WORK_DIR="${OUT_DIR}/work"
EVIDENCE_DIR="${OUT_DIR}/evidence"
FINAL_HAP_DIR="${OUT_DIR}/final_haplotypes"
PCAM_DIR="${ROOT_DIR}/data/pcam"
PPA_DIR="${ROOT_DIR}/data/ppa"
DETAILS_DIR="${ROOT_DIR}/data/details"
PPE_DIR="${ROOT_DIR}/data/ppe"
PAR_DIR="${ROOT_DIR}/data/par"
PDU_DIR="${ROOT_DIR}/data/pdu"
PAV_DIR="${ROOT_DIR}/data/pav"
PMU_DIR="${ROOT_DIR}/data/pmu"
FVE_DIR="${ROOT_DIR}/data/fve"
RCH_DIR="${ROOT_DIR}/data/rch"
PCOM_DIR="${ROOT_DIR}/data/pcom"
FINAL_DIR="${ROOT_DIR}/final_assembly"
FINAL_ASM_DIR="${ROOT_DIR}/final_assembly/haplotype_provenance"

PZI_FASTA="${FINAL_DIR}/pzi_2n32.chromosomes.fa"
PCAM_FASTA="${PCAM_DIR}/GCA_050084535.1_PCAv2.0_genomic.fna.gz"
PCAM_REPORT="${DETAILS_DIR}/pcam/GCA_050084535.1_PCAv2.0_assembly_report.txt"
PPA_FASTA="${PPA_DIR}/GCA_964059315.1_drPruPadu1.1_genomic.fna.gz"
PPA_REPORT="${DETAILS_DIR}/ppa/GCA_964059315.1_drPruPadu1.1_assembly_report.txt"
PPE_FASTA="${PPE_DIR}/GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fna.gz"
PPE_REPORT="${DETAILS_DIR}/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_assembly_report.txt"
PAR_FASTA="${PAR_DIR}/GCA_056632545.1_CAAS_ParmT2T_1.0_genomic.fna.gz"
PAR_REPORT="${DETAILS_DIR}/par/GCA_056632545.1_CAAS_ParmT2T_1.0_assembly_report.txt"
PDU_FASTA="${PDU_DIR}/GCF_902201215.1_ALMONDv2_genomic.fna.gz"
PDU_REPORT="${DETAILS_DIR}/pdu/GCF_902201215.1_ALMONDv2_assembly_report.txt"
PAV_FASTA="${PAV_DIR}/GCA_014155035.1_NE-ME-NP-PH_genomic.fna.gz"
PAV_REPORT="${DETAILS_DIR}/pav/GCA_014155035.1_NE-ME-NP-PH_assembly_report.txt"
PMU_FASTA="${PMU_DIR}/GCF_000346735.1_P.mume_V1.0_genomic.fna.gz"
PMU_REPORT="${DETAILS_DIR}/pmu/GCF_000346735.1_P.mume_V1.0_assembly_report.txt"
FVE_FASTA="${FVE_DIR}/GCF_000184155.1_FraVesHawaii_1.0_genomic.fna.gz"
FVE_REPORT="${DETAILS_DIR}/fve/GCF_000184155.1_FraVesHawaii_1.0_assembly_report.txt"
RCH_FASTA="${RCH_DIR}/GCF_002994745.2_RchiOBHm-V2_genomic.fna.gz"
RCH_REPORT="${DETAILS_DIR}/rch/GCF_002994745.2_RchiOBHm-V2_assembly_report.txt"
PCOM_FASTA="${PCOM_DIR}/GCF_963583255.1_drPyrComm1.1_genomic.fna.gz"
PCOM_REPORT="${DETAILS_DIR}/pcom/GCF_963583255.1_drPyrComm1.1_assembly_report.txt"
CHR_HAP_PROV="${FINAL_ASM_DIR}/pzi_2n32.chr_hap_provenance.tsv"
CORRECTED_HOMOLOG_PAIRS="${FINAL_ASM_DIR}/pzi_2n32.corrected_homolog_pairs.tsv"
PZI_VS_PCAM_PAF="${WORK_DIR}/pzi_vs_pcam.asm10.paf"
PZI_VS_PPA_PAF="${WORK_DIR}/pzi_vs_ppa.asm10.paf"
PZI_VS_PPE_PAF="${WORK_DIR}/pzi_vs_ppe.asm10.paf"
PZI_VS_PAR_PAF="${WORK_DIR}/pzi_vs_par.asm10.paf"
PZI_VS_PDU_PAF="${WORK_DIR}/pzi_vs_pdu.asm10.paf"
PZI_VS_PAV_PAF="${WORK_DIR}/pzi_vs_pav.asm10.paf"
PZI_VS_PMU_PAF="${WORK_DIR}/pzi_vs_pmu.asm10.paf"
PZI_VS_FVE_PAF="${WORK_DIR}/pzi_vs_fve.asm10.paf"
PZI_VS_RCH_PAF="${WORK_DIR}/pzi_vs_rch.asm10.paf"
PZI_VS_PCOM_PAF="${WORK_DIR}/pzi_vs_pcom.asm10.paf"
CHR32_SELF_PAF="${WORK_DIR}/pzi_chr32_self.asm10.paf"

CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
MINIMAP_ENV="minimap"
HELPER="${ROOT_DIR}/scripts/haploid_phase_helper.py"
FINAL_HELPER="${ROOT_DIR}/scripts/final_haplotype_assignment.py"
VALIDATE_FINAL_HELPER="${ROOT_DIR}/scripts/validate_final_haplotypes.py"
DOTPLOT_SCRIPT="${ROOT_DIR}/scripts/plot_minimap_synteny_dotplot.py"
HEATMAP_SCRIPT="${ROOT_DIR}/scripts/plot_minimap_synteny_heatmap.py"
THREADS=24
MINIMAP_PRESET="asm10"
SEGMENT_BIN_SIZE=500000

require_file() {
    local file_path="$1"
    local label="$2"
    if [ ! -s "${file_path}" ]; then
        echo "ERROR: ${label} does not exist or is empty: ${file_path}" >&2
        exit 1
    fi
}

mkdir -p "${OUT_DIR}" "${WORK_DIR}" "${EVIDENCE_DIR}" "${FINAL_HAP_DIR}" "${ROOT_DIR}/logs"

require_file "${PZI_FASTA}" "final pzi chromosomes FASTA"
require_file "${PCAM_FASTA}" "P. campanulata FASTA"
require_file "${PCAM_REPORT}" "P. campanulata assembly report"
require_file "${PPA_FASTA}" "P. padus FASTA"
require_file "${PPA_REPORT}" "P. padus assembly report"
require_file "${PPE_FASTA}" "P. persica FASTA"
require_file "${PPE_REPORT}" "P. persica assembly report"
require_file "${PAR_FASTA}" "P. armeniaca FASTA"
require_file "${PAR_REPORT}" "P. armeniaca assembly report"
require_file "${PDU_FASTA}" "P. dulcis FASTA"
require_file "${PDU_REPORT}" "P. dulcis assembly report"
require_file "${PAV_FASTA}" "P. avium FASTA"
require_file "${PAV_REPORT}" "P. avium assembly report"
require_file "${PMU_FASTA}" "P. mume FASTA"
require_file "${PMU_REPORT}" "P. mume assembly report"
require_file "${FVE_FASTA}" "Fragaria vesca FASTA"
require_file "${FVE_REPORT}" "Fragaria vesca assembly report"
require_file "${RCH_FASTA}" "Rosa chinensis FASTA"
require_file "${RCH_REPORT}" "Rosa chinensis assembly report"
require_file "${PCOM_FASTA}" "Pyrus communis FASTA"
require_file "${PCOM_REPORT}" "Pyrus communis assembly report"
require_file "${CHR_HAP_PROV}" "chromosome hap provenance table"
require_file "${CORRECTED_HOMOLOG_PAIRS}" "corrected homolog pairs"
require_file "${CONDA_SH}" "Conda init script"
require_file "${HELPER}" "haploid phase helper"
require_file "${FINAL_HELPER}" "final haplotype assignment helper"
require_file "${VALIDATE_FINAL_HELPER}" "final haplotype validation helper"
require_file "${DOTPLOT_SCRIPT}" "minimap dotplot script"
require_file "${HEATMAP_SCRIPT}" "minimap heatmap script"

source "${CONDA_SH}"
conda activate "${MINIMAP_ENV}"

if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi

echo "=== Haploid phase-candidate workflow ==="
echo "PZI FASTA   : ${PZI_FASTA}"
echo "PCAM FASTA  : ${PCAM_FASTA}"
echo "PPA FASTA   : ${PPA_FASTA}"
echo "N8 refs     : ppe par pdu pav pmu"
echo "Outgroups   : fve rch pcom"
echo "Fixed rule  : x=8 only; no chromosome-specific prior"
echo "Work dir    : ${WORK_DIR}"
echo "Evidence dir: ${EVIDENCE_DIR}"
echo "Final dir   : ${FINAL_HAP_DIR}"
echo "Threads     : ${THREADS}"
echo "Preset      : ${MINIMAP_PRESET}"

if [ ! -s "${PZI_VS_PCAM_PAF}" ]; then
    echo ">>> Align pzi chr32 vs P. campanulata"
    minimap2 -x "${MINIMAP_PRESET}" -t "${THREADS}" "${PCAM_FASTA}" "${PZI_FASTA}" > "${PZI_VS_PCAM_PAF}"
fi

if [ ! -s "${PZI_VS_PPA_PAF}" ]; then
    echo ">>> Align pzi chr32 vs P. padus"
    minimap2 -x "${MINIMAP_PRESET}" -t "${THREADS}" "${PPA_FASTA}" "${PZI_FASTA}" > "${PZI_VS_PPA_PAF}"
fi

if [ ! -s "${PZI_VS_PPE_PAF}" ]; then
    echo ">>> Align pzi chr32 vs P. persica"
    minimap2 -x "${MINIMAP_PRESET}" -t "${THREADS}" "${PPE_FASTA}" "${PZI_FASTA}" > "${PZI_VS_PPE_PAF}"
fi

if [ ! -s "${PZI_VS_PAR_PAF}" ]; then
    echo ">>> Align pzi chr32 vs P. armeniaca"
    minimap2 -x "${MINIMAP_PRESET}" -t "${THREADS}" "${PAR_FASTA}" "${PZI_FASTA}" > "${PZI_VS_PAR_PAF}"
fi

if [ ! -s "${PZI_VS_PDU_PAF}" ]; then
    echo ">>> Align pzi chr32 vs P. dulcis"
    minimap2 -x "${MINIMAP_PRESET}" -t "${THREADS}" "${PDU_FASTA}" "${PZI_FASTA}" > "${PZI_VS_PDU_PAF}"
fi

if [ ! -s "${PZI_VS_PAV_PAF}" ]; then
    echo ">>> Align pzi chr32 vs P. avium"
    minimap2 -x "${MINIMAP_PRESET}" -t "${THREADS}" "${PAV_FASTA}" "${PZI_FASTA}" > "${PZI_VS_PAV_PAF}"
fi

if [ ! -s "${PZI_VS_PMU_PAF}" ]; then
    echo ">>> Align pzi chr32 vs P. mume"
    minimap2 -x "${MINIMAP_PRESET}" -t "${THREADS}" "${PMU_FASTA}" "${PZI_FASTA}" > "${PZI_VS_PMU_PAF}"
fi

if [ ! -s "${PZI_VS_FVE_PAF}" ]; then
    echo ">>> Align pzi chr32 vs Fragaria vesca"
    minimap2 -x "${MINIMAP_PRESET}" -t "${THREADS}" "${FVE_FASTA}" "${PZI_FASTA}" > "${PZI_VS_FVE_PAF}"
fi

if [ ! -s "${PZI_VS_RCH_PAF}" ]; then
    echo ">>> Align pzi chr32 vs Rosa chinensis"
    minimap2 -x "${MINIMAP_PRESET}" -t "${THREADS}" "${RCH_FASTA}" "${PZI_FASTA}" > "${PZI_VS_RCH_PAF}"
fi

if [ ! -s "${PZI_VS_PCOM_PAF}" ]; then
    echo ">>> Align pzi chr32 vs Pyrus communis"
    minimap2 -x "${MINIMAP_PRESET}" -t "${THREADS}" "${PCOM_FASTA}" "${PZI_FASTA}" > "${PZI_VS_PCOM_PAF}"
fi

if [ ! -s "${CHR32_SELF_PAF}" ]; then
    echo ">>> Align pzi chr32 all-vs-all"
    minimap2 -x "${MINIMAP_PRESET}" -t "${THREADS}" "${PZI_FASTA}" "${PZI_FASTA}" > "${CHR32_SELF_PAF}"
fi

echo ">>> Build quartet and naming candidate tables"
"${PYTHON_BIN}" "${HELPER}" \
    --root-dir "${ROOT_DIR}" \
    --out-dir "${EVIDENCE_DIR}" \
    --pzi-chrom-fasta "${PZI_FASTA}" \
    --pcam-fasta "${PCAM_FASTA}" \
    --pcam-assembly-report "${PCAM_REPORT}" \
    --ppa-fasta "${PPA_FASTA}" \
    --ppa-assembly-report "${PPA_REPORT}" \
    --chr-hap-provenance "${CHR_HAP_PROV}" \
    --corrected-homolog-pairs "${CORRECTED_HOMOLOG_PAIRS}" \
    --pzi-vs-pcam-paf "${PZI_VS_PCAM_PAF}" \
    --pzi-vs-ppa-paf "${PZI_VS_PPA_PAF}" \
    --n8-reference "ppe|${PPE_FASTA}|${PPE_REPORT}|${PZI_VS_PPE_PAF}|8" \
    --n8-reference "par|${PAR_FASTA}|${PAR_REPORT}|${PZI_VS_PAR_PAF}|8" \
    --n8-reference "pdu|${PDU_FASTA}|${PDU_REPORT}|${PZI_VS_PDU_PAF}|8" \
    --n8-reference "pav|${PAV_FASTA}|${PAV_REPORT}|${PZI_VS_PAV_PAF}|8" \
    --n8-reference "pmu|${PMU_FASTA}|${PMU_REPORT}|${PZI_VS_PMU_PAF}|8" \
    --outgroup-reference "fve|${FVE_FASTA}|${FVE_REPORT}|${PZI_VS_FVE_PAF}|7" \
    --outgroup-reference "rch|${RCH_FASTA}|${RCH_REPORT}|${PZI_VS_RCH_PAF}|7" \
    --outgroup-reference "pcom|${PCOM_FASTA}|${PCOM_REPORT}|${PZI_VS_PCOM_PAF}|17" \
    --focus-chrs "" \
    --segment-bin-size "${SEGMENT_BIN_SIZE}" \
    --chr32-self-paf "${CHR32_SELF_PAF}" \
    --min-mapq 10 \
    --min-aln-len 10000

draw_reference_plots() {
    local label="$1"
    local paf="$2"
    local species_name="$3"

    echo ">>> Draw pzi-vs-${label} dotplot and heatmap"
    "${PYTHON_BIN}" "${DOTPLOT_SCRIPT}" \
        --paf "${paf}" \
        --out-svg "${EVIDENCE_DIR}/pzi_vs_${label}.asm10.dotplot.svg" \
        --title "P. zippeliana chr32 vs ${species_name}" \
        --axis-colors "${EVIDENCE_DIR}/pzi_vs_${label}_axis_colors.tsv" \
        --axis-segments "${EVIDENCE_DIR}/pzi_vs_${label}_axis_segments.tsv"
    "${PYTHON_BIN}" "${HEATMAP_SCRIPT}" \
        --summary "${EVIDENCE_DIR}/pzi_vs_${label}_alignment_summary.tsv" \
        --out-svg "${EVIDENCE_DIR}/pzi_vs_${label}.asm10.heatmap.svg" \
        --title "P. zippeliana chr32 vs ${species_name} aligned blocks" \
        --x-label "P. zippeliana chromosomes" \
        --y-label "${species_name} chromosomes" \
        --axis-colors "${EVIDENCE_DIR}/pzi_vs_${label}_axis_colors.tsv"
}

draw_outgroup_plots() {
    local label="$1"
    local paf="$2"
    local species_name="$3"

    echo ">>> Draw pzi-vs-${label} outgroup dotplot and heatmap"
    "${PYTHON_BIN}" "${DOTPLOT_SCRIPT}" \
        --paf "${paf}" \
        --out-svg "${EVIDENCE_DIR}/pzi_vs_${label}.asm10.dotplot.svg" \
        --title "P. zippeliana chr32 vs ${species_name}" \
        --min-alignment-len 10000 \
        --min-query-aggregate 100000 \
        --min-target-aggregate 100000 \
        --min-visible-px 2.4 \
        --axis-colors "${EVIDENCE_DIR}/pzi_vs_${label}_axis_colors.tsv" \
        --axis-segments "${EVIDENCE_DIR}/pzi_vs_${label}_axis_segments.tsv"
    "${PYTHON_BIN}" "${HEATMAP_SCRIPT}" \
        --summary "${EVIDENCE_DIR}/pzi_vs_${label}_alignment_summary.tsv" \
        --out-svg "${EVIDENCE_DIR}/pzi_vs_${label}.asm10.heatmap.svg" \
        --title "P. zippeliana chr32 vs ${species_name} aligned blocks" \
        --x-label "P. zippeliana chromosomes" \
        --y-label "${species_name} chromosomes" \
        --axis-colors "${EVIDENCE_DIR}/pzi_vs_${label}_axis_colors.tsv"
}

draw_reference_plots "pcam" "${PZI_VS_PCAM_PAF}" "P. campanulata"
draw_reference_plots "ppa" "${PZI_VS_PPA_PAF}" "P. padus"
draw_reference_plots "ppe" "${PZI_VS_PPE_PAF}" "P. persica"
draw_reference_plots "par" "${PZI_VS_PAR_PAF}" "P. armeniaca"
draw_reference_plots "pdu" "${PZI_VS_PDU_PAF}" "P. dulcis"
draw_reference_plots "pav" "${PZI_VS_PAV_PAF}" "P. avium"
draw_reference_plots "pmu" "${PZI_VS_PMU_PAF}" "P. mume"
draw_outgroup_plots "fve" "${PZI_VS_FVE_PAF}" "Fragaria vesca"
draw_outgroup_plots "rch" "${PZI_VS_RCH_PAF}" "Rosa chinensis"
draw_outgroup_plots "pcom" "${PZI_VS_PCOM_PAF}" "Pyrus communis"

echo ">>> Draw chr32 all-vs-all dotplot and heatmap"
"${PYTHON_BIN}" "${DOTPLOT_SCRIPT}" \
    --paf "${CHR32_SELF_PAF}" \
    --out-svg "${EVIDENCE_DIR}/pzi_chr32_self.asm10.dotplot.svg" \
    --title "P. zippeliana chr32 all-vs-all"
"${PYTHON_BIN}" "${HEATMAP_SCRIPT}" \
    --summary "${EVIDENCE_DIR}/chr32_self_alignment_summary.tsv" \
    --out-svg "${EVIDENCE_DIR}/pzi_chr32_self.asm10.heatmap.svg" \
    --title "P. zippeliana chr32 all-vs-all aligned blocks" \
    --x-label "P. zippeliana query chromosomes" \
    --y-label "P. zippeliana reference chromosomes"

echo ">>> Finalize A1/A2/B1/B2 haplotype assignment and export renamed FASTA"
"${PYTHON_BIN}" "${FINAL_HELPER}" \
    --root-dir "${ROOT_DIR}" \
    --evidence-dir "${EVIDENCE_DIR}" \
    --out-dir "${FINAL_HAP_DIR}" \
    --pzi-chrom-fasta "${PZI_FASTA}" \
    --chr-hap-provenance "${CHR_HAP_PROV}" \
    --corrected-homolog-pairs "${CORRECTED_HOMOLOG_PAIRS}"

echo ">>> Validate final x=8 haplotype outputs"
"${PYTHON_BIN}" "${VALIDATE_FINAL_HELPER}" \
    --final-hap-dir "${FINAL_HAP_DIR}"

echo "=== Haploid workflow outputs ==="
echo "Evidence tables and plots:"
echo "${EVIDENCE_DIR}/multi_reference_consensus_quartets.tsv"
echo "${EVIDENCE_DIR}/outgroup_chr_votes.tsv"
echo "${EVIDENCE_DIR}/pzi_vs_pcam.asm10.dotplot.svg"
echo "${EVIDENCE_DIR}/pzi_vs_ppa.asm10.dotplot.svg"
echo "${EVIDENCE_DIR}/pzi_vs_ppe.asm10.dotplot.svg"
echo "${EVIDENCE_DIR}/pzi_vs_par.asm10.dotplot.svg"
echo "${EVIDENCE_DIR}/pzi_vs_pdu.asm10.dotplot.svg"
echo "${EVIDENCE_DIR}/pzi_vs_pav.asm10.dotplot.svg"
echo "${EVIDENCE_DIR}/pzi_vs_pmu.asm10.dotplot.svg"
echo "${EVIDENCE_DIR}/pzi_vs_fve.asm10.dotplot.svg"
echo "${EVIDENCE_DIR}/pzi_vs_rch.asm10.dotplot.svg"
echo "${EVIDENCE_DIR}/pzi_vs_pcom.asm10.dotplot.svg"
echo "Final haplotype assignment:"
echo "${FINAL_HAP_DIR}/FINAL_HAPLOTYPE_OVERVIEW.tsv"
echo "${FINAL_HAP_DIR}/FINAL_REVIEW_NOTES.tsv"
echo "${FINAL_HAP_DIR}/FINAL_FILES.tsv"
echo "${FINAL_HAP_DIR}/FINAL_VALIDATION.tsv"
echo "${FINAL_HAP_DIR}/final_assignment_summary.tsv"
echo "${FINAL_HAP_DIR}/pzi_A1A2B1B2.rename_map.tsv"
echo "${FINAL_HAP_DIR}/pzi_hapA1.fa"
echo "${FINAL_HAP_DIR}/pzi_hapA2.fa"
echo "${FINAL_HAP_DIR}/pzi_hapB1.fa"
echo "${FINAL_HAP_DIR}/pzi_hapB2.fa"
echo "${FINAL_HAP_DIR}/pzi_A1A2B1B2.renamed.chr32.fa"
echo "${FINAL_HAP_DIR}/details/candidate_quartet_scores.tsv"
echo "${FINAL_HAP_DIR}/details/global_matching_solutions.tsv"
echo "${FINAL_HAP_DIR}/annotation_next_steps.tsv"
echo "${FINAL_HAP_DIR}/README_FINAL_HAPLOTYPES.md"
