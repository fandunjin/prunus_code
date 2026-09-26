#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# PSP transcript models through -e, and Helixer models through -c.
# Run Helixer first with: jsub < jobs/pzi_32_helixer.sh
# Then submit: jsub < jobs/pzi_annotation_eviann.sh

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/annotation}"
TRANS_DIR="${TRANS_DIR:-${ROOT_DIR}/transcriptome}"
GENOME_FA="${GENOME_FA:-${ROOT_DIR}/final_assembly/pzi_2n32.chrom/pzi_2n32.chromosomes.fa}"
if [[ ! -s "${GENOME_FA}" && -s "${ROOT_DIR}/final_assembly/pzi_2n32.chromosomes.fa" ]]; then GENOME_FA="${ROOT_DIR}/final_assembly/pzi_2n32.chromosomes.fa"; fi
HELIXER_GFF="${HELIXER_GFF:-${ANNOTATION_DIR}/genes/pzi32_helixer_edta_masked/pzi_2n32.chromosomes.edta_masked.helixer.gff3}"
TE_GFF3="${TE_GFF3:-${ANNOTATION_DIR}/chr32_edta/pzi_2n32.chromosomes.fa.mod.EDTA.TEanno.gff3}"
EVIDENCE_DIR="${ANNOTATION_DIR}/genes/pzi_eviann/evidence"
EVIDENCE_TRANSCRIPTS="${EVIDENCE_TRANSCRIPTS:-${EVIDENCE_DIR}/pzi_illumina_pacbio.transcripts.fa}"
EVIDENCE_GTF="${EVIDENCE_GTF:-${EVIDENCE_DIR}/pzi_illumina_pacbio.transcripts.gtf}"
READS_LIST="${READS_LIST:-${EVIDENCE_DIR}/pzi_eviann_reads.txt}"
EVIANN_TRANSCRIPTS="${EVIANN_TRANSCRIPTS:-${EVIDENCE_DIR}/pzi_target_plus_psp.transcripts.fa}"
OUT_DIR="${ANNOTATION_DIR}/genes/pzi_eviann"
POST_DIR="${OUT_DIR}/postprocess"
PREFIX="pzi32.eviann"
THREADS="${THREADS:-32}"
FORCE="${FORCE:-0}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
EVIANN_ENV="${EVIANN_ENV:-eviann}"
PY_ENV="${PY_ENV:-haphic}"
SWISSPROT="${SWISSPROT:-${HOME}/CX/database/SwissProt/uniprot_sprot.fasta}"
POSTPROCESS_SCRIPT="${POSTPROCESS_SCRIPT:-${ROOT_DIR}/scripts/helixer_postprocess_qc.py}"
GFF_STATS_SCRIPT="${GFF_STATS_SCRIPT:-${ROOT_DIR}/scripts/gff3_basic_stats.py}"

RELATED_TRANSCRIPTS="${RELATED_TRANSCRIPTS:-${ROOT_DIR}/data/psp/annotation/genes/psp32_eviann/evidence/psp_illumina.transcripts.fa}"
PROTEIN_EVIDENCE_DEFAULT="${EVIDENCE_DIR}/ppe_pcam.proteins.fa"
RELATED_PROTEINS="${RELATED_PROTEINS:-${PROTEIN_EVIDENCE_DEFAULT}}"
USE_PROTEINS="${USE_PROTEINS:-1}"
USE_SWISSPROT="${USE_SWISSPROT:-0}"
PLOIDY="${PLOIDY:-4}"

RELATED_PROTEIN_SOURCES=(
    "${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_protein.faa.gz"
    "${ROOT_DIR}/data/pcam/GCA_050084535.1_PCAv2.0_protein.faa.gz"
)
RELATED_PROTEIN_LABELS=(ppe pcam)

die() { echo "ERROR: $*" >&2; exit 1; }
require_file() { [[ -s "$1" ]] || die "$2 does not exist or is empty: $1"; }
require_cmd() { command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"; }

case "${USE_PROTEINS}" in
    0|1) ;;
    *) die "USE_PROTEINS must be 0 or 1 (observed: ${USE_PROTEINS})" ;;
esac
case "${USE_SWISSPROT}" in
    0|1) ;;
    *) die "USE_SWISSPROT must be 0 or 1 (observed: ${USE_SWISSPROT})" ;;
esac
[[ "${PLOIDY}" =~ ^[1-9][0-9]*$ ]] || die "PLOIDY must be a positive integer (observed: ${PLOIDY})"

append_fasta() {
    local source="$1"
    local target="$2"
    if [[ -s "${source}" ]]; then
        if [[ "${source}" == *.gz ]]; then
            gzip -dc "${source}" >> "${target}"
        else
            cat "${source}" >> "${target}"
        fi
        return 0
    fi
    return 1
}

append_prefixed_fasta() {
    local source="$1"
    local prefix="$2"
    local target="$3"
    if [[ "${source}" == *.gz ]]; then
        gzip -dc "${source}"
    else
        cat "${source}"
    fi | awk -v prefix="${prefix}" '/^>/{sub(/^>/, ">" prefix "|")} {print}' >> "${target}"
}

build_related_evidence() {
    if [[ "${USE_PROTEINS}" == "1" && "${RELATED_PROTEINS}" == "${PROTEIN_EVIDENCE_DEFAULT}" ]]; then
        : > "${RELATED_PROTEINS}"
        for i in "${!RELATED_PROTEIN_SOURCES[@]}"; do
            source="${RELATED_PROTEIN_SOURCES[$i]}"
            append_prefixed_fasta "${source}" "${RELATED_PROTEIN_LABELS[$i]}" "${RELATED_PROTEINS}" \
                || die "missing PPE/PCAM protein source: ${source}"
        done
    fi
    require_file "${RELATED_TRANSCRIPTS}" "PSP transcript evidence for PZI"
    require_file "${EVIDENCE_TRANSCRIPTS}" "PZI genome-guided transcript evidence"
    local combined_tmp="${EVIANN_TRANSCRIPTS}.tmp.$$"
    : > "${combined_tmp}"
    # Keep the assembled PZI Illumina/PacBio models and add reciprocal PSP
    # transcripts as the related-species evidence.
    append_prefixed_fasta "${EVIDENCE_TRANSCRIPTS}" "pzi" "${combined_tmp}"
    append_prefixed_fasta "${RELATED_TRANSCRIPTS}" "psp" "${combined_tmp}"
    mv -f "${combined_tmp}" "${EVIANN_TRANSCRIPTS}"
    require_file "${EVIANN_TRANSCRIPTS}" "combined PZI target plus PSP transcript evidence"
    if [[ "${USE_PROTEINS}" == "1" ]]; then
        require_file "${RELATED_PROTEINS}" "PPE+PCAM protein evidence"
    fi
}

build_reads_list() {
    local tmp="${READS_LIST}.tmp.$$"
    : > "${tmp}"
    for lib in 14-2_LIBY260601004066 14-2_LIBY260601004067 14-2_LIBY260601004068; do
        local r1="${TRANS_DIR}/01_clean_data/${lib}_clean_R1.fq.gz"
        local r2="${TRANS_DIR}/01_clean_data/${lib}_clean_R2.fq.gz"
        require_file "${r1}" "PZI cleaned Illumina R1 ${lib}"
        require_file "${r2}" "PZI cleaned Illumina R2 ${lib}"
        printf '%s %s\n' "${r1}" "${r2}" >> "${tmp}"
    done
    mv -f "${tmp}" "${READS_LIST}"
}

require_file "${CONDA_SH}" "Conda initialization script"
require_file "${GENOME_FA}" "PZI 2n=32 genome FASTA"
require_file "${HELIXER_GFF}" "Helixer GFF3"
require_file "${TE_GFF3}" "PZI EDTA TE annotation GFF3"
require_file "${POSTPROCESS_SCRIPT}" "Helixer postprocess/QC script"
require_file "${GFF_STATS_SCRIPT}" "GFF3 statistics script"
mkdir -p "${EVIDENCE_DIR}" "${OUT_DIR}" "${POST_DIR}" "${ROOT_DIR}/logs"

source "${CONDA_SH}"
conda activate "${EVIANN_ENV}"
require_cmd eviann.sh
build_related_evidence
build_reads_list

GENOME_BASENAME="$(basename "${GENOME_FA}")"
RAW_GFF="${OUT_DIR}/${GENOME_BASENAME}.pseudo_label.gff"
FUNCTIONAL_GFF="${OUT_DIR}/${GENOME_BASENAME}.functional_note.pseudo_label.gff"

echo "=== EviAnn PZI 2n=32 annotation ==="
echo "Genome       : ${GENOME_FA}"
echo "Reads list   : ${READS_LIST}"
echo "EviAnn tx     : ${EVIANN_TRANSCRIPTS}"
echo "Related tx   : ${RELATED_TRANSCRIPTS}"
echo "Protein      : ${RELATED_PROTEINS} (enabled=${USE_PROTEINS})"
echo "Helixer GFF3 : ${HELIXER_GFF}"
echo "Output       : ${OUT_DIR}"

if [[ ! -s "${RAW_GFF}" || "${FORCE}" == "1" ]]; then
    cd "${OUT_DIR}"
    EVIANN_CMD=(
        eviann.sh
        -g "${GENOME_FA}"
        -t "${THREADS}"
        -d "${PLOIDY}"
        -r "${READS_LIST}"
        -e "${EVIANN_TRANSCRIPTS}"
        -c "${HELIXER_GFF}"
        --verbose
    )
    if [[ "${USE_PROTEINS}" == "1" ]]; then
        EVIANN_CMD+=(-p "${RELATED_PROTEINS}")
    fi
    if [[ "${USE_PROTEINS}" == "1" && "${USE_SWISSPROT}" == "1" && -s "${SWISSPROT}" ]]; then
        EVIANN_CMD+=(-s "${SWISSPROT}")
    fi
    "${EVIANN_CMD[@]}"
else
    echo "Existing EviAnn GFF found, skipping: ${RAW_GFF}"
fi
require_file "${RAW_GFF}" "EviAnn raw GFF3"

FINAL_GFF="${RAW_GFF}"
if [[ -s "${FUNCTIONAL_GFF}" ]]; then
    FINAL_GFF="${FUNCTIONAL_GFF}"
fi

source "${CONDA_SH}"
conda activate "${PY_ENV}"
if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi

POSTPROCESS_CMD=(
    "${PYTHON_BIN}" "${POSTPROCESS_SCRIPT}"
    --genome "${GENOME_FA}"
    --gff3 "${FINAL_GFF}"
    --out-dir "${POST_DIR}"
    --prefix "${PREFIX}"
    --short-cds 300
    --gene-te-fraction 0.80
    --cds-te-fraction 0.50
    --te-gff3 "${TE_GFF3}"
)
"${POSTPROCESS_CMD[@]}"

FINAL_GFF3="${POST_DIR}/${PREFIX}.genes.gff3"
FINAL_CDS="${POST_DIR}/${PREFIX}.cds.fa"
FINAL_PROTEIN="${POST_DIR}/${PREFIX}.protein.fa"
require_file "${FINAL_GFF3}" "postprocessed EviAnn GFF3"
require_file "${FINAL_CDS}" "postprocessed EviAnn CDS"
require_file "${FINAL_PROTEIN}" "postprocessed EviAnn protein FASTA"
"${PYTHON_BIN}" "${GFF_STATS_SCRIPT}" --gff3 "${FINAL_GFF3}" --prefix "${PREFIX}" --out-dir "${POST_DIR}/stats"

{
    echo "PZI EviAnn annotation manifest"
    echo "Genome: ${GENOME_FA}"
    echo "EviAnn reads list: ${READS_LIST}"
    echo "EviAnn transcript FASTA: ${EVIANN_TRANSCRIPTS}"
    echo "Related transcript evidence: ${RELATED_TRANSCRIPTS}"
    echo "Genome-guided transcript FASTA (auxiliary): ${EVIDENCE_TRANSCRIPTS}"
    echo "Genome-guided transcript GTF (auxiliary): ${EVIDENCE_GTF}"
    echo "Related proteins: ${RELATED_PROTEINS}"
    echo "Protein evidence enabled: ${USE_PROTEINS}"
    echo "Swiss-Prot evidence enabled: ${USE_SWISSPROT}"
    echo "Protein source PPE: ${RELATED_PROTEIN_SOURCES[0]}"
    echo "Protein source PCAM: ${RELATED_PROTEIN_SOURCES[1]}"
    echo "Ploidy: ${PLOIDY}"
    echo "Helixer input GFF3: ${HELIXER_GFF}"
    echo "EviAnn raw GFF3: ${RAW_GFF}"
    echo "Final GFF3: ${FINAL_GFF3}"
    echo "Final CDS: ${FINAL_CDS}"
    echo "Final protein: ${FINAL_PROTEIN}"
    echo
    echo "EviAnn is the evidence-supported annotation; Helixer is supplied as complementary prediction evidence."
} > "${OUT_DIR}/pzi_eviann_manifest.txt"

echo "=== EviAnn annotation completed ==="
echo "Final GFF3   : ${FINAL_GFF3}"
echo "Final CDS    : ${FINAL_CDS}"
echo "Final protein: ${FINAL_PROTEIN}"
