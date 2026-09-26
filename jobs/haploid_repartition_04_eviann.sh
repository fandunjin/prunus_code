#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/haploid_repartition/inputs}"
GENE_ROOT="${GENE_ROOT:-${ROOT_DIR}/annotation/repartition/genes}"
EVIDENCE_ROOT="${EVIDENCE_ROOT:-${GENE_ROOT}/eviann_evidence}"
THREADS="${THREADS:-32}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
EVIANN_ENV="${EVIANN_ENV:-eviann}"
PY_ENV="${PY_ENV:-haphic}"
POSTPROCESS_SCRIPT="${ROOT_DIR}/scripts/helixer_postprocess_qc.py"
SWISSPROT="${SWISSPROT:-${HOME}/CX/database/SwissProt/uniprot_sprot.fasta}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing or empty $2: $1"; }
case "${SPECIES}:${LABEL}" in
    pzi:hapA1|pzi:hapA2|pzi:hapB1|pzi:hapB2|psp:hapA1|psp:hapA2|psp:hapB1|psp:hapB2|ppa:hapA|ppa:hapB) ;;
    *) die "Unsupported SPECIES/LABEL=${SPECIES}/${LABEL}" ;;
esac

EDTA_DIR="${EDTA_DIR:-${EDTA_INPUT_ROOT}/${SPECIES}/edta/${LABEL}}"
GENOME="${EDTA_DIR}/${SPECIES}_${LABEL}.fa.mod.MAKER.masked"
TE_GFF="${EDTA_DIR}/${SPECIES}_${LABEL}.fa.mod.EDTA.TEanno.gff3"
HELIXER_GFF="${GENE_ROOT}/${SPECIES}/${LABEL}_helixer/${SPECIES}_${LABEL}.edta_masked.helixer.gff3"
OUT_DIR="${GENE_ROOT}/${SPECIES}/${LABEL}_eviann"
POST_DIR="${OUT_DIR}/postprocess"
PREFIX="${SPECIES}_${LABEL}.eviann"
require_file "${GENOME}" "EDTA masked genome"
require_file "${TE_GFF}" "EDTA TE GFF3"
require_file "${HELIXER_GFF}" "Helixer GFF3"
require_file "${POSTPROCESS_SCRIPT}" "postprocess helper"
require_file "${CONDA_SH}" "conda initialization"

mkdir -p "${EVIDENCE_ROOT}" "${OUT_DIR}" "${POST_DIR}" "${ROOT_DIR}/logs"
RELATED_TRANSCRIPTS="${EVIDENCE_ROOT}/related_species_transcripts.fa"
RELATED_PROTEINS="${EVIDENCE_ROOT}/related_species_proteins.fa"

build_evidence() {
    if [[ ! -s "${RELATED_TRANSCRIPTS}" ]]; then
        : > "${RELATED_TRANSCRIPTS}"
        for source in \
            "${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_cds_from_genomic.fna.gz" \
            "${ROOT_DIR}/data/pcam/GCA_050084535.1_PCAv2.0_cds_from_genomic.fna.gz" \
            "${ROOT_DIR}/data/pmu/GCF_000346735.1_P.mume_V1.0_cds_from_genomic.fna.gz" \
            "${ROOT_DIR}/data/pspe/GCA_041154625.1_CerSpe_IZO01_v1.0_cds_from_genomic.fna.gz"; do
            [[ -s "${source}" ]] && gzip -dc "${source}" >> "${RELATED_TRANSCRIPTS}"
        done
    fi
    if [[ ! -s "${RELATED_PROTEINS}" ]]; then
        : > "${RELATED_PROTEINS}"
        for source in \
            "${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_protein.faa.gz" \
            "${ROOT_DIR}/data/pcam/GCA_050084535.1_PCAv2.0_protein.faa.gz" \
            "${ROOT_DIR}/data/pmu/GCF_000346735.1_P.mume_V1.0_protein.faa.gz" \
            "${ROOT_DIR}/data/pspe/GCA_041154625.1_CerSpe_IZO01_v1.0_protein.faa.gz"; do
            [[ -s "${source}" ]] && gzip -dc "${source}" >> "${RELATED_PROTEINS}"
        done
    fi
    require_file "${RELATED_TRANSCRIPTS}" "related CDS evidence"
    require_file "${RELATED_PROTEINS}" "related protein evidence"
}

source "${CONDA_SH}"
conda activate "${EVIANN_ENV}"
build_evidence
cd "${OUT_DIR}"
RAW_GFF="${OUT_DIR}/${SPECIES}_${LABEL}.fa.pseudo_label.gff"
if [[ ! -s "${RAW_GFF}" ]]; then
    if [[ "${SPECIES}" == "ppa" ]]; then
        # PPA structural annotation is kept separate from functional
        # annotation; do not invoke EviAnn's -f/my_maker_functional_gff.
        cmd=(eviann.sh -g "${GENOME}" -t "${THREADS}" -e "${RELATED_TRANSCRIPTS}" -p "${RELATED_PROTEINS}" -c "${HELIXER_GFF}" --verbose)
    else
        cmd=(eviann.sh -g "${GENOME}" -t "${THREADS}" -e "${RELATED_TRANSCRIPTS}" -p "${RELATED_PROTEINS}" -c "${HELIXER_GFF}" -f --verbose)
        [[ -s "${SWISSPROT}" ]] && cmd+=(-s "${SWISSPROT}")
    fi
    "${cmd[@]}"
else
    echo "[SKIP] Existing EviAnn pseudo-label GFF: ${RAW_GFF}"
fi
require_file "${RAW_GFF}" "EviAnn pseudo-label GFF"

FUNCTIONAL_GFF="${OUT_DIR}/${SPECIES}_${LABEL}.fa.functional_note.pseudo_label.gff"
POST_GFF="${RAW_GFF}"
[[ "${SPECIES}" != "ppa" && -s "${FUNCTIONAL_GFF}" ]] && POST_GFF="${FUNCTIONAL_GFF}"
source "${CONDA_SH}"
conda activate "${PY_ENV}"
python "${POSTPROCESS_SCRIPT}" --genome "${GENOME}" --gff3 "${POST_GFF}" --te-gff3 "${TE_GFF}" --out-dir "${POST_DIR}" --prefix "${PREFIX}" --short-cds 300 --gene-te-fraction 0.80 --cds-te-fraction 0.50

for output in "${POST_DIR}/${PREFIX}.genes.gff3" "${POST_DIR}/${PREFIX}.cds.fa" "${POST_DIR}/${PREFIX}.protein.fa" "${POST_DIR}/${PREFIX}.gene_model_qc.tsv"; do
    require_file "${output}" "EviAnn postprocess output"
done
echo "[OK] EviAnn output: ${POST_DIR}"
