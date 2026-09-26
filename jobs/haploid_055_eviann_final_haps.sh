#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
EVIANN_ENV="${EVIANN_ENV:-eviann}"
PY_ENV="${PY_ENV:-haphic}"
FINAL_HAP_DIR="${ROOT_DIR}/haploid/final_haplotypes"
HELIXER_ROOT="${ROOT_DIR}/annotation/genes"
EVIANN_ROOT="${ROOT_DIR}/annotation/genes/final_haplotypes_eviann"
EVIANN_PACKAGE="${ROOT_DIR}/final_haplotype_genomes_eviann"
COLLECT_SCRIPT="${ROOT_DIR}/scripts/collect_eviann_haplotype_genomes.py"
POSTPROCESS_SCRIPT="${ROOT_DIR}/scripts/helixer_postprocess_qc.py"

# No PZI RNA/transcript evidence is available. Following the EviAnn evidence
# route, use related-species CDS/transcripts and proteins as external evidence.
RELATED_TRANSCRIPTS="${RELATED_TRANSCRIPTS:-${EVIANN_ROOT}/evidence/related_species_transcripts.fa}"
RELATED_PROTEINS="${RELATED_PROTEINS:-${EVIANN_ROOT}/evidence/related_species_proteins.fa}"
SWISSPROT="${SWISSPROT:-${HOME}/CX/database/SwissProt/uniprot_sprot.fasta}"
RUN_COLLECT="${RUN_COLLECT:-1}"

TRANSCRIPT_EVIDENCE_SOURCES=(
    "${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_cds_from_genomic.fna.gz"
    "${ROOT_DIR}/data/pcam/GCA_050084535.1_PCAv2.0_cds_from_genomic.fna.gz"
    "${ROOT_DIR}/data/pmu/GCF_000346735.1_P.mume_V1.0_cds_from_genomic.fna.gz"
    "${ROOT_DIR}/data/pspe/GCA_041154625.1_CerSpe_IZO01_v1.0_cds_from_genomic.fna.gz"
)
PROTEIN_EVIDENCE_SOURCES=(
    "${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_protein.faa.gz"
    "${ROOT_DIR}/data/pcam/GCA_050084535.1_PCAv2.0_protein.faa.gz"
    "${ROOT_DIR}/data/pmu/GCF_000346735.1_P.mume_V1.0_protein.faa.gz"
    "${ROOT_DIR}/data/pspe/GCA_041154625.1_CerSpe_IZO01_v1.0_protein.faa.gz"
)

require_file() {
    local file_path="$1"
    local label="$2"
    if [ ! -s "${file_path}" ]; then
        echo "ERROR: ${label} does not exist or is empty: ${file_path}" >&2
        exit 1
    fi
}

append_fasta_if_exists() {
    local src="$1"
    local dst="$2"
    if [ -s "${src}" ]; then
        if [[ "${src}" == *.gz ]]; then
            gzip -dc "${src}" >> "${dst}"
        else
            cat "${src}" >> "${dst}"
        fi
        printf "%s\n" "${src}" >> "${dst}.sources.txt"
    fi
}

build_related_evidence() {
    mkdir -p "${EVIANN_ROOT}/evidence"
    if [ "${FORCE}" = "1" ]; then
        rm -f "${RELATED_TRANSCRIPTS}" "${RELATED_TRANSCRIPTS}.sources.txt" "${RELATED_PROTEINS}" "${RELATED_PROTEINS}.sources.txt"
    fi
    if [ ! -s "${RELATED_TRANSCRIPTS}" ]; then
        : > "${RELATED_TRANSCRIPTS}"
        : > "${RELATED_TRANSCRIPTS}.sources.txt"
        for src in "${TRANSCRIPT_EVIDENCE_SOURCES[@]}"; do
            append_fasta_if_exists "${src}" "${RELATED_TRANSCRIPTS}"
        done
    fi
    if [ ! -s "${RELATED_PROTEINS}" ]; then
        : > "${RELATED_PROTEINS}"
        : > "${RELATED_PROTEINS}.sources.txt"
        for src in "${PROTEIN_EVIDENCE_SOURCES[@]}"; do
            append_fasta_if_exists "${src}" "${RELATED_PROTEINS}"
        done
    fi
    require_file "${RELATED_TRANSCRIPTS}" "related-species transcript/CDS evidence for EviAnn -e"
    require_file "${RELATED_PROTEINS}" "related-species protein evidence for EviAnn -p"
}

run_one_hap() {
    local hap="$1"
    local genome="${FINAL_HAP_DIR}/pzi_${hap}.fa"
    local helixer_dir="${HELIXER_ROOT}/final_${hap}_helixer"
    local helixer_gff="${helixer_dir}/pzi_${hap}.helixer.gff3"
    local out_dir="${EVIANN_ROOT}/${hap}"
    local raw_gff="${out_dir}/pzi_${hap}.fa.pseudo_label.gff"
    local functional_gff="${out_dir}/pzi_${hap}.fa.functional_note.pseudo_label.gff"
    local post_dir="${out_dir}/postprocess"
    local prefix="pzi_${hap}.eviann"

    require_file "${genome}" "${hap} genome FASTA"
    require_file "${helixer_gff}" "${hap} Helixer GFF3"
    mkdir -p "${out_dir}" "${post_dir}"

    echo ">>> EviAnn ${hap}"
    echo "Genome             : ${genome}"
    echo "Helixer evidence   : ${helixer_gff}"
    echo "Transcript evidence: ${RELATED_TRANSCRIPTS}"
    echo "Protein evidence   : ${RELATED_PROTEINS}"
    if [ -s "${SWISSPROT}" ]; then
        echo "SwissProt          : ${SWISSPROT}"
    else
        echo "SwissProt          : not used; file not found: ${SWISSPROT}"
    fi
    echo "Out                : ${out_dir}"
    echo

    cd "${out_dir}"
    if [ ! -s "${raw_gff}" ] || [ "${FORCE}" = "1" ]; then
        EVIANN_CMD=(
            eviann.sh
            -g "${genome}"
            -t "${THREADS}"
            -e "${RELATED_TRANSCRIPTS}"
            -p "${RELATED_PROTEINS}"
            -c "${helixer_gff}"
            -f
            --verbose
        )
        if [ -s "${SWISSPROT}" ]; then
            EVIANN_CMD+=(-s "${SWISSPROT}")
        fi
        "${EVIANN_CMD[@]}"
    else
        echo "Existing EviAnn GFF found, skipping: ${raw_gff}"
    fi

    if [ ! -s "${raw_gff}" ]; then
        echo "ERROR: expected EviAnn output was not found: ${raw_gff}" >&2
        echo "Check ${out_dir}; if your EviAnn version writes another GFF name, update raw_gff in this script." >&2
        exit 1
    fi
    if [ -s "${functional_gff}" ]; then
        postprocess_gff="${functional_gff}"
    else
        postprocess_gff="${raw_gff}"
    fi
    echo "Postprocess GFF    : ${postprocess_gff}"

    source "${CONDA_SH}"
    conda activate "${PY_ENV}"
    if command -v python3 >/dev/null 2>&1; then
        PYTHON_BIN="python3"
    else
        PYTHON_BIN="python"
    fi
    "${PYTHON_BIN}" "${POSTPROCESS_SCRIPT}" \
        --genome "${genome}" \
        --gff3 "${postprocess_gff}" \
        --out-dir "${post_dir}" \
        --prefix "${prefix}" \
        --short-cds 300 \
        --gene-te-fraction 0.80 \
        --cds-te-fraction 0.50
    conda activate "${EVIANN_ENV}"

    require_file "${post_dir}/${prefix}.genes.gff3" "${hap} EviAnn postprocessed GFF3"
    require_file "${post_dir}/${prefix}.cds.fa" "${hap} EviAnn postprocessed CDS"
    require_file "${post_dir}/${prefix}.protein.fa" "${hap} EviAnn postprocessed protein"
}

case "${EVIANN_HAP_SET}" in
    all|hapA1|hapA2|hapB1|hapB2) ;;
    *)
        echo "ERROR: EVIANN_HAP_SET must be all, hapA1, hapA2, hapB1, or hapB2. Observed: ${EVIANN_HAP_SET}" >&2
        exit 1
        ;;
esac

require_file "${CONDA_SH}" "Conda initialization script"
require_file "${COLLECT_SCRIPT}" "EviAnn final haplotype collection script"
require_file "${POSTPROCESS_SCRIPT}" "postprocess/QC Python script"
mkdir -p "${EVIANN_ROOT}" "${ROOT_DIR}/logs"

source "${CONDA_SH}"
conda activate "${EVIANN_ENV}"

echo "=== EviAnn refinement for final PZI haplotypes ==="
echo "EVIANN_HAP_SET : ${EVIANN_HAP_SET}"
echo "Threads : ${THREADS}"
echo "Out root: ${EVIANN_ROOT}"
echo "Note    : this script ignores any inherited HAP_SET from Helixer/postprocess jobs."
echo

build_related_evidence

if [ "${EVIANN_HAP_SET}" = "all" ]; then
    for hap in hapA1 hapA2 hapB1 hapB2; do
        run_one_hap "${hap}"
    done
else
    run_one_hap "${EVIANN_HAP_SET}"
fi

if [ "${RUN_COLLECT}" = "1" ] && [ "${EVIANN_HAP_SET}" = "all" ]; then
    echo ">>> Collecting EviAnn-refined haplotypes into downstream package"
    conda activate "${PY_ENV}"
    if command -v python3 >/dev/null 2>&1; then
        PYTHON_BIN="python3"
    else
        PYTHON_BIN="python"
    fi
    "${PYTHON_BIN}" "${COLLECT_SCRIPT}" \
        --root-dir "${ROOT_DIR}" \
        --eviann-root "${EVIANN_ROOT}" \
        --out-dir "${EVIANN_PACKAGE}"
else
    echo "Skipping collection because RUN_COLLECT=${RUN_COLLECT} and EVIANN_HAP_SET=${EVIANN_HAP_SET}."
    echo "To collect only, run: jsub < jobs/haploid_056_collect_eviann_final_haps.sh"
fi

echo "=== EviAnn finished ==="
echo "EviAnn root: ${EVIANN_ROOT}"
echo "Package    : ${EVIANN_PACKAGE}"
