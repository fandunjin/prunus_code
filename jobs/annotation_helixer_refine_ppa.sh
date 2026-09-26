#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${PROJECT_ROOT}/data/ppa"
ANNOTATION_DIR="${DATA_DIR}/annotation"
PPA_EDTA_DIR="${PPA_EDTA_DIR:-${ROOT_DIR}/../prunus_old/analysis/prunus_pan_te/discovery/ppa}"
PPA_EDTA_TE_GFF3="${PPA_EDTA_TE_GFF3:-${PPA_EDTA_DIR}/ppa.fa.mod.EDTA.TEanno.gff3}"
# Use the unmasked PPA assembly as EviAnn's genome input. The EDTA TE GFF3 is
# retained for postprocess filtering and has matching coordinates/sequence IDs.
GENOME="${GENOME:-${PPA_EDTA_DIR}/ppa.fa}"
TE_GFF3="${TE_GFF3:-${PPA_EDTA_TE_GFF3}}"
HELIXER_GFF="${ANNOTATION_DIR}/helixer/ppa.helixer.gff3"
REFINE_DIR="${ANNOTATION_DIR}/eviann"
EVIDENCE_DIR="${REFINE_DIR}/evidence"
POST_DIR="${REFINE_DIR}/postprocess"
# Isolate EviAnn intermediates from previous PPA masked-genome runs.
EVIANN_RUN_DIR="${EVIANN_RUN_DIR:-${REFINE_DIR}/run_ppa_unmasked}"
PREFIX="ppa.eviann"

# Use both completed genome-guided transcript sets as the PPA external
# transcript evidence. Protein evidence is restricted to PPE and PCAM so that
# PZI, PSP and PPA use the same comparative protein database.
RELATED_TRANSCRIPTS="${EVIDENCE_DIR}/psp_pzi.transcripts.fa"
RELATED_PROTEINS="${EVIDENCE_DIR}/ppe_pcam.proteins.fa"
USE_PROTEINS="${USE_PROTEINS:-1}"
USE_SWISSPROT="${USE_SWISSPROT:-0}"
LIFTOVER_MODE="${LIFTOVER_MODE:-1}"
# PPA is diploid (2x=16); keep EviAnn's ploidy fixed at 2.
PLOIDY=2
RELATED_TRANSCRIPT_SOURCES=(
    "${ROOT_DIR}/data/psp/annotation/genes/psp32_eviann/evidence/psp_illumina.transcripts.fa"
    "${ROOT_DIR}/annotation/genes/pzi_eviann/evidence/pzi_illumina_pacbio.transcripts.fa"
)
RELATED_PROTEIN_SOURCES=(
    "${ROOT_DIR}/data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_protein.faa.gz"
    "${ROOT_DIR}/data/pcam/GCA_050084535.1_PCAv2.0_protein.faa.gz"
)
RELATED_PROTEIN_LABELS=(ppe pcam)

CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
EVIANN_ENV="${EVIANN_ENV:-eviann}"
PY_ENV="${PY_ENV:-haphic}"
THREADS="${THREADS:-24}"
SWISSPROT="${SWISSPROT:-${HOME}/CX/database/SwissProt/uniprot_sprot.fasta}"
POSTPROCESS_SCRIPT="${ROOT_DIR}/scripts/helixer_postprocess_qc.py"
GFF_STATS_SCRIPT="${ROOT_DIR}/scripts/gff3_basic_stats.py"

require_file() {
    local path="$1"
    local label="$2"
    if [ ! -s "${path}" ]; then
        echo "ERROR: ${label} does not exist or is empty: ${path}" >&2
        exit 1
    fi
}

resolve_eviann_gff() {
    [[ -s "${RAW_GFF}" ]] && return 0
    local candidates=()
    while IFS= read -r -d '' candidate; do
        candidates+=("${candidate}")
    done < <(find "${EVIANN_RUN_DIR}" -maxdepth 1 -type f -name '*.pseudo_label.gff' ! -name '*functional_note*' -size +0c -print0)
    if [[ "${#candidates[@]}" -eq 1 ]]; then
        RAW_GFF="${candidates[0]}"
        echo "WARNING: EviAnn used a nonstandard output basename; using ${RAW_GFF}"
        return 0
    fi
    echo "EviAnn output files in ${EVIANN_RUN_DIR}:" >&2
    find "${EVIANN_RUN_DIR}" -maxdepth 1 -type f -printf '  %f\t%s bytes\n' | sort >&2 || true
    return 1
}

append_fasta() {
    local src="$1" dst="$2"
    [[ -s "$src" ]] || return 1
    if [[ "$src" == *.gz ]]; then gzip -dc "$src"; else cat "$src"; fi \
        | awk -v prefix="$3" '/^>/{sub(/^>/, ">" prefix "|")} {print}' >> "$dst"
}

require_file "${GENOME}" "P. padus unmasked genome"
require_file "${TE_GFF3}" "P. padus EDTA TE annotation GFF3"
require_file "${HELIXER_GFF}" "P. padus Helixer GFF3"
require_file "${CONDA_SH}" "Conda initialization script"
require_file "${POSTPROCESS_SCRIPT}" "postprocess/QC script"
require_file "${GFF_STATS_SCRIPT}" "GFF3 statistics script"
mkdir -p "${EVIDENCE_DIR}" "${POST_DIR}" "${EVIANN_RUN_DIR}"

case "${USE_PROTEINS}" in 0|1) ;; *) echo "ERROR: USE_PROTEINS must be 0 or 1" >&2; exit 1 ;; esac
case "${USE_SWISSPROT}" in 0|1) ;; *) echo "ERROR: USE_SWISSPROT must be 0 or 1" >&2; exit 1 ;; esac
case "${LIFTOVER_MODE}" in 0|1) ;; *) echo "ERROR: LIFTOVER_MODE must be 0 or 1" >&2; exit 1 ;; esac
[[ "${PLOIDY}" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: PLOIDY must be a positive integer" >&2; exit 1; }
: > "${RELATED_TRANSCRIPTS}.tmp"
require_file "${RELATED_TRANSCRIPT_SOURCES[0]}" "PSP transcript evidence"
require_file "${RELATED_TRANSCRIPT_SOURCES[1]}" "PZI transcript evidence"
append_fasta "${RELATED_TRANSCRIPT_SOURCES[0]}" "${RELATED_TRANSCRIPTS}.tmp" psp
append_fasta "${RELATED_TRANSCRIPT_SOURCES[1]}" "${RELATED_TRANSCRIPTS}.tmp" pzi
mv -f "${RELATED_TRANSCRIPTS}.tmp" "${RELATED_TRANSCRIPTS}"
require_file "${RELATED_TRANSCRIPTS}" "combined PSP+PZI transcript evidence"
if [[ "${USE_PROTEINS}" == 1 ]]; then
    : > "${RELATED_PROTEINS}"
    for i in "${!RELATED_PROTEIN_SOURCES[@]}"; do
        source="${RELATED_PROTEIN_SOURCES[$i]}"
        append_fasta "$source" "${RELATED_PROTEINS}" "${RELATED_PROTEIN_LABELS[$i]}" \
            || { echo "ERROR: missing protein source: $source" >&2; exit 1; }
    done
fi
if [[ "${USE_PROTEINS}" == 1 ]]; then require_file "${RELATED_PROTEINS}" "PPE+PCAM protein evidence"; fi

source "${CONDA_SH}"
conda activate "${EVIANN_ENV}"
command -v eviann.sh >/dev/null 2>&1 || {
    echo "ERROR: eviann.sh not found in ${EVIANN_ENV}" >&2
    exit 1
}

# Some installed EviAnn releases do not implement -l.  Detect support rather
# than passing an option that would abort the complete PPA annotation.
LIFTOVER_USED=0
if [[ "${LIFTOVER_MODE}" == 1 ]] && eviann.sh --help 2>&1 | grep -Eq '(^|[[:space:]])-l([[:space:]]|$)'; then
    LIFTOVER_USED=1
fi

GENOME_BASENAME="$(basename "${GENOME}")"
RAW_GFF="${EVIANN_RUN_DIR}/${GENOME_BASENAME}.pseudo_label.gff"
FUNCTIONAL_GFF="${EVIANN_RUN_DIR}/${GENOME_BASENAME}.functional_note.pseudo_label.gff"
if [ ! -s "${RAW_GFF}" ]; then
    cd "${EVIANN_RUN_DIR}"
    EVIANN_CMD=(
        eviann.sh
        -g "${GENOME}"
        -t "${THREADS}"
        -d "${PLOIDY}"
        -e "${RELATED_TRANSCRIPTS}"
        -c "${HELIXER_GFF}"
        -f
        --verbose
    )
    if [[ "${LIFTOVER_USED}" == 1 ]]; then
        EVIANN_CMD+=(-l)
    fi
    if [[ "${USE_PROTEINS}" == 1 ]]; then
        EVIANN_CMD+=(-p "${RELATED_PROTEINS}")
    fi
    if [[ "${USE_PROTEINS}" == 1 && "${USE_SWISSPROT}" == 1 && -s "${SWISSPROT}" ]]; then
        EVIANN_CMD+=(-s "${SWISSPROT}")
    fi
    "${EVIANN_CMD[@]}"
else
    echo "Existing EviAnn GFF found; skipping refinement: ${RAW_GFF}"
fi
resolve_eviann_gff || { echo "ERROR: EviAnn did not produce a non-empty PPA pseudo_label GFF3 in ${EVIANN_RUN_DIR}" >&2; exit 1; }
require_file "${RAW_GFF}" "EviAnn GFF3"

FINAL_GFF="${RAW_GFF}"
if [ -s "${FUNCTIONAL_GFF}" ]; then
    FINAL_GFF="${FUNCTIONAL_GFF}"
fi

conda activate "${PY_ENV}"
if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi
POSTPROCESS_CMD=(
    "${PYTHON_BIN}" "${POSTPROCESS_SCRIPT}"
    --genome "${GENOME}"
    --gff3 "${FINAL_GFF}"
    --out-dir "${POST_DIR}"
    --prefix "${PREFIX}"
    --short-cds 300
    --gene-te-fraction 0.80
    --cds-te-fraction 0.50
    --te-gff3 "${TE_GFF3}"
)
"${POSTPROCESS_CMD[@]}"

require_file "${POST_DIR}/${PREFIX}.genes.gff3" "final P. padus GFF3"
require_file "${POST_DIR}/${PREFIX}.cds.fa" "final P. padus CDS"
require_file "${POST_DIR}/${PREFIX}.protein.fa" "final P. padus protein"
"${PYTHON_BIN}" "${GFF_STATS_SCRIPT}" --gff3 "${POST_DIR}/${PREFIX}.genes.gff3" --prefix "${PREFIX}" --out-dir "${POST_DIR}/stats"

echo "Final GFF3   : ${POST_DIR}/${PREFIX}.genes.gff3"
echo "Final CDS    : ${POST_DIR}/${PREFIX}.cds.fa"
echo "Final protein: ${POST_DIR}/${PREFIX}.protein.fa"
