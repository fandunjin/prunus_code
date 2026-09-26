#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

# Align the unmodified PPA primary assembly against the unmodified PZI copy1
# assembly. PPA chromosome IDs are recorded only to filter and summarize PAF.
ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/analysis/figure2/01_ppa_pzi_copy1_minimap}"
PPA_PRIMARY_GENOME="${PPA_PRIMARY_GENOME:-${ROOT_DIR}/data/ppa/GCA_964059315.1_drPruPadu1.1_genomic.fna.gz}"
PPA_RENAME_MAP="${PPA_RENAME_MAP:-${WORK_DIR}/ppa_primary_chr16.rename_map.tsv}"
PZI_COPY1_FASTA="${PZI_COPY1_FASTA:-${ROOT_DIR}/final_genome/pzi_copy1.fa}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
MINIMAP_ENV="${MINIMAP_ENV:-minimap}"
DOTPLOT_SCRIPT="${DOTPLOT_SCRIPT:-${ROOT_DIR}/scripts/plot_minimap_synteny_dotplot.py}"
HEATMAP_SCRIPT="${HEATMAP_SCRIPT:-${ROOT_DIR}/scripts/plot_minimap_synteny_heatmap.py}"
PREPARE_PPA_SCRIPT="${PREPARE_PPA_SCRIPT:-${ROOT_DIR}/scripts/prepare_ppa_primary_chr16.py}"
SUMMARIZE_PAF_SCRIPT="${SUMMARIZE_PAF_SCRIPT:-${ROOT_DIR}/scripts/summarize_minimap_paf_pairs.py}"
THREADS="${THREADS:-24}"
MINIMAP_PRESET="${MINIMAP_PRESET:-asm5}"
MIN_MAPQ="${MIN_MAPQ:-20}"
MIN_ALN_LEN="${MIN_ALN_LEN:-10000}"
# Default to a clean rerun because earlier outputs used SyRI-derived inputs.
# Set FORCE_RERUN=0 only after this original-input result has completed.
FORCE_RERUN="${FORCE_RERUN:-1}"

require_file() {
    local path="$1"
    local label="$2"
    [[ -s "${path}" ]] || { echo "[ERROR] Missing ${label}: ${path}" >&2; exit 2; }
}

remove_if_forced() {
    if [[ "${FORCE_RERUN}" == "1" ]]; then
        rm -f "$@"
    fi
}

require_file "${PPA_PRIMARY_GENOME}" "PPA annotated primary genome"
require_file "${PZI_COPY1_FASTA}" "PZI raw copy1 FASTA"
require_file "${CONDA_SH}" "conda initialization script"
require_file "${DOTPLOT_SCRIPT}" "dotplot script"
require_file "${HEATMAP_SCRIPT}" "heatmap script"
require_file "${PREPARE_PPA_SCRIPT}" "PPA primary chromosome extractor"
require_file "${SUMMARIZE_PAF_SCRIPT}" "primary PAF coverage summarizer"

mkdir -p "${WORK_DIR}" "${ROOT_DIR}/logs"

source "${CONDA_SH}"
conda activate "${MINIMAP_ENV}"

command -v minimap2 >/dev/null 2>&1 || { echo "[ERROR] minimap2 is unavailable in ${MINIMAP_ENV}" >&2; exit 2; }
command -v python3 >/dev/null 2>&1 && PYTHON_BIN=python3 || PYTHON_BIN=python

if [[ ! -s "${PPA_RENAME_MAP}" || "${FORCE_RERUN}" == "1" ]]; then
    "${PYTHON_BIN}" "${PREPARE_PPA_SCRIPT}" \
        --genome "${PPA_PRIMARY_GENOME}" \
        --out-map "${PPA_RENAME_MAP}"
fi
require_file "${PPA_RENAME_MAP}" "PPA primary chromosome ID map"

prefix="ppa_vs_pzi_copy1_chr16"
paf="${WORK_DIR}/${prefix}.${MINIMAP_PRESET}.paf"
filtered_paf="${WORK_DIR}/${prefix}.${MINIMAP_PRESET}.mapq${MIN_MAPQ}.min${MIN_ALN_LEN}.paf"
summary="${WORK_DIR}/${prefix}.${MINIMAP_PRESET}.alignment_summary.tsv"
best="${WORK_DIR}/${prefix}.${MINIMAP_PRESET}.best_pzi_copy1_hit_per_ppa_superscaffold.tsv"
dotplot="${WORK_DIR}/${prefix}.${MINIMAP_PRESET}.alignment_dotplot.svg"
heatmap="${WORK_DIR}/${prefix}.${MINIMAP_PRESET}.alignment_heatmap.svg"
readme="${WORK_DIR}/README.txt"

remove_if_forced "${paf}" "${filtered_paf}" "${summary}" "${best}" "${dotplot}" "${heatmap}" "${readme}"

echo "=== Figure 2: PPA vs PZI copy1 minimap2 macro-synteny ==="
echo "PPA target (raw primary assembly): ${PPA_PRIMARY_GENOME}"
echo "PZI query  (raw copy1 assembly) : ${PZI_COPY1_FASTA}"
echo "Output directory : ${WORK_DIR}"
echo "Preset           : ${MINIMAP_PRESET}"
echo "Minimum MAPQ     : ${MIN_MAPQ}"
echo "Minimum alignment: ${MIN_ALN_LEN} bp"

if [[ ! -s "${paf}" ]]; then
    minimap2 -x "${MINIMAP_PRESET}" -t "${THREADS}" "${PPA_PRIMARY_GENOME}" "${PZI_COPY1_FASTA}" > "${paf}"
else
    echo "[INFO] Reusing existing PAF: ${paf}"
fi

"${PYTHON_BIN}" "${SUMMARIZE_PAF_SCRIPT}" \
    --paf "${paf}" \
    --ppa-map "${PPA_RENAME_MAP}" \
    --out-filtered-paf "${filtered_paf}" \
    --out-summary "${summary}" \
    --out-best "${best}" \
    --min-mapq "${MIN_MAPQ}" \
    --min-alignment-len "${MIN_ALN_LEN}"
require_file "${filtered_paf}" "filtered minimap2 PAF"

"${PYTHON_BIN}" "${DOTPLOT_SCRIPT}" \
    --paf "${filtered_paf}" \
    --out-svg "${dotplot}" \
    --title "P. padus primary assembly vs P. zippeliana copy1 minimap2 alignment" \
    --min-alignment-len "${MIN_ALN_LEN}" \
    --min-query-aggregate 500000 \
    --min-target-aggregate 500000

"${PYTHON_BIN}" "${HEATMAP_SCRIPT}" \
    --summary "${summary}" \
    --out-svg "${heatmap}" \
    --title "P. padus primary chromosomes vs P. zippeliana copy1 primary-chain coverage" \
    --x-label "PZI copy1 chromosomes" \
    --y-label "PPA copy1 superscaffolds"

cat > "${readme}" <<EOF
Figure 2 PPA-PZI copy1 minimap2 macro-synteny
===============================================

Inputs
------
- PPA target, unmodified annotated primary genome: ${PPA_PRIMARY_GENOME}
- PZI query, unmodified copy1 genome: ${PZI_COPY1_FASTA}
- PPA chromosome-ID map used only to retain 16 primary chromosomes: ${PPA_RENAME_MAP}

Filters
-------
- minimap2 preset: ${MINIMAP_PRESET}
- MAPQ >= ${MIN_MAPQ}
- alignment length >= ${MIN_ALN_LEN} bp

Outputs
-------
- ${paf}: raw minimap2 alignments.
- ${filtered_paf}: filtered alignments used for plotting and summaries.
- ${summary}: PPA-by-PZI non-overlapping primary-chain coverage matrix.
- ${best}: best and second-best PZI copy1 matches for each PPA primary chromosome.
- ${dotplot}: coordinate-based macro-synteny dotplot.
- ${heatmap}: PPA-by-PZI chromosome-pair heatmap.

Interpretation boundary
-----------------------
This workflow aligns the two original genome files; no SyRI-derived or
sequence-renamed FASTA is used as a minimap2 input. It identifies structural
correspondence between PPA primary chromosomes and PZI copy1 chromosomes. It
does not infer direct parental origin without gene-level synteny and
multi-gene phylogenetic confirmation.
EOF

echo "[OK] PAF: ${paf}"
echo "[OK] Summary: ${summary}"
echo "[OK] Best hits: ${best}"
echo "[OK] Dotplot: ${dotplot}"
echo "[OK] Heatmap: ${heatmap}"
