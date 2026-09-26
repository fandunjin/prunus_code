#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
WORK_DIR="${OUT_DIR}/work"

PZI_FASTA="${ROOT_DIR}/final_assembly/pzi_2n32.chromosomes.fa"
PCAM_FASTA="${ROOT_DIR}/data/pcam/GCA_050084535.1_PCAv2.0_genomic.fna.gz"
PPA_FASTA="${ROOT_DIR}/data/ppa/GCA_964059315.1_drPruPadu1.1_genomic.fna.gz"

PZI_VS_PCAM_PAF="${ROOT_DIR}/haploid/work/pzi_vs_pcam.asm5.paf"
PZI_SELF_PAF="${ROOT_DIR}/haploid/work/pzi_chr32_self.asm5.paf"
PPA_VS_PZI_PAF="${WORK_DIR}/ppa_vs_pzi.asm5.paf"

CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
MINIMAP_ENV="minimap"
HELPER="${ROOT_DIR}/scripts/trispecies_rearrangement_review.py"
DOTPLOT_SCRIPT="${ROOT_DIR}/scripts/plot_minimap_synteny_dotplot.py"
HEATMAP_SCRIPT="${ROOT_DIR}/scripts/plot_minimap_synteny_heatmap.py"

THREADS=24
MINIMAP_PRESET="asm10"

require_file() {
    local file_path="$1"
    local label="$2"
    if [ ! -s "${file_path}" ]; then
        echo "ERROR: ${label} does not exist or is empty: ${file_path}" >&2
        exit 1
    fi
}

mkdir -p "${OUT_DIR}" "${WORK_DIR}" "${ROOT_DIR}/logs"

require_file "${PZI_FASTA}" "PZI FASTA"
require_file "${PCAM_FASTA}" "P. campanulata FASTA"
require_file "${PPA_FASTA}" "P. persica FASTA"
require_file "${PZI_VS_PCAM_PAF}" "PZI vs PCAM PAF"
require_file "${PZI_SELF_PAF}" "PZI self PAF"
require_file "${CONDA_SH}" "Conda init script"
require_file "${HELPER}" "tri-species review helper"
require_file "${DOTPLOT_SCRIPT}" "dotplot script"
require_file "${HEATMAP_SCRIPT}" "heatmap script"

source "${CONDA_SH}"
conda activate "${MINIMAP_ENV}"

if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    PYTHON_BIN="python"
fi

if [ ! -s "${PPA_VS_PZI_PAF}" ]; then
    minimap2 -x "${MINIMAP_PRESET}" -t "${THREADS}" "${PZI_FASTA}" "${PPA_FASTA}" > "${PPA_VS_PZI_PAF}"
fi

"${PYTHON_BIN}" "${HELPER}" \
    --pzi-vs-pcam-paf "${PZI_VS_PCAM_PAF}" \
    --pzi-self-paf "${PZI_SELF_PAF}" \
    --ppa-vs-pzi-paf "${PPA_VS_PZI_PAF}" \
    --pzi-fasta "${PZI_FASTA}" \
    --pcam-fasta "${PCAM_FASTA}" \
    --ppa-fasta "${PPA_FASTA}" \
    --ppa-main-count 16 \
    --out-dir "${OUT_DIR}"

"${PYTHON_BIN}" "${DOTPLOT_SCRIPT}" \
    --paf "${PZI_VS_PCAM_PAF}" \
    --out-svg "${OUT_DIR}/pzi_vs_pcam.colorband.dotplot.svg" \
    --title "P. zippeliana vs P. campanulata with axis bands" \
    --axis-colors "${OUT_DIR}/axis_colors.tsv"

"${PYTHON_BIN}" "${HEATMAP_SCRIPT}" \
    --summary "${OUT_DIR}/pzi_vs_pcam_chr_alignment.tsv" \
    --out-svg "${OUT_DIR}/pzi_vs_pcam.colorband.heatmap.svg" \
    --title "P. zippeliana vs P. campanulata with axis bands" \
    --x-label "P. zippeliana chromosomes" \
    --y-label "P. campanulata chromosomes" \
    --axis-colors "${OUT_DIR}/axis_colors.tsv"

cat > "${OUT_DIR}/README_TRI_SPECIES_REARRANGEMENT.md" <<'EOF'
# Trispecies rearrangement review

This folder contains chromosome-level review tables for:

- P. zippeliana
- P. campanulata
- P. persica

The current output is for review only.
It adds axis color bands to synteny plots and classifies each P. campanulata
chromosome into provisional rearrangement categories:

- balanced_quartet
- possible_break_or_split
- possible_fusion_or_nested_fusion
- complex_rearrangement

The labels are provisional and must be checked against the dotplot, heatmap,
and the P. padus corroboration table.
EOF

echo "Wrote tri-species rearrangement review into ${OUT_DIR}"
