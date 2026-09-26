#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
RM_ROOT="${RM_ROOT:-${PAN_ROOT}/mask}"
STATS_SCRIPT="${ROOT_DIR}/scripts/prunus_pan_te_followup_stats.py"
DENSITY_SCRIPT="${ROOT_DIR}/scripts/prunus_pan_te_te_trf_density.py"
PLOT_SCRIPT="${ROOT_DIR}/scripts/plot_prunus_pan_te_followup.py"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
PYTHON_ENV="${PYTHON_ENV:-haphic}"
PYTHON="${PYTHON:-python}"
OUT="${PAN_ROOT}/followup_v2"
FIG="${PAN_ROOT}/figures_followup_v2"

require_file() {
    [ -s "$1" ] || { echo "ERROR: missing or empty $2: $1" >&2; exit 1; }
}

for path in \
    "${STATS_SCRIPT}" "${DENSITY_SCRIPT}" "${PLOT_SCRIPT}" "${CONDA_SH}" \
    "${PAN_ROOT}/postprocess_v2/te_driver_families.tsv" \
    "${PAN_ROOT}/postprocess_v2/te_repeatmasker_divergence.tsv" \
    "${PAN_ROOT}/postprocess_v2/te_intact_ltr.tsv" \
    "${PAN_ROOT}/followup_sololtr_v2/sololtr_summary/sololtr_by_target_lineage.tsv" \
    "${PAN_ROOT}/figure1_tandem/score50/tandem_top_arrays.tsv" \
    "${PAN_ROOT}/followup_orthosite_asm20_flank5k/orthosite/combined/orthosite_site_evidence.tsv" \
    "${PAN_ROOT}/figure1_tandem/score50/tandem_window_density.tsv" \
    "${PAN_ROOT}/postprocess_v2/te_composition.tsv" \
    "${ROOT_DIR}/downstream_analysis/phylogenetics/09_mcmctree_species14_plant_outgroup_fossil55_35_r2/03_mcmctree_chain1/species14_fossil55_racemose35.mcmctree.out" \
    "${ROOT_DIR}/downstream_analysis/phylogenetics/09_mcmctree_subgenome17_outgroup_fossil55_35_astral_r2/03_mcmctree_chain1/subgenome17_fossil55_racemose35.mcmctree.out"; do
    require_file "${path}" "required input"
done

[ "$(find "${RM_ROOT}" -type f -name '*.out' -print -quit)" ] || {
    echo "ERROR: no RepeatMasker .out files found below ${RM_ROOT}" >&2
    exit 1
}

source "${CONDA_SH}"
conda activate "${PYTHON_ENV}"
export MPLBACKEND=Agg
mkdir -p "${OUT}" "${FIG}"

echo "=== Follow-up statistics ==="
"${PYTHON}" "${STATS_SCRIPT}" \
    --drivers "${PAN_ROOT}/postprocess_v2/te_driver_families.tsv" \
    --divergence "${PAN_ROOT}/postprocess_v2/te_repeatmasker_divergence.tsv" \
    --intact-ltr "${PAN_ROOT}/postprocess_v2/te_intact_ltr.tsv" \
    --solo-by-target-lineage "${PAN_ROOT}/followup_sololtr_v2/sololtr_summary/sololtr_by_target_lineage.tsv" \
    --arrays "${PAN_ROOT}/figure1_tandem/score50/tandem_top_arrays.tsv" \
    --orthosite "${PAN_ROOT}/followup_orthosite_asm20_flank5k/orthosite/combined/orthosite_site_evidence.tsv" \
    --repeatmasker-root "${RM_ROOT}" \
    --species-mcmctree "${ROOT_DIR}/downstream_analysis/phylogenetics/09_mcmctree_species14_plant_outgroup_fossil55_35_r2/03_mcmctree_chain1/species14_fossil55_racemose35.mcmctree.out" \
    --subgenome-mcmctree "${ROOT_DIR}/downstream_analysis/phylogenetics/09_mcmctree_subgenome17_outgroup_fossil55_35_astral_r2/03_mcmctree_chain1/subgenome17_fossil55_racemose35.mcmctree.out" \
    --mu-time-unit-mya 100 \
    --output-dir "${OUT}"

echo "=== 1-Mb TE/TRF density ==="
"${PYTHON}" "${DENSITY_SCRIPT}" \
    --repeatmasker-root "${RM_ROOT}" \
    --trf-density "${PAN_ROOT}/figure1_tandem/score50/tandem_window_density.tsv" \
    --te-composition "${PAN_ROOT}/postprocess_v2/te_composition.tsv" \
    --output-dir "${OUT}"

echo "=== Follow-up figures ==="
"${PYTHON}" "${PLOT_SCRIPT}" --stats-dir "${OUT}" --output-dir "${FIG}"

echo "=== TE/TRF follow-up completed ==="
