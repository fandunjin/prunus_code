#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Local WSL runner for the reference-style Fig. 2b PUG mapping.
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping}"
INPUT="${INPUT:-${BASE}/00_family_selection_alltaxa_syntenic_strict_v2}"
TREE_STAGE="${TREE_STAGE:-${BASE}/02_family_trees_reference_default_alltaxa_strict}"
OUT="${OUT:-${BASE}/03_pug_reference_default_alltaxa_strict}"
PUG_DIR="${PUG_DIR:-${ROOT_DIR}/software/PUG}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
PERL_BIN="${PERL_BIN:-/usr/bin/perl}"
PUG_SPECIES_TREE="${PUG_SPECIES_TREE:-}"
PUG_EXPECTED_LABELS="${PUG_EXPECTED_LABELS:-}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }

require_file "${INPUT}/selected_collinear_wgd_families.tsv"
[[ -e "${TREE_STAGE}/COMPLETE" ]] || die "Missing completion marker: ${TREE_STAGE}/COMPLETE"
[[ -d "${TREE_STAGE}/trees" ]] || die "Missing formal tree directory: ${TREE_STAGE}/trees"
[[ -f "${PUG_DIR}/PUG.pl" ]] || die "Missing PUG.pl: ${PUG_DIR}"
[[ -x "${PERL_BIN}" ]] || die "Perl executable unavailable: ${PERL_BIN}"
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || die "python is unavailable"
"${PERL_BIN}" -MBio::TreeIO -e 1 >/dev/null 2>&1 || die "BioPerl TreeIO unavailable in ${PERL_BIN}"
[[ ! -e "${OUT}/COMPLETE" ]] || die "Output is already complete: ${OUT}"

mkdir -p "${OUT}"
PREPARE_ARGS=(
    --families "${INPUT}/selected_collinear_wgd_families.tsv"
    --trees "${TREE_STAGE}/trees"
    --outdir "${OUT}/prepared"
    --force
)
if [[ -n "${PUG_SPECIES_TREE}" ]]; then
    PREPARE_ARGS+=(--species-tree "${PUG_SPECIES_TREE}")
fi
if [[ -n "${PUG_EXPECTED_LABELS}" ]]; then
    IFS=',' read -r -a EXPECTED_LABEL_ARRAY <<< "${PUG_EXPECTED_LABELS}"
    for label in "${EXPECTED_LABEL_ARRAY[@]}"; do
        PREPARE_ARGS+=(--expected-label "${label}")
    done
fi
"${PYTHON_BIN}" "${ROOT_DIR}/scripts/fig2_prepare_pug_inputs.py" "${PREPARE_ARGS[@]}"

for event in shared_candidate pzi_ab psp_ab ppa_ab; do
    mkdir -p "${OUT}/${event}"
    (
        cd "${OUT}/${event}"
        "${PERL_BIN}" "${PUG_DIR}/PUG.pl" \
            --paralogs "${OUT}/prepared/paralogs_${event}.tsv" \
            --trees "${OUT}/prepared/trees" \
            --outgroups rch \
            --species_tree "${OUT}/prepared/species_tree_for_pug.tre" \
            --name "fig2b_${event}" \
            --tree_type Bayesian \
            > "${OUT}/${event}.log" 2>&1
    )
done

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/fig2_summarize_pug.py" \
    --pug-root "${OUT}" \
    --out "${OUT}/pug_node_counts.tsv" \
    --cutoff 95

{
    echo "families=$(awk 'NR>1 {n++} END {print n+0}' "${INPUT}/selected_collinear_wgd_families.tsv")"
    echo "tree_stage=${TREE_STAGE}"
    echo "pug_dir=${PUG_DIR}"
    echo "support_cutoff=95"
    echo "support_note=FastTree_default local support in 0-1; PUG Bayesian mode used only for x100 scaling, not posterior interpretation"
    echo "outgroup=rch"
    echo "counting=default PUG unique LCA per gene tree; no --all_pairs"
    [[ -n "${PUG_SPECIES_TREE}" ]] && echo "species_tree=${PUG_SPECIES_TREE}"
    [[ -n "${PUG_EXPECTED_LABELS}" ]] && echo "expected_family_tree_labels=${PUG_EXPECTED_LABELS}"
} > "${OUT}/software_and_parameters.txt"
touch "${OUT}/COMPLETE"
echo "[OK] Fig. 2b PUG mapping completed: ${OUT}"
