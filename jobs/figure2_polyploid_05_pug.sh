#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

# Archived HPC equivalent of jobs/figure2_polyploid_05_pug_local.sh.
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/downstream_analysis/figure2_polyploid/strict_wgd_mapping}"
INPUT="${INPUT:-${BASE}/00_family_selection_alltaxa_syntenic_strict_v2}"
TREE_STAGE="${TREE_STAGE:-${BASE}/02_family_trees_reference_default_alltaxa_strict}"
OUT="${OUT:-${BASE}/03_pug_reference_default_alltaxa_strict}"
PUG_DIR="${PUG_DIR:-${ROOT_DIR}/software/PUG}"
PYTHON_BIN="${PYTHON_BIN:-${CONDA_ROOT}/envs/phylo/bin/python}"
PERL_BIN="${PERL_BIN:-/usr/bin/perl}"

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
"${PYTHON_BIN}" "${ROOT_DIR}/scripts/fig2_prepare_pug_inputs.py" \
    --families "${INPUT}/selected_collinear_wgd_families.tsv" \
    --trees "${TREE_STAGE}/trees" \
    --outdir "${OUT}/prepared" \
    --force

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
    --pug-root "${OUT}" --out "${OUT}/pug_node_counts.tsv" --cutoff 95

{
    echo "tree_stage=${TREE_STAGE}"
    echo "support_cutoff=95"
    echo "support_note=FastTree_default local support in 0-1; PUG Bayesian mode used only for x100 scaling, not posterior interpretation"
    echo "outgroup=rch"
    echo "counting=default PUG unique LCA per gene tree; no --all_pairs"
} > "${OUT}/software_and_parameters.txt"
touch "${OUT}/COMPLETE"
echo "[OK] Fig. 2b PUG mapping completed: ${OUT}"
