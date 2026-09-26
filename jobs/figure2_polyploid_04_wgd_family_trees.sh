#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

# Figure 2b: reference Fig. 3b-style family-tree reconstruction.
# The preceding selection stage has already restricted families to collinear
# A/B pairs in PZI, PSP and PPA plus GTR and RCH outgroups.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/downstream_analysis/figure2_polyploid/strict_wgd_mapping}"
INPUT="${INPUT:-${BASE}/00_family_selection_alltaxa_syntenic_strict_v2}"
OUT="${OUT:-${BASE}/02_family_trees_reference_default_alltaxa_strict}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
TREE_ENV="${TREE_ENV:-phylo}"
THREADS="${THREADS:-48}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }

require_file "${INPUT}/selected_collinear_wgd_families.tsv"
[[ -d "${INPUT}/family_fastas" ]] || die "Missing family FASTAs: ${INPUT}/family_fastas"
[[ ! -e "${OUT}" ]] || die "Output already exists: ${OUT}"
mkdir -p "${OUT}/mafft" "${OUT}/trees" "${OUT}/logs"

set +u
source "${CONDA_SH}"
conda activate "${TREE_ENV}"
set -u

for executable in mafft parallel; do
    command -v "${executable}" >/dev/null 2>&1 || die "${executable} unavailable in ${TREE_ENV}"
done
FASTTREE_BIN=""
for candidate in FastTree FastTreeMP fasttree; do
    if command -v "${candidate}" >/dev/null 2>&1; then
        FASTTREE_BIN="$(command -v "${candidate}")"
        break
    fi
done
[[ -n "${FASTTREE_BIN}" ]] || die "FastTree unavailable in ${TREE_ENV}; install FastTree before this Fig. 3b replication stage."

find "${INPUT}/family_fastas" -maxdepth 1 -type f -name 'OG*.fa' -print0 | sort -z > "${OUT}/family_fastas.list0"
EXPECTED="$(tr -cd '\0' < "${OUT}/family_fastas.list0" | wc -c)"
[[ "${EXPECTED}" -gt 0 ]] || die "No selected family FASTAs found"

export FASTTREE_BIN OUT
cat "${OUT}/family_fastas.list0" | parallel -0 -j "${THREADS}" --halt now,fail=1 \
    'base=$(basename {} .fa); mafft {} > "${OUT}/mafft/${base}.fa" 2> "${OUT}/logs/${base}.mafft.log"; "${FASTTREE_BIN}" "${OUT}/mafft/${base}.fa" > "${OUT}/trees/${base}.tre" 2> "${OUT}/logs/${base}.fasttree.log"'

OBSERVED="$(find "${OUT}/trees" -maxdepth 1 -type f -name 'OG*.tre' -size +0c | wc -l)"
[[ "${OBSERVED}" -eq "${EXPECTED}" ]] || die "Expected ${EXPECTED} family trees, observed ${OBSERVED}"
find "${OUT}/trees" -maxdepth 1 -type f -name 'OG*.tre' -print0 | sort -z | xargs -0 -r cat > "${OUT}/all_family_trees.tre"
require_file "${OUT}/all_family_trees.tre"

{
    echo "mafft=$(mafft --version 2>&1 | head -n 1)"
    echo "fasttree=$(${FASTTREE_BIN} 2>&1 | head -n 1 || true)"
    echo "families_selected=${EXPECTED}"
    echo "alignment=MAFFT_default"
    echo "tree_inference=FastTree_default"
    echo "support_filter_for_mapping=deferred_to_PUG_reconciliation; reference_caption_requires_mapped_events_ge_95"
    echo "method_boundary=tree_reconstruction_only; PUG_or_explicit_reconciliation_required_before_branch_counts"
} > "${OUT}/software_and_parameters.txt"

touch "${OUT}/COMPLETE"
echo "[OK] Fig. 2b family trees completed: ${OUT}"
echo "[NEXT] Apply the reference >=95 support cutoff only after PUG/reconciliation maps duplication events."
