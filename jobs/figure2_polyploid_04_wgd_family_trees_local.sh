#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Local WSL runner for Figure 2b family-tree reconstruction.
# Usage from the repository root:
#   bash jobs/figure2_polyploid_04_wgd_family_trees_local.sh

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping}"
INPUT="${INPUT:-${BASE}/00_family_selection_alltaxa_syntenic_strict_v2}"
OUT="${OUT:-${BASE}/02_family_trees_reference_default_alltaxa_strict}"
THREADS="${THREADS:-12}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }

require_file "${INPUT}/selected_collinear_wgd_families.tsv"
[[ -d "${INPUT}/family_fastas" ]] || die "Missing family FASTAs: ${INPUT}/family_fastas"
[[ ! -e "${OUT}/COMPLETE" ]] || die "Output is already complete: ${OUT}"
command -v mafft >/dev/null 2>&1 || die "mafft is unavailable on PATH"
command -v xargs >/dev/null 2>&1 || die "xargs is unavailable on PATH"

FASTTREE_BIN=""
for candidate in FastTree FastTreeMP fasttree; do
    if command -v "${candidate}" >/dev/null 2>&1; then
        FASTTREE_BIN="$(command -v "${candidate}")"
        break
    fi
done
[[ -n "${FASTTREE_BIN}" ]] || die "FastTree is unavailable on PATH"

mkdir -p "${OUT}/mafft" "${OUT}/trees" "${OUT}/logs"
find "${INPUT}/family_fastas" -maxdepth 1 -type f -name 'OG*.fa' -print0 | sort -z > "${OUT}/family_fastas.list0"
EXPECTED="$(tr -cd '\0' < "${OUT}/family_fastas.list0" | wc -c)"
[[ "${EXPECTED}" -gt 0 ]] || die "No selected family FASTAs found"

export FASTTREE_BIN OUT

run_one_family() {
    input="$1"
    base="$(basename "${input}" .fa)"
    # Reference Methods: MAFFT with its default algorithm/parameters.
    mafft "${input}" \
        > "${OUT}/mafft/${base}.fa" \
        2> "${OUT}/logs/${base}.mafft.log"
    # Reference Methods: FastTree with default values.
    "${FASTTREE_BIN}" "${OUT}/mafft/${base}.fa" \
        > "${OUT}/trees/${base}.tre" \
        2> "${OUT}/logs/${base}.fasttree.log"
}
export -f run_one_family

if command -v parallel >/dev/null 2>&1; then
    cat "${OUT}/family_fastas.list0" | parallel -0 -j "${THREADS}" --halt now,fail=1 \
        'run_one_family "{}"'
else
    # xargs is the standard WSL fallback when GNU parallel is not installed.
    cat "${OUT}/family_fastas.list0" | xargs -0 -r -n 1 -P "${THREADS}" \
        bash -c 'run_one_family "$1"' _
fi

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
    echo "analysis_scope=primary_A1_B1_only"
    echo "method_boundary=tree_reconstruction_only; PUG_or_explicit_reconciliation_required_before_branch_counts"
} > "${OUT}/software_and_parameters.txt"

touch "${OUT}/COMPLETE"
echo "[OK] Fig. 2b family trees completed: ${OUT}"
echo "[NEXT] Apply the reference >=95 support cutoff only after PUG/reconciliation maps duplication events."
