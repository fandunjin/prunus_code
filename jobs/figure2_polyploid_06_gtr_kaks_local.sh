#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Local WSL runner for Figure 2a GTR-comparison Ks profiles.
# Uses the existing all-taxon JCVI anchors and the project-standard YN Ka/Ks
# filters. It does not recompute within-species A-B comparisons.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
JCVI_DIR="${JCVI_DIR:-${ROOT_DIR}/downstream_analysis/ks/jcvi}"
PAIR_TABLE="${PAIR_TABLE:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/alltaxa_fractionation_states/summary/gtr_alltaxa_syntenic_pairs.tsv}"
OUT="${OUT:-${ROOT_DIR}/downstream_analysis/ks/kaks/gtr_alltaxa_collinear_all_pairs_YN}"
THREADS="${THREADS:-12}"
KAKS_BIN="${KAKS_BIN:-}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }

require_file "${ROOT_DIR}/scripts/downstream_wgd_kaks_from_jcvi.py"
require_file "${PAIR_TABLE}"
[[ ! -e "${OUT}" ]] || die "Output directory already exists: ${OUT}"
for exe in python3 mafft pal2nal.pl; do
    command -v "${exe}" >/dev/null 2>&1 || die "Required executable unavailable: ${exe}"
done
if [[ -z "${KAKS_BIN}" ]]; then
    if command -v KaKs_Calculator >/dev/null 2>&1; then
        KAKS_BIN="$(command -v KaKs_Calculator)"
    elif command -v KaKs >/dev/null 2>&1; then
        KAKS_BIN="$(command -v KaKs)"
    elif [[ -x "${ROOT_DIR}/software/KaKs_Calculator3.0_linux/bin/KaKs" ]]; then
        KAKS_BIN="${ROOT_DIR}/software/KaKs_Calculator3.0_linux/bin/KaKs"
    else
        die "KaKs_Calculator/KaKs unavailable; run jobs/install_kaks_calculator_local.sh"
    fi
fi
for label in pzi_hapA1 pzi_hapA2 pzi_hapB1 pzi_hapB2 psp_hapA1 psp_hapA2 psp_hapB1 psp_hapB2 ppa_hapA ppa_hapB pcam ppe gtr; do
    require_file "${JCVI_DIR}/inputs/${label}/${label}.pep"
    require_file "${JCVI_DIR}/inputs/${label}/${label}.cds"
done

python3 "${ROOT_DIR}/scripts/downstream_wgd_kaks_from_jcvi.py" \
    --analysis-dir "${ROOT_DIR}" \
    --jcvi-dir "${JCVI_DIR}" \
    --pair-table "${PAIR_TABLE}" \
    --outdir "${OUT}" \
    --threads "${THREADS}" \
    --pair-mode collinear_all_pairs \
    --include-comparison pzi_hapA1_vs_gtr \
    --include-comparison pzi_hapB1_vs_gtr \
    --include-comparison psp_hapA1_vs_gtr \
    --include-comparison psp_hapB1_vs_gtr \
    --include-comparison ppa_hapA_vs_gtr \
    --include-comparison ppa_hapB_vs_gtr \
    --include-comparison pcam_vs_gtr \
    --include-comparison ppe_vs_gtr \
    --min-aa 50 \
    --min-codons 50 \
    --min-4d-sites 4 \
    --min-reference-block-anchors 20 \
    --kaks-cmd "${KAKS_BIN}" \
    --kaks-method YN

require_file "${OUT}/summary/all_merged_4dtv_kaks.tsv"
require_file "${OUT}/summary/wgd_kaks_comparison_summary.tsv"
{
    echo "pair_table=${PAIR_TABLE}"
    echo "pair_mode=collinear_all_pairs"
    echo "kaks_method=YN"
    echo "kaks_cmd=${KAKS_BIN}"
    echo "min_aa=50"
    echo "min_codons=50"
    echo "min_4d_sites=4"
    echo "min_reference_block_anchors=20"
    echo "analysis_role=Figure2a_external_GTR_context"
} > "${OUT}/software_and_parameters.txt"
touch "${OUT}/COMPLETE"
echo "[OK] Figure 2a all-taxon GTR Ks completed: ${OUT}"
