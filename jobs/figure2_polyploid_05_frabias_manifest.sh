#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

# Archive HPC equivalent of the local Fig. 2c manifest stage.
# It deliberately does not launch LAST/Quota/FraBias until those executables
# have been verified in a dedicated environment.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/downstream_analysis/ks}"
PAIR_TABLE="${PAIR_TABLE:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/alltaxa_fractionation_states/summary/gtr_alltaxa_syntenic_pairs.tsv}"
OUT="${OUT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/03_frabias_manifest}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
PY_ENV="${PY_ENV:-akr}"
SCRIPT="${ROOT_DIR}/scripts/fig2_prepare_frabias_manifest.py"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }

require_file "${SCRIPT}"
require_file "${PAIR_TABLE}"
[[ -d "${KS_ROOT}/jcvi/inputs" ]] || die "Missing JCVI inputs: ${KS_ROOT}/jcvi/inputs"
[[ ! -e "${OUT}" ]] || die "Output already exists: ${OUT}"
mkdir -p "${ROOT_DIR}/downstream_analysis/figure2_polyploid/logs"

set +u
source "${CONDA_SH}"
conda activate "${PY_ENV}"
set -u
command -v python >/dev/null 2>&1 || die "python unavailable in ${PY_ENV}"

python "${SCRIPT}" \
    --jcvi-inputs "${KS_ROOT}/jcvi/inputs" \
    --pairs "${PAIR_TABLE}" \
    --outdir "${OUT}"

for file in \
    eight_track_manifest.tsv \
    selected_pzi_homoeologous_pairs.tsv \
    selected_pzi_chromosome_to_all_gtr_chromosomes.tsv \
    method_note.txt; do
    require_file "${OUT}/${file}"
done
{
    echo "python=$(python --version 2>&1)"
    echo "analysis_scope=primary_A1_B1_only"
    echo "tracks=PZI_A,PZI_B,PSP_A,PSP_B,PPA_A,PPA_B,PCAM,PPE"
    echo "selected_pairs=C06_A1-C06_B1,C02_A1-C02_B1,C03_A1-C03_B1"
    echo "next_step=LAST_then_QuotaAlignMerge_20genes_then_QuotaAlign_8to1_then_FraBias_100genes"
    echo "quota_adaptation=eight_merged_query_tracks_vs_GTR"
} > "${OUT}/software_and_parameters.txt"
touch "${OUT}/COMPLETE"
echo "[OK] Fig. 2c FraBias manifest completed: ${OUT}"
