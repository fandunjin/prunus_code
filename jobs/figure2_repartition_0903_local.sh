#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Re-run the formal five-panel Ks workflow with the corrected 0903 inputs.
# Analysis code and software remain under the project root; only the input
# repartition and all generated results are redirected into 0903.

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
RUN_ROOT="${RUN_ROOT:-${PROJECT_ROOT}/0903}"
REPARTITION_ROOT="${REPARTITION_ROOT:-${RUN_ROOT}/haploid_repartition}"
KS_BASE="${KS_BASE:-${RUN_ROOT}/downstream_analysis/figure2_repartition_ks}"
CANONICAL_ROOT="${CANONICAL_ROOT:-${PROJECT_ROOT}/downstream_analysis/ks/stage/canonical}"
GTR_ROOT="${GTR_ROOT:-${PROJECT_ROOT}/data/gtr}"
THREADS="${THREADS:-12}"

die() { echo "[ERROR] $*" >&2; exit 2; }
need() { [[ -s "$1" ]] || die "Missing input: $1"; }

need "${REPARTITION_ROOT}/inputs/ANNOTATIONS_COMPLETE"
need "${CANONICAL_ROOT}/pcam/primary.gff3"
need "${CANONICAL_ROOT}/ppe/primary.gff3"
need "${GTR_ROOT}/Gtr.NCBI.gene.gff3"
mkdir -p "${RUN_ROOT}/downstream_analysis"

cd "${PROJECT_ROOT}"

echo "[1/7] Preparing normalized 0903 Ks inputs"
ROOT_DIR="${PROJECT_ROOT}" OUT_DIR="${KS_BASE}" \
  REPARTITION_ROOT="${REPARTITION_ROOT}" \
  CANONICAL_ROOT="${CANONICAL_ROOT}" GTR_ROOT="${GTR_ROOT}" \
  FORCE=1 bash jobs/figure2_repartition_02_prepare_ks_inputs_local.sh

echo "[2/7] Building JCVI external and A-B pairs"
ROOT_DIR="${PROJECT_ROOT}" BASE="${KS_BASE}" THREADS="${THREADS}" FORCE=1 \
  bash jobs/figure2_repartition_03_jcvi_pairs_local.sh

echo "[3/7] Building i-ADHoRe self-synteny pairs"
ROOT_DIR="${PROJECT_ROOT}" BASE="${KS_BASE}" THREADS="${THREADS}" FORCE=1 \
  bash jobs/figure2_repartition_04_iadhore_self_pairs_local.sh

echo "[4/7] Checking the formal CodeML environment"
ROOT_DIR="${PROJECT_ROOT}" BASE="${KS_BASE}" FORCE=1 \
  bash jobs/figure2_repartition_05_codeml_preflight_local.sh

echo "[5/7] Preparing pairwise CodeML inputs"
ROOT_DIR="${PROJECT_ROOT}" BASE="${KS_BASE}" FORCE=1 \
  bash jobs/figure2_repartition_06_codeml_prepare_local.sh

echo "[6/7] Estimating Ks for all planned comparisons"
ROOT_DIR="${PROJECT_ROOT}" BASE="${KS_BASE}" GROUP=all THREADS="${THREADS}" \
  bash jobs/figure2_repartition_07_codeml_run_local.sh

echo "[7/7] Combining and plotting the five Ks panels"
ROOT_DIR="${PROJECT_ROOT}" BASE="${KS_BASE}" FORCE=1 \
  bash jobs/figure2_repartition_08_combine_plot_codeml_local.sh

echo "[OK] 0903 formal Ks workflow completed: ${KS_BASE}"
