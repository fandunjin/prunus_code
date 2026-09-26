#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Strict-reference MUSCLE/PAL2NAL/PAML CodeML Ks calculation for one A-B dataset.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CONDA_BASE="${CONDA_BASE:-${HOME}/miniconda3}"
WGD_ENV="${WGD_ENV:-wgd2-reference-clean}"
TOOLS_ENV="${TOOLS_ENV:-wgd2-reference}"
THREADS="${THREADS:-12}"
DATASET="${DATASET:-}"

die() { echo "[ERROR] $*" >&2; exit 2; }

case "${DATASET}" in
  PZI|PSP|PPA) ;;
  *) die "Set DATASET to exactly one of: PZI, PSP, PPA" ;;
esac

WGD_BIN="${CONDA_BASE}/envs/${WGD_ENV}/bin"
TOOLS_BIN="${CONDA_BASE}/envs/${TOOLS_ENV}/bin"
export PATH="${WGD_BIN}:${TOOLS_BIN}:/usr/bin:/bin"
export MPLBACKEND=Agg
export QT_QPA_PLATFORM=offscreen
export LC_ALL=C

for exe in python3 wgd muscle pal2nal.pl codeml; do
  command -v "${exe}" >/dev/null 2>&1 || die "Required executable unavailable: ${exe}"
done

BASE="${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/06_panel_a_reference_exact"
FAMILIES="${BASE}/04_cds_preflight_upper/${DATASET}/${DATASET}_AB.families.cds_valid.tsv"
CDS="${BASE}/01_prepared_inputs/upper/${DATASET}_AB.cds.fa"
OUT="${BASE}/05_wgd_ks_upper/${DATASET}"
WGD_OUT="${OUT}/wgd_ksd"
AUDIT_OUT="${OUT}/audit"
WORK_DIR="${OUT}/work"
SUMMARIZER="${ROOT_DIR}/scripts/fig2_summarize_wgd_pair_ks.py"
KS_TSV="${WGD_OUT}/$(basename "${FAMILIES}").ks.tsv"

[[ -s "${FAMILIES}" ]] || die "Missing CDS-valid families: ${FAMILIES}"
[[ -s "${CDS}" ]] || die "Missing CDS FASTA: ${CDS}"
[[ -s "${SUMMARIZER}" ]] || die "Missing Ks summarizer: ${SUMMARIZER}"
if [[ -e "${OUT}/COMPLETE" ]]; then
  die "WGD Ks calculation is already complete: ${OUT}"
fi
mkdir -p "${WGD_OUT}" "${AUDIT_OUT}" "${WORK_DIR}"

RUN_TMP="/tmp/prunus_fig2a_wgd/${DATASET}_$$"
mkdir -p "${RUN_TMP}"
RESOLVED_TMP="$(realpath -m "${RUN_TMP}")"
case "${RESOLVED_TMP}" in
  /tmp/prunus_fig2a_wgd/${DATASET}_*) ;;
  *) die "Unsafe temporary path: ${RESOLVED_TMP}" ;;
esac

START_EPOCH="$(date +%s)"
if [[ ! -s "${KS_TSV}" ]]; then
  echo "[INFO] Running ${DATASET} pairwise MUSCLE/PAL2NAL/CodeML on $(($(wc -l < "${FAMILIES}") - 1)) families."
  cd "${WORK_DIR}"
  wgd ksd \
    "${FAMILIES}" \
    "${CDS}" \
    --outdir "${WGD_OUT}" \
    --tmpdir "${RUN_TMP}" \
    --nthreads "${THREADS}" \
    --cds \
    --pairwise \
    --aligner muscle \
    --tree_method cluster \
    2>&1 | tee "${OUT}/wgd_ksd.log"
else
  echo "[INFO] Reusing existing WGD Ks table: ${KS_TSV}"
fi
[[ -s "${KS_TSV}" ]] || die "Missing WGD Ks output: ${KS_TSV}"

"${WGD_BIN}/python3" "${SUMMARIZER}" \
  --ks "${KS_TSV}" \
  --families "${FAMILIES}" \
  --dataset "${DATASET}_AB" \
  --max-ks 5 \
  --outdir "${AUDIT_OUT}"

END_EPOCH="$(date +%s)"
WGD_COMMIT="$(git -C "${ROOT_DIR}/software/wgd_v2" rev-parse HEAD 2>/dev/null || echo d2c06a543bc272b2cb5dc469f600380e7dc0771a)"
{
  echo -e "parameter\tvalue"
  echo -e "dataset\t${DATASET}_AB"
  echo -e "wgd_version\t2.0.38"
  echo -e "wgd_commit\t${WGD_COMMIT}"
  echo -e "alignment\tMUSCLE_3.8.1551_protein_then_PAL2NAL"
  echo -e "codeml_runmode\t-2"
  echo -e "codeml_CodonFreq\t2_F3X4"
  echo -e "codeml_cleandata\t1"
  echo -e "pairwise\ttrue"
  echo -e "tree_method\tcluster_two_sequence_family"
  echo -e "Ks_filter\t0_lt_Ks_lt_5"
  echo -e "threads\t${THREADS}"
  echo -e "elapsed_seconds\t$((END_EPOCH - START_EPOCH))"
  echo -e "families\t${FAMILIES}"
  echo -e "cds\t${CDS}"
} > "${OUT}/run_parameters.tsv"

rm -rf -- "${RESOLVED_TMP}"
printf 'Strict-reference WGD pairwise Ks calculation completed.\n' > "${OUT}/COMPLETE"
echo "[OK] Figure 2a ${DATASET} A-B WGD Ks completed: ${OUT}"
cat "${AUDIT_OUT}/${DATASET}_AB.ks_summary.tsv"
