#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Recalculate Figure 2 Ks with the old panel-c CodeML contract, one comparison
# at a time: MUSCLE 3.8.1551 -> PAL2NAL -> PAML CodeML (F3X4, pairwise).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/downstream_analysis/figure2_repartition_ks}"
INPUT_ROOT="${INPUT_ROOT:-${BASE}/04_codeml_inputs}"
OUT="${OUT:-${BASE}/05_codeml_results}"
CONDA_BASE="${CONDA_BASE:-${HOME}/miniconda3}"
WGD_ENV="${WGD_ENV:-wgd2-reference-clean}"
TOOLS_ENV="${TOOLS_ENV:-wgd2-reference}"
GROUP="${GROUP:-all}"
COMPARISON="${COMPARISON:-}"
THREADS="${THREADS:-12}"

die() { echo "[ERROR] $*" >&2; exit 2; }
need() { [[ -s "$1" ]] || die "Missing input: $1"; }
require_exec() { command -v "$1" >/dev/null 2>&1 || die "Required executable unavailable: $1"; }

case "${GROUP}" in all|ab|external|self) ;; *) die "GROUP must be all, ab, external, or self" ;; esac
[[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] || die "THREADS must be a positive integer"

WGD_BIN="${CONDA_BASE}/envs/${WGD_ENV}/bin"
TOOLS_BIN="${CONDA_BASE}/envs/${TOOLS_ENV}/bin"
PYTHON_BIN="${WGD_BIN}/python3"
[[ -x "${PYTHON_BIN}" ]] || die "Missing WGD Python: ${PYTHON_BIN}"
export PATH="${WGD_BIN}:${TOOLS_BIN}:${HOME}/.local/bin:/usr/bin:/bin"
export MPLBACKEND=Agg
export QT_QPA_PLATFORM=offscreen
export LC_ALL=C

for exe in wgd muscle pal2nal.pl codeml tee awk realpath; do require_exec "${exe}"; done
for script in fig2_validate_wgd_pair_families.py fig2_summarize_wgd_pair_ks.py; do need "${ROOT_DIR}/scripts/${script}"; done
need "${INPUT_ROOT}/COMPLETE"

run_one() {
  local group="$1"
  local comparison="$2"
  local source_dir="${INPUT_ROOT}/${group}/${comparison}"
  local result_dir="${OUT}/${group}/${comparison}"
  local validation_dir="${result_dir}/validation"
  local wgd_dir="${result_dir}/wgd_ksd"
  local audit_dir="${result_dir}/audit"
  local work_dir="${result_dir}/work"
  local families="${source_dir}/families.tsv"
  local cds="${source_dir}/combined.cds.fa"
  local valid_families="${validation_dir}/${comparison}.families.cds_valid.tsv"
  local ks_tsv="${wgd_dir}/$(basename "${valid_families}").ks.tsv"
  local run_tmp="/tmp/prunus_fig2_repartition_codeml/${group}/${comparison}_$$"
  local resolved_tmp
  local start_epoch end_epoch wgd_status

  need "${families}"
  need "${cds}"
  if [[ -e "${result_dir}/COMPLETE" ]]; then
    echo "[INFO] Reusing completed CodeML result: ${group}/${comparison}"
    return
  fi
  mkdir -p "${validation_dir}" "${wgd_dir}" "${audit_dir}" "${work_dir}" "${run_tmp}"

  if [[ ! -e "${validation_dir}/COMPLETE" ]]; then
    "${PYTHON_BIN}" "${ROOT_DIR}/scripts/fig2_validate_wgd_pair_families.py" \
      --families "${families}" --cds "${cds}" --dataset "${comparison}" --outdir "${validation_dir}"
  fi
  need "${valid_families}"

  resolved_tmp="$(realpath -m "${run_tmp}")"
  case "${resolved_tmp}" in
    /tmp/prunus_fig2_repartition_codeml/${group}/${comparison}_*) ;;
    *) die "Unsafe temporary path: ${resolved_tmp}" ;;
  esac
  start_epoch="$(date +%s)"
  if [[ ! -s "${ks_tsv}" ]]; then
    echo "[INFO] ${group}/${comparison}: CodeML for $(($(wc -l < "${valid_families}") - 1)) two-gene families."
    set +e
    (
      cd "${work_dir}"
      wgd ksd \
        "${valid_families}" "${cds}" \
        --outdir "${wgd_dir}" --tmpdir "${resolved_tmp}" --nthreads "${THREADS}" \
        --cds --pairwise --aligner muscle --tree_method cluster
    ) 2>&1 | tee "${result_dir}/wgd_ksd.log"
    wgd_status="${PIPESTATUS[0]}"
    set -e
    # WGD may return non-zero only while writing its optional plotting files.
    # A non-empty pairwise table is still required before the result is accepted.
    if [[ "${wgd_status}" -ne 0 ]]; then
      echo "[WARN] wgd exited ${wgd_status}; checking for a complete Ks table."
    fi
  else
    echo "[INFO] Reusing existing CodeML Ks table: ${ks_tsv}"
  fi
  need "${ks_tsv}"
  "${PYTHON_BIN}" "${ROOT_DIR}/scripts/fig2_summarize_wgd_pair_ks.py" \
    --ks "${ks_tsv}" --families "${valid_families}" --dataset "${comparison}" \
    --max-ks 5 --outdir "${audit_dir}"
  need "${audit_dir}/${comparison}.ks_0_5.tsv"
  end_epoch="$(date +%s)"

  {
    printf 'parameter\tvalue\n'
    printf 'comparison\t%s\n' "${comparison}"
    printf 'input_pairs\tphase_specific_corrected_JCVI_or_iADHoRe\n'
    printf 'alignment\tMUSCLE_3.8.1551_protein_then_PAL2NAL\n'
    printf 'codeml_runmode\t-2\n'
    printf 'codeml_CodonFreq\t2_F3X4\n'
    printf 'codeml_cleandata\t1\n'
    printf 'pairwise\ttrue\n'
    printf 'tree_method\tcluster_two_sequence_family\n'
    printf 'Ks_filter\t0_lt_Ks_lt_5\n'
    printf 'threads\t%s\n' "${THREADS}"
    printf 'elapsed_seconds\t%s\n' "$((end_epoch - start_epoch))"
    printf 'families\t%s\n' "${families}"
    printf 'cds\t%s\n' "${cds}"
  } > "${result_dir}/run_parameters.tsv"
  rm -rf -- "${resolved_tmp}"
  printf 'Pairwise CodeML Ks completed.\n' > "${result_dir}/COMPLETE"
  echo "[OK] CodeML completed: ${group}/${comparison}"
}

groups=()
if [[ "${GROUP}" == "all" ]]; then groups=(ab external self); else groups=("${GROUP}"); fi
for group_name in "${groups[@]}"; do
  manifest="${INPUT_ROOT}/${group_name}/comparisons.tsv"
  need "${manifest}"
  matched=0
  while IFS=$'\t' read -r comparison group_value query_label ref_label comparison_type panel curve families unique_cds source; do
    [[ "${comparison}" == "comparison" ]] && continue
    [[ -n "${comparison}" ]] || continue
    if [[ -n "${COMPARISON}" && "${comparison}" != "${COMPARISON}" ]]; then continue; fi
    matched=1
    run_one "${group_name}" "${comparison}"
  done < "${manifest}"
  if [[ -n "${COMPARISON}" && "${matched}" -eq 0 ]]; then
    die "COMPARISON is not present in ${group_name}: ${COMPARISON}"
  fi
done

all_complete=1
for required_group in ab external self; do
  required_manifest="${INPUT_ROOT}/${required_group}/comparisons.tsv"
  need "${required_manifest}"
  while IFS=$'\t' read -r required_comparison _; do
    [[ "${required_comparison}" == "comparison" ]] && continue
    [[ -n "${required_comparison}" ]] || continue
    if [[ ! -e "${OUT}/${required_group}/${required_comparison}/COMPLETE" ]]; then
      all_complete=0
      break
    fi
  done < "${required_manifest}"
  [[ "${all_complete}" -eq 1 ]] || break
done
if [[ "${all_complete}" -eq 1 ]]; then
  printf 'All planned Figure 2 CodeML comparisons completed.\n' > "${OUT}/COMPLETE"
  echo "[OK] All 48 planned CodeML comparisons are complete."
else
  echo "[INFO] Requested comparisons completed; the formal 48-comparison set is not complete yet."
fi
echo "[OK] Requested CodeML run completed: GROUP=${GROUP} COMPARISON=${COMPARISON:-all}"
