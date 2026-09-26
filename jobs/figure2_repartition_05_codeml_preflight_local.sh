#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Read-only preflight for the formal, old-panel-c-compatible Figure 2 Ks run.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/downstream_analysis/figure2_repartition_ks}"
CONDA_BASE="${CONDA_BASE:-${HOME}/miniconda3}"
WGD_ENV="${WGD_ENV:-wgd2-reference-clean}"
TOOLS_ENV="${TOOLS_ENV:-wgd2-reference}"
OUT="${OUT:-${BASE}/00_codeml_preflight}"
FORCE="${FORCE:-0}"

die() { echo "[ERROR] $*" >&2; exit 2; }
need() { [[ -s "$1" ]] || die "Missing input: $1"; }
require_exec() { command -v "$1" >/dev/null 2>&1 || die "Required executable unavailable: $1"; }

WGD_BIN="${CONDA_BASE}/envs/${WGD_ENV}/bin"
TOOLS_BIN="${CONDA_BASE}/envs/${TOOLS_ENV}/bin"
PYTHON_BIN="${WGD_BIN}/python3"
[[ -x "${PYTHON_BIN}" ]] || die "Missing WGD Python: ${PYTHON_BIN}"
export PATH="${WGD_BIN}:${TOOLS_BIN}:${HOME}/.local/bin:/usr/bin:/bin"
export MPLBACKEND=Agg
export QT_QPA_PLATFORM=offscreen
export LC_ALL=C

need "${BASE}/COMPLETE"
need "${BASE}/01_jcvi_pairs/COMPLETE"
need "${BASE}/02_iadhore_self/COMPLETE"
for file in \
  "${BASE}/plans/ab_pair_plan.tsv" \
  "${BASE}/plans/external_pair_plan.tsv" \
  "${BASE}/plans/self_pair_plan.tsv" \
  "${ROOT_DIR}/software/wgd_v2/wgd/codeml.py"; do
  need "${file}"
done
for exe in wgd muscle pal2nal.pl codeml; do require_exec "${exe}"; done

[[ "$(awk 'NR > 1 {n++} END {print n + 0}' "${BASE}/plans/ab_pair_plan.tsv")" -eq 5 ]] \
  || die "Expected 5 A-B comparisons"
[[ "$(awk 'NR > 1 {n++} END {print n + 0}' "${BASE}/plans/external_pair_plan.tsv")" -eq 30 ]] \
  || die "Expected 30 external comparisons"
[[ "$(awk 'NR > 1 {n++} END {print n + 0}' "${BASE}/plans/self_pair_plan.tsv")" -eq 13 ]] \
  || die "Expected 13 self-synteny comparisons"
grep -q "'runmode': -2" "${ROOT_DIR}/software/wgd_v2/wgd/codeml.py" \
  || die "WGD source does not expose CodeML runmode=-2"
grep -q "'CodonFreq': 2" "${ROOT_DIR}/software/wgd_v2/wgd/codeml.py" \
  || die "WGD source does not expose CodeML CodonFreq=2"
grep -q "'cleandata': 1" "${ROOT_DIR}/software/wgd_v2/wgd/codeml.py" \
  || die "WGD source does not expose CodeML cleandata=1"

if [[ -e "${OUT}" ]] && [[ -n "$(find "${OUT}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  if [[ "${FORCE}" != "1" ]]; then
    die "Existing preflight output: ${OUT}. Inspect it, or set FORCE=1 to archive and rebuild it."
  fi
  mv "${OUT}" "${OUT}.replaced.$(date +%Y%m%d_%H%M%S)"
fi
mkdir -p "${OUT}"
{
  printf 'parameter\tvalue\n'
  printf 'formal_method\tMUSCLE_3.8.1551_then_PAL2NAL_then_PAML_CodeML_pairwise\n'
  printf 'codeml_runmode\t-2\n'
  printf 'codeml_CodonFreq\t2_F3X4\n'
  printf 'codeml_cleandata\t1\n'
  printf 'ab_comparisons\t5\nexternal_comparisons\t30\nself_comparisons\t13\ntotal_comparisons\t48\n'
  printf 'wgd\t%s\n' "$(command -v wgd)"
  printf 'muscle\t%s\n' "$(command -v muscle)"
  printf 'pal2nal\t%s\n' "$(command -v pal2nal.pl)"
  printf 'codeml\t%s\n' "$(command -v codeml)"
  printf 'wgd_python\t%s\n' "${PYTHON_BIN}"
} > "${OUT}/tool_preflight.tsv"
muscle -version > "${OUT}/muscle_version.txt" 2>&1
wgd ksd --help > "${OUT}/wgd_ksd_help.txt" 2>&1
printf 'Formal Figure 2 CodeML preflight passed.\n' > "${OUT}/COMPLETE"
echo "[OK] Formal CodeML preflight: ${OUT}"
