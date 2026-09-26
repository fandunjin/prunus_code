#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Prepare separate two-gene CodeML inputs from the corrected repartition pairs.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/downstream_analysis/figure2_repartition_ks}"
OUT="${OUT:-${BASE}/04_codeml_inputs}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
FORCE="${FORCE:-0}"

die() { echo "[ERROR] $*" >&2; exit 2; }
need() { [[ -s "$1" ]] || die "Missing input: $1"; }

need "${BASE}/COMPLETE"
need "${BASE}/01_jcvi_pairs/COMPLETE"
need "${BASE}/02_iadhore_self/COMPLETE"
need "${ROOT_DIR}/scripts/prepare_fig2_repartition_codeml_inputs.py"
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || die "Python unavailable: ${PYTHON_BIN}"

for group in ab external self; do
  group_out="${OUT}/${group}"
  if [[ -e "${group_out}/COMPLETE" && "${FORCE}" != "1" ]]; then
    echo "[INFO] Reusing prepared ${group} CodeML inputs: ${group_out}"
    continue
  fi
  if [[ -e "${group_out}" ]] && [[ -n "$(find "${group_out}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
    if [[ "${FORCE}" != "1" ]]; then
      die "Incomplete non-empty output: ${group_out}. Inspect it before rerunning, then use FORCE=1 to archive and rebuild it."
    fi
    mv "${group_out}" "${group_out}.replaced.$(date +%Y%m%d_%H%M%S)"
  fi
  "${PYTHON_BIN}" "${ROOT_DIR}/scripts/prepare_fig2_repartition_codeml_inputs.py" \
    --base "${BASE}" --group "${group}" --outdir "${group_out}"
done

for group in ab external self; do need "${OUT}/${group}/comparisons.tsv"; done
printf 'All corrected CodeML inputs prepared.\n' > "${OUT}/COMPLETE"
echo "[OK] Figure 2 CodeML inputs: ${OUT}"
