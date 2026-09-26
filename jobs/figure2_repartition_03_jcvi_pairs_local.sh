#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Build corrected, phase-specific JCVI syntenic anchor pairs for Figure 2 Ks.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/downstream_analysis/figure2_repartition_ks}"
TRACKS="${TRACKS:-${BASE}/inputs}"
PLAN="${PLAN:-${BASE}/plans/jcvi_pair_plan.tsv}"
OUT="${OUT:-${BASE}/01_jcvi_pairs}"
WORK="${OUT}/work"
THREADS="${THREADS:-12}"
CSCORE="${CSCORE:-0.70}"
MIN_SIZE="${MIN_SIZE:-4}"
MINSPAN="${MINSPAN:-30}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
FALLBACK_JCVI_PYTHON="${FALLBACK_JCVI_PYTHON:-/opt/plant-wakr-venv/bin/python3}"
FORCE="${FORCE:-0}"

die() { echo "[ERROR] $*" >&2; exit 2; }
need() { [[ -s "$1" ]] || die "Missing input: $1"; }

need "${BASE}/COMPLETE"
need "${PLAN}"
need "${ROOT_DIR}/scripts/ks_jcvi_anchors_to_pairs.py"
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || die "Python unavailable: ${PYTHON_BIN}"
command -v diamond >/dev/null 2>&1 || die "diamond unavailable on PATH"
if ! "${PYTHON_BIN}" -m jcvi.compara.catalog ortholog --help >/dev/null 2>&1; then
  if [[ -x "${FALLBACK_JCVI_PYTHON}" ]] && "${FALLBACK_JCVI_PYTHON}" -m jcvi.compara.catalog ortholog --help >/dev/null 2>&1; then
    echo "[INFO] Switching to the local JCVI Python: ${FALLBACK_JCVI_PYTHON}"
    PYTHON_BIN="${FALLBACK_JCVI_PYTHON}"
  else
    die "JCVI catalog unavailable in ${PYTHON_BIN}; set PYTHON_BIN to the Python environment containing JCVI"
  fi
fi

if [[ -e "${OUT}/COMPLETE" && "${FORCE}" != "1" ]]; then
  echo "[OK] Reusing corrected JCVI pair table: ${OUT}"
  cat "${OUT}/jcvi_pair_status.tsv"
  exit 0
fi
if [[ -e "${OUT}" && "${FORCE}" == "1" ]]; then
  mv "${OUT}" "${OUT}.replaced.$(date +%Y%m%d_%H%M%S)"
fi
mkdir -p "${WORK}"

while IFS=$'\t' read -r comparison query ref comparison_type panel curve; do
  [[ "${comparison}" == "comparison" ]] && continue
  [[ -n "${comparison}" ]] || continue
  for label in "${query}" "${ref}"; do
    for ext in bed pep cds; do
      need "${TRACKS}/${label}/${label}.${ext}"
    done
  done
  pair_work="${WORK}/${comparison}"
  mkdir -p "${pair_work}"
  for label in "${query}" "${ref}"; do
    for ext in bed pep cds; do
      ln -sfn "${TRACKS}/${label}/${label}.${ext}" "${pair_work}/${label}.${ext}"
    done
  done
  raw="${pair_work}/${query}.${ref}.anchors"
  screened="${pair_work}/${query}.${ref}.anchors.new"
  if [[ ! -s "${raw}" ]]; then
    echo "[INFO] JCVI catalog: ${comparison}"
    (
      cd "${pair_work}"
      "${PYTHON_BIN}" -m jcvi.compara.catalog ortholog "${query}" "${ref}" \
        --dbtype=prot --align_soft=diamond_blastp --cscore="${CSCORE}" \
        --min_size="${MIN_SIZE}" --cpus="${THREADS}" --no_strip_names --no_dotplot
    )
  fi
  need "${raw}"
  if [[ ! -s "${screened}" ]]; then
    echo "[INFO] JCVI screen: ${comparison}"
    "${PYTHON_BIN}" -m jcvi.compara.synteny screen \
      --qbed="${pair_work}/${query}.bed" --sbed="${pair_work}/${ref}.bed" \
      --minspan="${MINSPAN}" --simple "${raw}" "${screened}"
  fi
  need "${screened}"
done < "${PLAN}"

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/ks_jcvi_anchors_to_pairs.py" \
  --pair-plan "${PLAN}" --pair-root "${WORK}" \
  --out-pairs "${OUT}/jcvi_syntenic_pairs.tsv" \
  --out-blocks "${OUT}/jcvi_syntenic_blocks.tsv" \
  --out-status "${OUT}/jcvi_pair_status.tsv"

need "${OUT}/jcvi_syntenic_pairs.tsv"
[[ "$(awk 'NR > 1 {n++} END {print n + 0}' "${OUT}/jcvi_pair_status.tsv")" -eq 35 ]] || die "Expected 35 corrected JCVI comparisons"
printf 'Corrected phase-specific JCVI anchor pairs completed.\n' > "${OUT}/COMPLETE"
printf 'parameter\tvalue\ncscore\t%s\nmin_size\t%s\nminspan\t%s\n' "${CSCORE}" "${MIN_SIZE}" "${MINSPAN}" > "${OUT}/run_parameters.tsv"
echo "[OK] Corrected Figure 2 JCVI pairs: ${OUT}"
