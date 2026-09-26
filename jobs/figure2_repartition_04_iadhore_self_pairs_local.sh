#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Build independent i-ADHoRe cross-chromosome self-synteny pairs for five Ks panels.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE="${BASE:-${ROOT_DIR}/downstream_analysis/figure2_repartition_ks}"
TRACKS="${TRACKS:-${BASE}/inputs}"
PLAN="${PLAN:-${BASE}/plans/self_pair_plan.tsv}"
OUT="${OUT:-${BASE}/02_iadhore_self}"
THREADS="${THREADS:-12}"
EVALUE="${EVALUE:-1e-10}"
INFLATION="${INFLATION:-2.0}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
PROJECT_IADHORE_BIN="${ROOT_DIR}/software/i-ADHoRe/install/bin/i-adhore"
IADHORE_BIN="${IADHORE_BIN:-${PROJECT_IADHORE_BIN}}"
IADHORE_ENV="${IADHORE_ENV:-figure2_iadhore}"
CONDA_BASE="${CONDA_BASE:-${HOME}/miniconda3}"
FORCE="${FORCE:-0}"

die() { echo "[ERROR] $*" >&2; exit 2; }
need() { [[ -s "$1" ]] || die "Missing input: $1"; }

need "${BASE}/COMPLETE"
need "${PLAN}"
for script in fig2_build_self_iadhore_input.py fig2_extract_self_iadhore_pairs.py combine_figure2_repartition_self_pairs.py; do
  need "${ROOT_DIR}/scripts/${script}"
done
IADHORE_ENV_BIN="${CONDA_BASE}/envs/${IADHORE_ENV}/bin"
resolve_tool() {
  local requested="$1"
  local fallback="$2"
  if command -v "${requested}" >/dev/null 2>&1; then
    command -v "${requested}"
  elif [[ -x "${fallback}" ]]; then
    printf '%s\n' "${fallback}"
  else
    return 1
  fi
}
MAKEBLASTDB_BIN="$(resolve_tool makeblastdb "${IADHORE_ENV_BIN}/makeblastdb")" \
  || die "makeblastdb unavailable; run jobs/install_figure2_iadhore_local.sh"
BLASTP_BIN="$(resolve_tool blastp "${IADHORE_ENV_BIN}/blastp")" \
  || die "blastp unavailable; run jobs/install_figure2_iadhore_local.sh"
MCL_BIN="$(resolve_tool mcl "${IADHORE_ENV_BIN}/mcl")" \
  || die "mcl unavailable; run jobs/install_figure2_iadhore_local.sh"
IADHORE_BIN="$(resolve_tool "${IADHORE_BIN}" "${IADHORE_ENV_BIN}/i-adhore")" \
  || die "i-adhore unavailable; run jobs/install_figure2_iadhore_local.sh"
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || die "Python unavailable: ${PYTHON_BIN}"

if [[ -e "${OUT}/COMPLETE" && "${FORCE}" != "1" ]]; then
  echo "[OK] Reusing corrected i-ADHoRe self-synteny pairs: ${OUT}"
  cat "${OUT}/self_pair_status.tsv"
  exit 0
fi
if [[ -e "${OUT}" && "${FORCE}" == "1" ]]; then
  mv "${OUT}" "${OUT}.replaced.$(date +%Y%m%d_%H%M%S)"
fi
mkdir -p "${OUT}"

while IFS=$'\t' read -r label panel curve; do
  [[ "${label}" == "label" ]] && continue
  [[ -n "${label}" ]] || continue
  for ext in bed pep cds; do need "${TRACKS}/${label}/${label}.${ext}"; done
  base="${OUT}/${label}"
  mcl_dir="${base}/01_blastp_mcl"
  iad_input="${base}/02_iadhore_input"
  pair_dir="${base}/03_pairs"
  mkdir -p "${mcl_dir}/db"
  db="${mcl_dir}/db/${label}_self"
  blast="${mcl_dir}/${label}_all_vs_all.tsv.gz"
  mcl_out="${mcl_dir}/${label}.mcl.I${INFLATION}.txt"
  if [[ ! -s "${db}.pin" ]]; then
    "${MAKEBLASTDB_BIN}" -in "${TRACKS}/${label}/${label}.pep" -dbtype prot -out "${db}" -title "${label}_self"
  fi
  if [[ ! -s "${blast}" ]]; then
    echo "[INFO] ${label}: BLASTP all-versus-all"
    "${BLASTP_BIN}" -query "${TRACKS}/${label}/${label}.pep" -db "${db}" -evalue "${EVALUE}" \
      -max_target_seqs 1000000 -max_hsps 1 -outfmt '6 qseqid sseqid evalue bitscore' \
      -num_threads "${THREADS}" | gzip -c > "${blast}.partial"
    mv "${blast}.partial" "${blast}"
  fi
  if [[ ! -s "${mcl_out}" ]]; then
    echo "[INFO] ${label}: MCL inflation=${INFLATION}"
    gzip -dc "${blast}" | awk 'BEGIN {OFS="\t"} NF >= 4 {print $1,$2,$4}' | \
      "${MCL_BIN}" - --abc -I "${INFLATION}" -te "${THREADS}" -o "${mcl_out}.partial"
    mv "${mcl_out}.partial" "${mcl_out}"
  fi
  if [[ ! -e "${iad_input}/COMPLETE" ]]; then
    "${PYTHON_BIN}" "${ROOT_DIR}/scripts/fig2_build_self_iadhore_input.py" \
      --label "${label}" --bed "${TRACKS}/${label}/${label}.bed" \
      --mcl "${mcl_out}" --outdir "${iad_input}" --threads "${THREADS}"
  fi
  if [[ ! -s "${iad_input}/iadhore_output/multiplicons.txt" || ! -s "${iad_input}/iadhore_output/anchorpoints.txt" ]]; then
    echo "[INFO] ${label}: i-ADHoRe self-synteny"
    "${IADHORE_BIN}" "${iad_input}/iadhore.ini" > "${base}/iadhore.stdout.log" 2> "${base}/iadhore.stderr.log"
  fi
  if [[ ! -e "${pair_dir}/COMPLETE" ]]; then
    "${PYTHON_BIN}" "${ROOT_DIR}/scripts/fig2_extract_self_iadhore_pairs.py" \
      --label "${label}" --iadhore-input "${iad_input}" --outdir "${pair_dir}"
  fi
done < "${PLAN}"

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/combine_figure2_repartition_self_pairs.py" \
  --self-plan "${PLAN}" --self-root "${OUT}" \
  --out-pairs "${OUT}/self_syntenic_pairs.tsv" --out-status "${OUT}/self_pair_status.tsv"
need "${OUT}/self_syntenic_pairs.tsv"
[[ "$(awk 'NR > 1 {n++} END {print n + 0}' "${OUT}/self_pair_status.tsv")" -eq 13 ]] || die "Expected 13 corrected self-synteny tracks"
printf 'Corrected independent i-ADHoRe self-synteny pairs completed.\n' > "${OUT}/COMPLETE"
{
  printf 'parameter\tvalue\n'
  printf 'blastp_evalue\t%s\nmcl_inflation\t%s\n' "${EVALUE}" "${INFLATION}"
  printf 'iadhore\tgap_size=30;cluster_gap=35;q_value=0.75;prob_cutoff=0.01;anchor_points=3;gg2;FDR\n'
} > "${OUT}/run_parameters.tsv"
echo "[OK] Corrected Figure 2 i-ADHoRe self pairs: ${OUT}"
