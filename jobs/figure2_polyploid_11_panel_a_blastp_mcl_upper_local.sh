#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Strict-reference BLASTP/MCL family screen for one Figure 2a A-B dataset.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CONDA_BASE="${CONDA_BASE:-${HOME}/miniconda3}"
WGD_ENV="${WGD_ENV:-wgd2-reference-clean}"
TOOLS_ENV="${TOOLS_ENV:-wgd2-reference}"
THREADS="${THREADS:-12}"
DATASET="${DATASET:-}"
EVALUE="${EVALUE:-1e-10}"
INFLATION="${INFLATION:-2.0}"

die() { echo "[ERROR] $*" >&2; exit 2; }

case "${DATASET}" in
  PZI|PSP|PPA) ;;
  *) die "Set DATASET to exactly one of: PZI, PSP, PPA" ;;
esac

WGD_BIN="${CONDA_BASE}/envs/${WGD_ENV}/bin"
TOOLS_BIN="${CONDA_BASE}/envs/${TOOLS_ENV}/bin"
export PATH="${WGD_BIN}:${TOOLS_BIN}:/usr/bin:/bin"
export LC_ALL=C

for exe in python3 makeblastdb blastp mcl gzip awk; do
  command -v "${exe}" >/dev/null 2>&1 || die "Required executable unavailable: ${exe}"
done

BASE="${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/06_panel_a_reference_exact"
PREP="${BASE}/01_prepared_inputs/upper"
OUT="${BASE}/03_blastp_mcl/${DATASET}"
PEP="${PREP}/${DATASET}_AB.pep.fa"
PAIRS="${PREP}/${DATASET}_AB.iadhore_pairs_pre_mcl.tsv"
DB_DIR="${OUT}/db"
DB_PREFIX="${DB_DIR}/${DATASET}_AB"
BLAST_GZ="${OUT}/${DATASET}_AB.blastp_all_vs_all.tsv.gz"
MCL_OUT="${OUT}/${DATASET}_AB.mcl.I${INFLATION}.txt"
FILTER_DIR="${OUT}/filtered_pairs"

[[ -s "${PEP}" ]] || die "Missing peptide FASTA: ${PEP}"
[[ -s "${PAIRS}" ]] || die "Missing i-ADHoRe pair table: ${PAIRS}"
if [[ -e "${OUT}/COMPLETE" ]]; then
  die "BLASTP/MCL screen is already complete: ${OUT}"
fi
mkdir -p "${DB_DIR}" "${FILTER_DIR}"

if [[ ! -s "${DB_PREFIX}.pin" || ! -e "${DB_DIR}/MAKEBLASTDB.COMPLETE" ]]; then
  makeblastdb -in "${PEP}" -dbtype prot -out "${DB_PREFIX}" -title "${DATASET}_AB"
  printf 'makeblastdb completed.\n' > "${DB_DIR}/MAKEBLASTDB.COMPLETE"
fi

if [[ ! -s "${BLAST_GZ}" || ! -e "${OUT}/BLASTP.COMPLETE" ]]; then
  TMP_GZ="${BLAST_GZ}.partial"
  rm -f "${TMP_GZ}"
  echo "[INFO] Running ${DATASET} all-against-all BLASTP; this is the long step."
  blastp \
    -query "${PEP}" \
    -db "${DB_PREFIX}" \
    -evalue "${EVALUE}" \
    -max_target_seqs 1000000 \
    -max_hsps 1 \
    -outfmt '6 qseqid sseqid evalue bitscore' \
    -num_threads "${THREADS}" \
    | gzip -c > "${TMP_GZ}"
  mv "${TMP_GZ}" "${BLAST_GZ}"
  printf 'blastp completed.\n' > "${OUT}/BLASTP.COMPLETE"
else
  echo "[INFO] Reusing completed BLASTP output: ${BLAST_GZ}"
fi

if [[ ! -s "${MCL_OUT}" || ! -e "${OUT}/MCL.COMPLETE" ]]; then
  TMP_MCL="${MCL_OUT}.partial"
  rm -f "${TMP_MCL}"
  gzip -dc "${BLAST_GZ}" \
    | awk 'BEGIN {OFS="\t"} NF >= 4 {print $1, $2, $4}' \
    | mcl - --abc -I "${INFLATION}" -te "${THREADS}" -o "${TMP_MCL}"
  mv "${TMP_MCL}" "${MCL_OUT}"
  printf 'mcl completed.\n' > "${OUT}/MCL.COMPLETE"
else
  echo "[INFO] Reusing completed MCL output: ${MCL_OUT}"
fi

python3 "${ROOT_DIR}/scripts/fig2_filter_iadhore_pairs_by_mcl.py" \
  --pairs "${PAIRS}" \
  --mcl "${MCL_OUT}" \
  --dataset "${DATASET}" \
  --outdir "${FILTER_DIR}"

BLAST_ROWS="$(gzip -dc "${BLAST_GZ}" | wc -l)"
MCL_FAMILIES="$(wc -l < "${MCL_OUT}")"
{
  echo -e "parameter\tvalue"
  echo -e "dataset\t${DATASET}"
  echo -e "blast_program\tBLASTP"
  echo -e "blast_evalue\t${EVALUE}"
  echo -e "blast_max_target_seqs\t1000000"
  echo -e "blast_max_hsps\t1"
  echo -e "mcl_inflation\t${INFLATION}"
  echo -e "threads\t${THREADS}"
  echo -e "blast_rows\t${BLAST_ROWS}"
  echo -e "mcl_families\t${MCL_FAMILIES}"
  echo -e "peptide_fasta\t${PEP}"
  echo -e "anchor_pairs\t${PAIRS}"
} > "${OUT}/run_parameters.tsv"

printf 'Strict-reference BLASTP/MCL A-B family screen completed.\n' > "${OUT}/COMPLETE"
echo "[OK] Figure 2a ${DATASET} BLASTP/MCL screen completed: ${OUT}"
cat "${FILTER_DIR}/${DATASET}_AB.mcl_pair_summary.tsv"
