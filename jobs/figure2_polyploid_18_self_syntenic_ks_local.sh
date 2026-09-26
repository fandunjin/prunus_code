#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Genome-native or independently run subgenome self-syntenic Ks.
#
# Run one dataset at a time so the local machine does not contend for memory
# during all-against-all BLASTP, MCL and CodeML.  PZI/PSP/PPA are the prior
# representative A+B composite analyses; PZI_A/B, PSP_A/B and PPA_A/B are
# independent subgenome analyses and must not be confused with post hoc
# subsets of the composite results.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CONDA_BASE="${CONDA_BASE:-${HOME}/miniconda3}"
WGD_ENV="${WGD_ENV:-wgd2-reference-clean}"
TOOLS_ENV="${TOOLS_ENV:-wgd2-reference}"
DATASET="${DATASET:-}"
THREADS="${THREADS:-12}"
EVALUE="${EVALUE:-1e-10}"
INFLATION="${INFLATION:-2.0}"
IADHORE_BIN="${IADHORE_BIN:-${ROOT_DIR}/software/i-ADHoRe/install/bin/i-adhore}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }
require_exec() { command -v "$1" >/dev/null 2>&1 || die "Required executable unavailable: $1"; }

case "${DATASET}" in
  PZI|PSP|PPA|PZI_A|PZI_B|PSP_A|PSP_B|PPA_A|PPA_B|PPE|PCAM|GTR) ;;
  *) die "Set DATASET to PZI, PSP, PPA, PZI_A, PZI_B, PSP_A, PSP_B, PPA_A, PPA_B, PPE, PCAM or GTR" ;;
esac

WGD_BIN="${CONDA_BASE}/envs/${WGD_ENV}/bin"
TOOLS_BIN="${CONDA_BASE}/envs/${TOOLS_ENV}/bin"
PYTHON_BIN="${WGD_BIN}/python3"
[[ -x "${PYTHON_BIN}" ]] || die "Missing WGD Python: ${PYTHON_BIN}"
[[ -x "${IADHORE_BIN}" ]] || die "Required executable unavailable: ${IADHORE_BIN}"
export PATH="${WGD_BIN}:${TOOLS_BIN}:${HOME}/.local/bin:/usr/bin:/bin"
export MPLBACKEND=Agg
export QT_QPA_PLATFORM=offscreen
export LC_ALL=C

for exe in makeblastdb blastp mcl wgd muscle pal2nal.pl codeml gzip awk; do require_exec "${exe}"; done

BASE="${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/09_within_genome_ks"
TRACKS="${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/04_fractbias_8track_reference_exact/01_catalog/tracks"
PREP="${BASE}/01_prepared/${DATASET}"
MCL_ROOT="${BASE}/02_blastp_mcl/${DATASET}"
IAD_ROOT="${BASE}/03_iadhore/${DATASET}"
PAIR_ROOT="${IAD_ROOT}/02_pairs"
VALID_ROOT="${BASE}/04_cds_validation/${DATASET}"
WGD_ROOT="${BASE}/05_wgd_ks/${DATASET}"

PREP_SCRIPT="${ROOT_DIR}/scripts/fig2_prepare_self_ks_inputs.py"
IAD_INPUT_SCRIPT="${ROOT_DIR}/scripts/fig2_build_self_iadhore_input.py"
EXTRACT_SCRIPT="${ROOT_DIR}/scripts/fig2_extract_self_iadhore_pairs.py"
VALIDATE_SCRIPT="${ROOT_DIR}/scripts/fig2_validate_wgd_pair_families.py"
SUMMARIZE_SCRIPT="${ROOT_DIR}/scripts/fig2_summarize_wgd_pair_ks.py"
for file in "${PREP_SCRIPT}" "${IAD_INPUT_SCRIPT}" "${EXTRACT_SCRIPT}" "${VALIDATE_SCRIPT}" "${SUMMARIZE_SCRIPT}"; do require_file "${file}"; done

COMPONENTS=()
case "${DATASET}" in
  PZI) COMPONENTS=(PZI_A PZI_B) ;;
  PSP) COMPONENTS=(PSP_A PSP_B) ;;
  PPA) COMPONENTS=(PPA_A PPA_B) ;;
  PZI_A|PZI_B|PSP_A|PSP_B|PPA_A|PPA_B) COMPONENTS=("${DATASET}") ;;
  PPE|PCAM|GTR) COMPONENTS=("${DATASET}") ;;
esac
case "${DATASET}" in
  PZI_A|PZI_B|PSP_A|PSP_B|PPA_A|PPA_B) ANALYSIS_SCOPE="independent_subgenome_self_synteny" ;;
  *) ANALYSIS_SCOPE="genome_native_or_representative_composite_self_synteny" ;;
esac
for component in "${COMPONENTS[@]}"; do
  for ext in bed pep cds; do require_file "${TRACKS}/${component}.${ext}"; done
done

if [[ -e "${WGD_ROOT}/COMPLETE" ]]; then
  echo "[OK] Self-syntenic Ks is already complete: ${WGD_ROOT}"
  cat "${WGD_ROOT}/audit/${DATASET}_SELF.ks_summary.tsv"
  exit 0
fi
mkdir -p "${BASE}"

if [[ ! -e "${PREP}/COMPLETE" ]]; then
  PREP_ARGS=(--label "${DATASET}" --outdir "${PREP}")
  for component in "${COMPONENTS[@]}"; do
    PREP_ARGS+=(--bed "${TRACKS}/${component}.bed" --pep "${TRACKS}/${component}.pep" --cds "${TRACKS}/${component}.cds")
  done
  "${PYTHON_BIN}" "${PREP_SCRIPT}" "${PREP_ARGS[@]}"
else
  echo "[INFO] Reusing prepared ${DATASET} inputs: ${PREP}"
fi
require_file "${PREP}/${DATASET}.bed"
require_file "${PREP}/${DATASET}.pep.fa"
require_file "${PREP}/${DATASET}.cds.fa"

DB_DIR="${MCL_ROOT}/db"
DB_PREFIX="${DB_DIR}/${DATASET}_self"
BLAST_GZ="${MCL_ROOT}/${DATASET}_self.blastp_all_vs_all.tsv.gz"
MCL_OUT="${MCL_ROOT}/${DATASET}_self.mcl.I${INFLATION}.txt"
mkdir -p "${DB_DIR}"
if [[ ! -s "${DB_PREFIX}.pin" ]]; then
  makeblastdb -in "${PREP}/${DATASET}.pep.fa" -dbtype prot -out "${DB_PREFIX}" -title "${DATASET}_self"
fi
if [[ ! -s "${BLAST_GZ}" ]]; then
  echo "[INFO] ${DATASET}: running all-against-all BLASTP (E<=${EVALUE})."
  blastp \
    -query "${PREP}/${DATASET}.pep.fa" \
    -db "${DB_PREFIX}" \
    -evalue "${EVALUE}" \
    -max_target_seqs 1000000 \
    -max_hsps 1 \
    -outfmt '6 qseqid sseqid evalue bitscore' \
    -num_threads "${THREADS}" \
    | gzip -c > "${BLAST_GZ}.partial"
  mv "${BLAST_GZ}.partial" "${BLAST_GZ}"
else
  echo "[INFO] Reusing BLASTP output: ${BLAST_GZ}"
fi
if [[ ! -s "${MCL_OUT}" ]]; then
  echo "[INFO] ${DATASET}: clustering BLASTP graph with MCL I=${INFLATION}."
  gzip -dc "${BLAST_GZ}" \
    | awk 'BEGIN {OFS="\t"} NF >= 4 {print $1, $2, $4}' \
    | mcl - --abc -I "${INFLATION}" -te "${THREADS}" -o "${MCL_OUT}.partial"
  mv "${MCL_OUT}.partial" "${MCL_OUT}"
else
  echo "[INFO] Reusing MCL output: ${MCL_OUT}"
fi
printf 'BLASTP/MCL completed.\n' > "${MCL_ROOT}/COMPLETE"

IAD_INPUT="${IAD_ROOT}/01_input"
if [[ ! -e "${IAD_INPUT}/COMPLETE" ]]; then
  mkdir -p "${IAD_ROOT}"
  "${PYTHON_BIN}" "${IAD_INPUT_SCRIPT}" \
    --label "${DATASET}" \
    --bed "${PREP}/${DATASET}.bed" \
    --mcl "${MCL_OUT}" \
    --outdir "${IAD_INPUT}" \
    --threads "${THREADS}"
else
  echo "[INFO] Reusing i-ADHoRe input: ${IAD_INPUT}"
fi
if [[ ! -s "${IAD_INPUT}/iadhore_output/multiplicons.txt" || ! -s "${IAD_INPUT}/iadhore_output/anchorpoints.txt" ]]; then
  echo "[INFO] ${DATASET}: running i-ADHoRe 3.0 self-synteny."
  "${IADHORE_BIN}" "${IAD_INPUT}/iadhore.ini" \
    > "${IAD_ROOT}/iadhore.stdout.log" 2> "${IAD_ROOT}/iadhore.stderr.log"
else
  echo "[INFO] Reusing i-ADHoRe output: ${IAD_INPUT}/iadhore_output"
fi
require_file "${IAD_INPUT}/iadhore_output/multiplicons.txt"
require_file "${IAD_INPUT}/iadhore_output/anchorpoints.txt"
if [[ ! -e "${PAIR_ROOT}/COMPLETE" ]]; then
  "${PYTHON_BIN}" "${EXTRACT_SCRIPT}" \
    --label "${DATASET}" \
    --iadhore-input "${IAD_INPUT}" \
    --outdir "${PAIR_ROOT}"
else
  echo "[INFO] Reusing self-syntenic pairs: ${PAIR_ROOT}"
fi

FAMILIES="${PAIR_ROOT}/${DATASET}_SELF.families.tsv"
if [[ ! -e "${VALID_ROOT}/COMPLETE" ]]; then
  "${PYTHON_BIN}" "${VALIDATE_SCRIPT}" \
    --families "${FAMILIES}" \
    --cds "${PREP}/${DATASET}.cds.fa" \
    --dataset "${DATASET}_SELF" \
    --outdir "${VALID_ROOT}"
else
  echo "[INFO] Reusing CDS validation: ${VALID_ROOT}"
fi
VALID_FAMILIES="${VALID_ROOT}/${DATASET}_SELF.families.cds_valid.tsv"
require_file "${VALID_FAMILIES}"

WGD_OUT="${WGD_ROOT}/wgd_ksd"
AUDIT_OUT="${WGD_ROOT}/audit"
WORK_DIR="${WGD_ROOT}/work"
KS_TSV="${WGD_OUT}/$(basename "${VALID_FAMILIES}").ks.tsv"
mkdir -p "${WGD_OUT}" "${AUDIT_OUT}" "${WORK_DIR}"
RUN_TMP="/tmp/prunus_fig2_self_ks/${DATASET}_$$"
mkdir -p "${RUN_TMP}"
if [[ ! -s "${KS_TSV}" ]]; then
  echo "[INFO] ${DATASET}: pairwise MUSCLE/PAL2NAL/CodeML for $(($(wc -l < "${VALID_FAMILIES}") - 1)) self-syntenic pairs."
  cd "${WORK_DIR}"
  wgd ksd \
    "${VALID_FAMILIES}" \
    "${PREP}/${DATASET}.cds.fa" \
    --outdir "${WGD_OUT}" \
    --tmpdir "${RUN_TMP}" \
    --nthreads "${THREADS}" \
    --cds \
    --pairwise \
    --aligner muscle \
    --tree_method cluster \
    2>&1 | tee "${WGD_ROOT}/wgd_ksd.log"
else
  echo "[INFO] Reusing WGD Ks table: ${KS_TSV}"
fi
require_file "${KS_TSV}"
"${PYTHON_BIN}" "${SUMMARIZE_SCRIPT}" \
  --ks "${KS_TSV}" \
  --families "${VALID_FAMILIES}" \
  --dataset "${DATASET}_SELF" \
  --max-ks 5 \
  --outdir "${AUDIT_OUT}"

{
  printf 'parameter\tvalue\n'
  printf 'dataset\t%s_SELF\n' "${DATASET}"
  printf 'analysis_scope\t%s\n' "${ANALYSIS_SCOPE}"
  printf 'anchor_source\ti-ADHoRe_3.0_cross_chromosome_real_anchorpoints\n'
  printf 'blast_program\tBLASTP\nblast_evalue\t%s\n' "${EVALUE}"
  printf 'mcl_inflation\t%s\n' "${INFLATION}"
  printf 'iadhore_parameters\tgap_size=30;cluster_gap=35;q_value=0.75;prob_cutoff=0.01;anchor_points=3;gg2;FDR\n'
  printf 'alignment\tMUSCLE_3.8.1551_protein_then_PAL2NAL\n'
  printf 'codeml\tNielsen-Yang_pairwise_F3X4_cleandata_1\n'
  printf 'Ks_filter\t0<Ks<5\nthreads\t%s\n' "${THREADS}"
} > "${WGD_ROOT}/run_parameters.tsv"
printf '%s completed.\n' "${ANALYSIS_SCOPE}" > "${WGD_ROOT}/COMPLETE"
echo "[OK] ${DATASET} self-syntenic Ks completed: ${WGD_ROOT}"
cat "${AUDIT_OUT}/${DATASET}_SELF.ks_summary.tsv"
