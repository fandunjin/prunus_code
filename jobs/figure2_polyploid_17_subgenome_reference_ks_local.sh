#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

# Build the missing PSP/PPA-versus-PCAM/PPE JCVI pairs, calculate YN Ks, and
# plot all 18 independent A/B subgenome-versus-reference distributions.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
JCVI_DIR="${JCVI_DIR:-${ROOT_DIR}/downstream_analysis/ks/jcvi}"
PLAN="${PLAN:-${ROOT_DIR}/config/figure2_subgenome_reference_missing_plan.tsv}"
BASE_OUT="${BASE_OUT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/subgenome_reference_ks_YN}"
PAIR_ROOT="${BASE_OUT}/00_jcvi_pairs"
PAIR_WORK="${PAIR_ROOT}/work"
NEW_PAIR_TABLE="${PAIR_ROOT}/missing_reference_syntenic_pairs.tsv"
NEW_PAIR_BLOCKS="${PAIR_ROOT}/missing_reference_blocks.tsv"
NEW_PAIR_STATUS="${PAIR_ROOT}/missing_reference_status.tsv"
KAKS_OUT="${BASE_OUT}/01_missing_reference_kaks_YN"
SOURCE_OUT="${BASE_OUT}/02_combined_source_data"
PLOT_OUT="${ROOT_DIR}/downstream_analysis/figure2_polyploid/fig2a_subgenome_reference"
THREADS="${THREADS:-12}"
CSCORE="${CSCORE:-0.70}"
MIN_ANCHORS="${MIN_ANCHORS:-4}"
MINSPAN="${MINSPAN:-30}"
MIN_REFERENCE_BLOCK_ANCHORS="${MIN_REFERENCE_BLOCK_ANCHORS:-20}"
KAKS_BIN="${KAKS_BIN:-${ROOT_DIR}/software/KaKs_Calculator3.0_linux/bin/KaKs}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
FORCE_PLOT="${FORCE_PLOT:-0}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input: $1"; }

for file in \
  "${ROOT_DIR}/scripts/ks_jcvi_anchors_to_pairs.py" \
  "${ROOT_DIR}/scripts/downstream_wgd_kaks_from_jcvi.py" \
  "${ROOT_DIR}/scripts/fig2_prepare_subgenome_reference_ks.py" \
  "${ROOT_DIR}/scripts/plot_fig2_subgenome_reference_ks.py" \
  "${PLAN}" \
  "${ROOT_DIR}/downstream_analysis/ks/kaks/gtr_alltaxa_collinear_all_pairs_YN/summary/all_merged_4dtv_kaks.tsv" \
  "${ROOT_DIR}/downstream_analysis/ks/kaks/collinear_all_pairs_YN/summary/all_merged_4dtv_kaks.tsv"; do
  require_file "${file}"
done

for label in pzi_hapA1 pzi_hapB1 psp_hapA1 psp_hapB1 ppa_hapA ppa_hapB pcam ppe gtr; do
  require_file "${JCVI_DIR}/inputs/${label}/${label}.pep"
  require_file "${JCVI_DIR}/inputs/${label}/${label}.cds"
  require_file "${JCVI_DIR}/inputs/${label}/${label}.bed"
done
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || die "Python unavailable: ${PYTHON_BIN}"
command -v jcvi >/dev/null 2>&1 || die "jcvi unavailable on PATH"
command -v diamond >/dev/null 2>&1 || die "diamond unavailable on PATH"
command -v mafft >/dev/null 2>&1 || die "mafft unavailable on PATH"
command -v pal2nal.pl >/dev/null 2>&1 || die "pal2nal.pl unavailable on PATH"
[[ -x "${KAKS_BIN}" ]] || die "KaKs executable unavailable: ${KAKS_BIN}"

mkdir -p "${PAIR_WORK}" "${KAKS_OUT}" "${SOURCE_OUT}" "${PLOT_OUT}"

if [[ ! -e "${PAIR_ROOT}/COMPLETE" ]]; then
  while IFS=$'\t' read -r comparison query ref comparison_type analysis_role; do
    [[ "${comparison}" == "comparison" ]] && continue
    [[ -n "${comparison}" ]] || continue
    work="${PAIR_WORK}/${comparison}"
    mkdir -p "${work}"
    for label in "${query}" "${ref}"; do
      for ext in bed pep cds; do
        require_file "${JCVI_DIR}/inputs/${label}/${label}.${ext}"
        ln -sfn "${JCVI_DIR}/inputs/${label}/${label}.${ext}" "${work}/${label}.${ext}"
      done
    done
    cd "${work}"
    raw="${query}.${ref}.anchors"
    screened="${query}.${ref}.anchors.new"
    if [[ ! -s "${raw}" ]]; then
      echo "[INFO] JCVI catalog: ${comparison}"
      "${PYTHON_BIN}" -m jcvi.compara.catalog ortholog "${query}" "${ref}" \
        --dbtype=prot --align_soft=diamond_blastp --cscore="${CSCORE}" \
        --min_size="${MIN_ANCHORS}" --cpus="${THREADS}" --no_strip_names --no_dotplot
    else
      echo "[INFO] Reusing JCVI anchors: ${comparison}"
    fi
    require_file "${raw}"
    if [[ ! -s "${screened}" ]]; then
      "${PYTHON_BIN}" -m jcvi.compara.synteny screen \
        --qbed="${query}.bed" --sbed="${ref}.bed" \
        --minspan="${MINSPAN}" --simple "${raw}" "${screened}"
    fi
    require_file "${screened}"
    cd "${ROOT_DIR}"
  done < "${PLAN}"

  "${PYTHON_BIN}" "${ROOT_DIR}/scripts/ks_jcvi_anchors_to_pairs.py" \
    --pair-plan "${PLAN}" \
    --pair-root "${PAIR_WORK}" \
    --out-pairs "${NEW_PAIR_TABLE}" \
    --out-blocks "${NEW_PAIR_BLOCKS}" \
    --out-status "${NEW_PAIR_STATUS}"
  printf 'Missing PCAM/PPE syntenic comparisons completed.\n' > "${PAIR_ROOT}/COMPLETE"
else
  echo "[OK] Reusing missing-reference JCVI pairs: ${PAIR_ROOT}"
fi

if [[ ! -e "${KAKS_OUT}/COMPLETE" ]]; then
  "${PYTHON_BIN}" "${ROOT_DIR}/scripts/downstream_wgd_kaks_from_jcvi.py" \
    --analysis-dir "${ROOT_DIR}" \
    --jcvi-dir "${JCVI_DIR}" \
    --pair-table "${NEW_PAIR_TABLE}" \
    --outdir "${KAKS_OUT}" \
    --threads "${THREADS}" \
    --pair-mode collinear_all_pairs \
    --min-aa 50 \
    --min-codons 50 \
    --min-4d-sites 4 \
    --min-reference-block-anchors "${MIN_REFERENCE_BLOCK_ANCHORS}" \
    --kaks-cmd "${KAKS_BIN}" \
    --kaks-method YN
  require_file "${KAKS_OUT}/summary/all_merged_4dtv_kaks.tsv"
  {
    echo -e "parameter\tvalue"
    echo -e "pair_table\t${NEW_PAIR_TABLE}"
    echo -e "pair_mode\tcollinear_all_pairs"
    echo -e "kaks_method\tYN"
    echo -e "kaks_cmd\t${KAKS_BIN}"
    echo -e "min_aa\t50"
    echo -e "min_codons\t50"
    echo -e "min_4d_sites\t4"
    echo -e "min_reference_block_anchors\t${MIN_REFERENCE_BLOCK_ANCHORS}"
  } > "${KAKS_OUT}/software_and_parameters.txt"
  printf 'Missing PSP/PPA PCAM/PPE YN Ks calculations completed.\n' > "${KAKS_OUT}/COMPLETE"
else
  echo "[OK] Reusing missing-reference YN Ks: ${KAKS_OUT}"
fi

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/fig2_prepare_subgenome_reference_ks.py" \
  --gtr "${ROOT_DIR}/downstream_analysis/ks/kaks/gtr_alltaxa_collinear_all_pairs_YN/summary/all_merged_4dtv_kaks.tsv" \
  --pzi-reference "${ROOT_DIR}/downstream_analysis/ks/kaks/collinear_all_pairs_YN/summary/all_merged_4dtv_kaks.tsv" \
  --new-reference "${KAKS_OUT}/summary/all_merged_4dtv_kaks.tsv" \
  --outdir "${SOURCE_OUT}" --overwrite

if [[ "${FORCE_PLOT}" == "1" || ! -e "${PLOT_OUT}/COMPLETE" ]]; then
  "${PYTHON_BIN}" "${ROOT_DIR}/scripts/plot_fig2_subgenome_reference_ks.py" \
    --input "${SOURCE_OUT}/subgenome_reference_ks.tsv" \
    --outdir "${PLOT_OUT}" --bandwidth 0.04 --xmax 1.0 --overwrite
else
  echo "[OK] Reusing three-reference plot: ${PLOT_OUT}"
fi

for file in \
  "${PLOT_OUT}/figure2a_subgenome_reference_ks.png" \
  "${PLOT_OUT}/figure2a_subgenome_reference_ks.svg" \
  "${PLOT_OUT}/figure2a_subgenome_reference_ks.pdf" \
  "${PLOT_OUT}/figure2a_subgenome_reference_ks.tiff"; do
  require_file "${file}"
done
echo "[OK] Three-reference A/B Ks calculation and visualization completed: ${PLOT_OUT}"
