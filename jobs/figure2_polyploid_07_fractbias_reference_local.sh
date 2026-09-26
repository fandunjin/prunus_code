#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Exact local adaptation of the reference Fig. 3c workflow:
# LAST (CoGe-style CDS search) -> tandem filtering -> DAGChainer gene-order
# blocks -> Quota Align Merge -> Quota Align depth -> FraBias 100-gene windows.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
INPUT_ROOT="${INPUT_ROOT:-${ROOT_DIR}/downstream_analysis/ks/jcvi/inputs}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/downstream_analysis/figure2_polyploid/strict_wgd_mapping/04_fractbias_8track_reference_exact}"
THREADS="${THREADS:-12}"
LAST_THREADS="${LAST_THREADS:-2}"
LAST_BATCH="${LAST_BATCH:-256M}"

FIG2_ENV="${FIG2_ENV:-${CONDA_ROOT}/envs/fig2-reference}"
QA_ENV="${QA_ENV:-${CONDA_ROOT}/envs/quota-align-reference}"
PYTHON_BIN="${PYTHON_BIN:-${FIG2_ENV}/bin/python}"
QA_PYTHON="${QA_PYTHON:-${QA_ENV}/bin/python}"
LASTDB_BIN="${LASTDB_BIN:-${FIG2_ENV}/bin/lastdb}"
LASTAL_BIN="${LASTAL_BIN:-${FIG2_ENV}/bin/lastal}"
DAGCHAINER_BIN="${DAGCHAINER_BIN:-${FIG2_ENV}/bin/run_DAG_chainer.pl}"
GLPSOL_BIN="${GLPSOL_BIN:-${QA_ENV}/bin/glpsol}"
QA_ROOT="${QA_ROOT:-${ROOT_DIR}/software/quota-alignment}"
BLAST_TO_RAW="${BLAST_TO_RAW:-${QA_ROOT}/scripts/blast_to_raw.py}"
CLUSTER_UTILS="${CLUSTER_UTILS:-${QA_ROOT}/cluster_utils.py}"
QUOTA_ALIGN="${QUOTA_ALIGN:-${QA_ROOT}/quota_align.py}"

CATALOG_SCRIPT="${ROOT_DIR}/scripts/fig2_build_fractionation_catalogs.py"
RAW_TO_DAG_SCRIPT="${ROOT_DIR}/scripts/fig2_quota_raw_to_dag.py"
PARSE_SCRIPT="${ROOT_DIR}/scripts/fig2_parse_quota_anchors.py"
FRACTBIAS_SCRIPT="${ROOT_DIR}/scripts/fig2_fractbias_local_input.py"
PLOT_SCRIPT="${ROOT_DIR}/scripts/plot_fig2_panel_c_fractbias_reference_style.py"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing input/output: $1"; }
require_exec() { [[ -x "$1" ]] || die "Required executable unavailable: $1"; }

[[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] || die "THREADS must be a positive integer"
[[ "${LAST_THREADS}" =~ ^[1-9][0-9]*$ ]] || die "LAST_THREADS must be a positive integer"
[[ ! -e "${OUT_ROOT}" ]] || die "Output directory already exists: ${OUT_ROOT}"
for script in "${CATALOG_SCRIPT}" "${RAW_TO_DAG_SCRIPT}" "${PARSE_SCRIPT}" "${FRACTBIAS_SCRIPT}" "${PLOT_SCRIPT}" \
              "${BLAST_TO_RAW}" "${CLUSTER_UTILS}" "${QUOTA_ALIGN}"; do
    require_file "${script}"
done
for executable in "${PYTHON_BIN}" "${QA_PYTHON}" "${LASTDB_BIN}" "${LASTAL_BIN}" \
                  "${DAGCHAINER_BIN}" "${GLPSOL_BIN}"; do
    require_exec "${executable}"
done

CATALOG_DIR="${OUT_ROOT}/01_catalog"
LAST_DIR="${OUT_ROOT}/02_last"
DAG_DIR="${OUT_ROOT}/03_dagchainer"
QUOTA_DIR="${OUT_ROOT}/04_quota_align"
FRACTBIAS_DIR="${OUT_ROOT}/05_fractbias"
PANEL_DIR="${OUT_ROOT}/06_panel_c"
mkdir -p "${LAST_DIR}" "${DAG_DIR}" "${QUOTA_DIR}"

"${PYTHON_BIN}" "${CATALOG_SCRIPT}" --inputs-root "${INPUT_ROOT}" --outdir "${CATALOG_DIR}"
require_file "${CATALOG_DIR}/query8.pep"
require_file "${CATALOG_DIR}/tracks/GTR.pep"

cp "${CATALOG_DIR}/query8.bed" "${LAST_DIR}/query8.bed"
cp "${CATALOG_DIR}/tracks/GTR.bed" "${LAST_DIR}/GTR.bed"

LAST_DB="${LAST_DIR}/GTR_pep"
LAST_OUT="${LAST_DIR}/query8_vs_GTR.last.blasttab"
"${LASTDB_BIN}" -p -c "${LAST_DB}" "${CATALOG_DIR}/tracks/GTR.pep"
for track in PZI_A PZI_B PSP_A PSP_B PPA_A PPA_B PCAM PPE; do
    track_output="${LAST_DIR}/${track}_vs_GTR.last.blasttab"
    "${LASTAL_BIN}" -p BL62 -P "${LAST_THREADS}" -i"${LAST_BATCH}" -f BlastTab \
        "${LAST_DB}" "${CATALOG_DIR}/tracks/${track}.pep" \
        > "${track_output}" 2> "${LAST_DIR}/${track}.lastal.stderr.log"
    require_file "${track_output}"
done
awk '!/^#/' "${LAST_DIR}"/*_vs_GTR.last.blasttab > "${LAST_OUT}"
require_file "${LAST_OUT}"

PATH="${QA_ENV}/bin:${PATH}" "${QA_PYTHON}" "${BLAST_TO_RAW}" \
    "${LAST_OUT}" \
    --qbed="${LAST_DIR}/query8.bed" \
    --sbed="${LAST_DIR}/GTR.bed" \
    --no_strip_names \
    --tandem_Nmax=10 \
    --write-filtered-blast \
    > "${LAST_DIR}/query8_vs_GTR.last.tandem_filtered.blasttab" \
    2> "${LAST_DIR}/blast_to_raw.stderr.log"

RAW_HITS="${LAST_DIR}/query8_vs_GTR.last.raw"
QUERY_REDUCED_BED="${LAST_DIR}/query8.nolocaldups.bed"
TARGET_REDUCED_BED="${LAST_DIR}/GTR.nolocaldups.bed"
require_file "${RAW_HITS}"
require_file "${QUERY_REDUCED_BED}"
require_file "${TARGET_REDUCED_BED}"

DAG_INPUT="${DAG_DIR}/query8_vs_GTR.geneorder.dag"
"${PYTHON_BIN}" "${RAW_TO_DAG_SCRIPT}" \
    --raw "${RAW_HITS}" \
    --query-bed "${QUERY_REDUCED_BED}" \
    --target-bed "${TARGET_REDUCED_BED}" \
    --out "${DAG_INPUT}" \
    --stats "${DAG_DIR}/dag_input_stats.tsv"

"${DAGCHAINER_BIN}" -i "${DAG_INPUT}" -D 20 -g 10 -A 5 \
    > "${DAG_DIR}/dagchainer.stdout.log" 2> "${DAG_DIR}/dagchainer.stderr.log"
DAG_BLOCKS="${DAG_INPUT}.aligncoords"
require_file "${DAG_BLOCKS}"

INITIAL_QA="${QUOTA_DIR}/query8_vs_GTR.initial.qa"
PATH="${QA_ENV}/bin:${PATH}" "${QA_PYTHON}" "${CLUSTER_UTILS}" \
    --format=dag --log_evalue "${DAG_BLOCKS}" "${INITIAL_QA}" \
    > "${QUOTA_DIR}/cluster_utils.stdout.log" 2> "${QUOTA_DIR}/cluster_utils.stderr.log"
require_file "${INITIAL_QA}"

PATH="${QA_ENV}/bin:${PATH}" "${QA_PYTHON}" "${QUOTA_ALIGN}" \
    --Dm=20 --merge "${INITIAL_QA}" \
    > "${QUOTA_DIR}/quota_merge.stdout.log" 2> "${QUOTA_DIR}/quota_merge.stderr.log"
MERGED_QA="${INITIAL_QA}.merged"
require_file "${MERGED_QA}"

read -r QA_LEFT QA_LEFT_POS QA_RIGHT QA_RIGHT_POS QA_SCORE < <(awk '!/^#/ {print; exit}' "${MERGED_QA}")
if [[ "${QA_LEFT}" == GTR__* && "${QA_RIGHT}" != GTR__* ]]; then
    FILE_AXIS_QUOTA="1:8"
elif [[ "${QA_LEFT}" != GTR__* && "${QA_RIGHT}" == GTR__* ]]; then
    FILE_AXIS_QUOTA="8:1"
else
    die "Cannot identify GTR/query axes in merged QA: ${QA_LEFT} versus ${QA_RIGHT}"
fi
PATH="${QA_ENV}/bin:${PATH}" "${QA_PYTHON}" "${QUOTA_ALIGN}" \
    --Nm=40 --quota="${FILE_AXIS_QUOTA}" --solver=GLPK "${MERGED_QA}" \
    > "${QUOTA_DIR}/quota_depth.stdout.log" 2> "${QUOTA_DIR}/quota_depth.stderr.log"
FILTERED_QA="${MERGED_QA}.filtered"
require_file "${FILTERED_QA}"

NORMALIZED_DIR="${QUOTA_DIR}/normalized"
"${PYTHON_BIN}" "${PARSE_SCRIPT}" \
    --quota-output "${FILTERED_QA}" \
    --catalog "${CATALOG_DIR}" \
    --query-bed "${QUERY_REDUCED_BED}" \
    --target-bed "${TARGET_REDUCED_BED}" \
    --outdir "${NORMALIZED_DIR}"
require_file "${NORMALIZED_DIR}/quota_filtered_anchors.tsv"

"${PYTHON_BIN}" "${FRACTBIAS_SCRIPT}" \
    --anchors "${NORMALIZED_DIR}/quota_filtered_anchors.tsv" \
    --catalog "${CATALOG_DIR}" \
    --window-size 100 \
    --outdir "${FRACTBIAS_DIR}"
require_file "${FRACTBIAS_DIR}/fractbias_windows.tsv"

"${PYTHON_BIN}" "${PLOT_SCRIPT}" \
    --fractbias "${FRACTBIAS_DIR}/fractbias_windows.tsv" \
    --outdir "${PANEL_DIR}"

{
    echo "alignment=protein-LAST on frozen gene catalogs (CoGe SynMap protein-feature mode)"
    echo "lastdb=-p -c; lastal=-p BL62 -P ${LAST_THREADS} -i${LAST_BATCH} -f BlastTab with LAST default score threshold; eight query tracks run separately and concatenated"
    echo "local_duplicate_filter=tandem_Nmax 10"
    echo "dagchainer=gene order; D 20; g 10; A 5; default maximum E-value 1e-5"
    echo "quota_align_merge=Dm 20"
    echo "quota_align_depth=logical query:GTR 8:1; file-axis ${FILE_AXIS_QUOTA}; Nm 40; GLPK"
    echo "fractbias=all GTR target genes; overlapping windows of 100 genes"
    echo "tracks=PZI_A,PZI_B,PSP_A,PSP_B,PPA_A,PPA_B,PCAM,PPE"
    echo "display_pairs=C06,C02,C03"
    echo "boundary=quota-filtered anchor retention; absence is not a validated gene loss"
} > "${OUT_ROOT}/software_and_parameters.txt"

{
    "${LASTAL_BIN}" -V 2>&1 | head -1
    "${DAGCHAINER_BIN}" -h 2>&1 | head -1 || true
    "${GLPSOL_BIN}" --version 2>&1 | head -1
} > "${OUT_ROOT}/tool_versions.txt"
touch "${OUT_ROOT}/COMPLETE"
echo "[OK] Exact reference-adapted Fig. 2c completed: ${OUT_ROOT}"
