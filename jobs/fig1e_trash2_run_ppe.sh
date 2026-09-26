#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# fig1e_trash2_run_node_ppe.sh
# TRASH_2 for ppe on the node node pool (TARGET hard-coded inside).
#   cd ${PROJECT_ROOT}/fig1e_trash2_run_node_ppe.sh
# Optional: CLEAN=1 removes an existing (failed/partial) output dir first.
INPUT_DIR="${FIG_ROOT}/inputs/genomes"
OUT_ROOT="${FIG_ROOT}/trash"
TARGET="ppe"
THREADS="${THREADS:-24}"
CLEAN="${CLEAN:-0}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"

die() { echo "[ERROR] $*" >&2; exit 2; }
require_file() { [[ -s "$1" ]] || die "Missing or empty $2: $1"; }

INPUT="${INPUT_DIR}/${TARGET}.fa"
OUT_DIR="${OUT_ROOT}/${TARGET}"
TRASH_SRC="${FIG_ROOT}/tools/TRASH_2/src"

require_file "${CONDA_SH}" "conda initialization"
require_file "${INPUT}" "${TARGET} genome FASTA"
require_file "${TRASH_SRC}/TRASH.R" "TRASH.R (run fig1e_trash2_setup.sh first)"
mkdir -p "${OUT_ROOT}" "${ROOT_DIR}/logs"

if [[ -s "${OUT_DIR}/trash.completed.tsv" && "${CLEAN}" != "1" ]]; then
    echo "[SKIP] ${TARGET} already completed: ${OUT_DIR}/trash.completed.tsv"
    exit 0
fi

if [[ "${CLEAN}" == "1" && -d "${OUT_DIR}" ]]; then
    echo "[CLEAN] removing previous output dir: ${OUT_DIR}"
    rm -rf "${OUT_DIR}"
fi
mkdir -p "${OUT_DIR}"

source "${CONDA_SH}"
conda activate fig1e
grep -q doMC "${TRASH_SRC}/main.R" || die "TRASH_2 main.R is missing the doMC patch"

echo "=== fig1e TRASH_2 (node): ${TARGET} ==="
echo "Input : ${INPUT}"
echo "Output: ${OUT_DIR}"
echo "Cores : ${THREADS}"

cd "${TRASH_SRC}"
Rscript ./TRASH.R -f "${INPUT}" -o "${OUT_DIR}" -p "${THREADS}"

for suffix in ".fa_repeats.csv" ".fa_repeats_with_seq.csv" ".fa_arrays.csv" ".fa_repeats.gff" ".fa_arrays.gff" ".fa_run_time.csv"; do
    require_file "${OUT_DIR}/${TARGET}${suffix}" "${TARGET} TRASH output ${suffix}"
done

printf 'target\tinput\trepeats_csv\trepeats_with_seq\tarrays_csv\trepeats_gff\trun_time\tstatus\n' > "${OUT_DIR}/trash.completed.tsv"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\tpass\n' "${TARGET}" "${INPUT}" \
    "${OUT_DIR}/${TARGET}.fa_repeats.csv" \
    "${OUT_DIR}/${TARGET}.fa_repeats_with_seq.csv" \
    "${OUT_DIR}/${TARGET}.fa_arrays.csv" \
    "${OUT_DIR}/${TARGET}.fa_repeats.gff" \
    "${OUT_DIR}/${TARGET}.fa_run_time.csv" >> "${OUT_DIR}/trash.completed.tsv"

echo "[OK] ${TARGET} TRASH_2 complete: ${OUT_DIR}"
