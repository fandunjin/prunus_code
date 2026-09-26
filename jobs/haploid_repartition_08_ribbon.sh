#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/downstream_analysis/repartition/ribbon}"
JCVI_ROOT="${JCVI_ROOT:-${ROOT_DIR}/downstream_analysis/repartition/work/jcvi}"
INPUT_ROOT="${INPUT_ROOT:-${ROOT_DIR}/downstream_analysis/repartition/work/jcvi_inputs}"
SCRIPT="${ROOT_DIR}/scripts/prepare_repartitioned_ribbon.py"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
JCVI_ENV="${JCVI_ENV:-jcvi}"
FIGSIZE="${FIGSIZE:-12x12}"

die() { echo "[ERROR] $*" >&2; exit 2; }
req() { [[ -s "$1" ]] || die "Missing or empty $2: $1"; }
req "${SCRIPT}" "ribbon preparation script"
req "${CONDA_SH}" "conda initialization"
mkdir -p "${OUT_ROOT}" "${ROOT_DIR}/logs"
source "${CONDA_SH}"
conda activate "${JCVI_ENV}"
command -v python >/dev/null 2>&1 || die "python is unavailable in ${JCVI_ENV}"
python "${SCRIPT}" --jcvi-root "${JCVI_ROOT}" --input-root "${INPUT_ROOT}" --out-root "${OUT_ROOT}"

cd "${OUT_ROOT}"
mkdir -p figures
for format in pdf svg png; do
    python -m jcvi.graphics.karyotype \
        repartitioned.seqids repartitioned.layout \
        --figsize="${FIGSIZE}" --format="${format}" --notex --nocircles \
        --keep-chrlabels --chrstyle=roundrect --shadestyle=curve \
        --outfile="figures/pzi_repartition_ribbon.${format}"
done
for output in figures/pzi_repartition_ribbon.pdf figures/pzi_repartition_ribbon.svg; do
    req "${OUT_ROOT}/${output}" "JCVI ribbon output"
done
echo "[OK] Native JCVI ribbon: ${OUT_ROOT}/figures"
