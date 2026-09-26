#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/haploid_repartition/ppa_copy1_minimap}"
PPA_GENOME="${PPA_GENOME:-${ROOT_DIR}/data/ppa/GCA_964059315.1_drPruPadu1.1_genomic.fna.gz}"
PZI_COPY1="${PZI_COPY1:-${ROOT_DIR}/final_genome/pzi_copy1.fa}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
MINIMAP_ENV="${MINIMAP_ENV:-minimap}"
PREPARE_PPA_SCRIPT="${ROOT_DIR}/scripts/prepare_ppa_primary_chr16.py"
SUMMARIZE_SCRIPT="${ROOT_DIR}/scripts/summarize_minimap_paf_pairs.py"
THREADS="${THREADS:-24}"
MINIMAP_PRESET="${MINIMAP_PRESET:-asm5}"
MIN_MAPQ="${MIN_MAPQ:-20}"
MIN_ALN_LEN="${MIN_ALN_LEN:-10000}"

die() { echo "[ERROR] $*" >&2; exit 2; }
req() { [[ -s "$1" ]] || die "Missing or empty $2: $1"; }

req "${PPA_GENOME}" "PPA primary genome"
req "${PZI_COPY1}" "current PZI copy1 FASTA"
req "${CONDA_SH}" "conda initialization"
req "${PREPARE_PPA_SCRIPT}" "PPA primary chromosome preparation script"
req "${SUMMARIZE_SCRIPT}" "PAF summarization script"

mkdir -p "${WORK_DIR}" "${ROOT_DIR}/logs"
source "${CONDA_SH}"
conda activate "${MINIMAP_ENV}"
command -v minimap2 >/dev/null 2>&1 || die "minimap2 is unavailable"
command -v python3 >/dev/null 2>&1 && PYTHON_BIN=python3 || PYTHON_BIN=python

PPA_MAP="${WORK_DIR}/ppa_primary_chr16.rename_map.tsv"
PAF="${WORK_DIR}/ppa_vs_pzi_copy1.asm5.paf"
FILTERED_PAF="${WORK_DIR}/ppa_vs_pzi_copy1.asm5.mapq${MIN_MAPQ}.min${MIN_ALN_LEN}.paf"
SUMMARY="${WORK_DIR}/ppa_vs_pzi_copy1.asm5.alignment_summary.tsv"
BEST_HITS="${WORK_DIR}/ppa_vs_pzi_copy1.asm5.best_pzi_copy1_hit_per_ppa.tsv"

"${PYTHON_BIN}" "${PREPARE_PPA_SCRIPT}" \
    --genome "${PPA_GENOME}" \
    --out-map "${PPA_MAP}"

echo "=== New PPA vs PZI copy1 minimap2 alignment ==="
echo "PPA target : ${PPA_GENOME}"
echo "PZI query  : ${PZI_COPY1}"
echo "Output     : ${WORK_DIR}"
echo "Preset     : ${MINIMAP_PRESET}"

minimap2 -x "${MINIMAP_PRESET}" -t "${THREADS}" \
    "${PPA_GENOME}" "${PZI_COPY1}" > "${PAF}.tmp"
mv "${PAF}.tmp" "${PAF}"

"${PYTHON_BIN}" "${SUMMARIZE_SCRIPT}" \
    --paf "${PAF}" \
    --ppa-map "${PPA_MAP}" \
    --out-filtered-paf "${FILTERED_PAF}" \
    --out-summary "${SUMMARY}" \
    --out-best "${BEST_HITS}" \
    --min-mapq "${MIN_MAPQ}" \
    --min-alignment-len "${MIN_ALN_LEN}"

for output in "${PPA_MAP}" "${PAF}" "${FILTERED_PAF}" "${SUMMARY}" "${BEST_HITS}"; do
    req "${output}" "new PPA-PZI output"
done
echo "[OK] New PPA-PZI copy1 outputs: ${WORK_DIR}"
