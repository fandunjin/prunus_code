#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/analysis}"
SCRIPTS_DIR="${SCRIPTS_DIR:-${ROOT_DIR}/scripts}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
WGD_ENV="${WGD_ENV:-wgd_kaks}"
THREADS="${THREADS:-32}"
JCVI_DIR="${ANALYSIS_DIR}/karyotype/jcvi"

# Publication default: one representative gene pair per JCVI syntenic block.
# Use PAIR_MODE=all_pairs only for a much heavier sensitivity analysis.
PAIR_MODE="${PAIR_MODE:-block_representative}"
MAX_PAIRS_PER_COMPARISON="${MAX_PAIRS_PER_COMPARISON:-0}"
MIN_AA="${MIN_AA:-50}"
MIN_CODONS="${MIN_CODONS:-50}"
MIN_4D_SITES="${MIN_4D_SITES:-4}"
KAKS_METHOD="${KAKS_METHOD:-YN}"
WGD_ROOT="${ANALYSIS_DIR}/wgd_kaks"
WGD_RUN_NAME="${WGD_RUN_NAME:-${PAIR_MODE}_${KAKS_METHOD}}"
WGD_DIR="${WGD_ROOT}/${WGD_RUN_NAME}"

mkdir -p "${WGD_ROOT}" "${WGD_DIR}" "${ANALYSIS_DIR}/logs"

set +u
source "${CONDA_SH}"
conda activate "${WGD_ENV}"
set -u

for exe in python mafft pal2nal.pl KaKs_Calculator; do
    if ! command -v "${exe}" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found in WGD_ENV=${WGD_ENV}: ${exe}" >&2
        exit 2
    fi
done

if [ ! -s "${SCRIPTS_DIR}/downstream_wgd_kaks_from_jcvi.py" ]; then
    echo "[ERROR] Missing helper script: ${SCRIPTS_DIR}/downstream_wgd_kaks_from_jcvi.py" >&2
    exit 2
fi

for table in \
    "${JCVI_DIR}/summary/jcvi_syntenic_pairs.tsv" \
    "${JCVI_DIR}/summary/jcvi_pair_status.tsv"; do
    if [ ! -s "${table}" ]; then
        echo "[ERROR] Missing JCVI summary table: ${table}" >&2
        echo "[ERROR] Run downstream_09_karyotype_jcvi.sh successfully before 07." >&2
        exit 2
    fi
done

python "${SCRIPTS_DIR}/downstream_wgd_kaks_from_jcvi.py" \
    --analysis-dir "${ANALYSIS_DIR}" \
    --jcvi-dir "${JCVI_DIR}" \
    --outdir "${WGD_DIR}" \
    --threads "${THREADS}" \
    --pair-mode "${PAIR_MODE}" \
    --max-pairs-per-comparison "${MAX_PAIRS_PER_COMPARISON}" \
    --min-aa "${MIN_AA}" \
    --min-codons "${MIN_CODONS}" \
    --min-4d-sites "${MIN_4D_SITES}" \
    --kaks-method "${KAKS_METHOD}"

cat > "${WGD_DIR}/README.txt" <<'EOF'
downstream_07_wgd_kaks.sh performs synteny-based WGD/4DTV/KaKs analysis from
downstream_09 JCVI outputs.

Input source:
  analysis/karyotype/jcvi/summary/jcvi_syntenic_pairs.tsv
  analysis/karyotype/jcvi/inputs/<label>/<label>.pep/.cds

Default strategy:
  PAIR_MODE=block_representative
    Use one representative syntenic gene pair from each JCVI block. This avoids
    overweighting dense blocks and makes the first publication-grade run
    tractable. Use PAIR_MODE=all_pairs for a slower sensitivity analysis.

Alignment and Ka/Ks strategy:
  1. Align each protein pair with MAFFT.
  2. Back-translate protein alignment to codon alignment with PAL2NAL -nogap.
  3. Calculate pairwise 4DTV from fourfold-degenerate codon columns.
  4. Write one AXT file per comparison.
  5. Run KaKs_Calculator with the YN method by default.

Main outputs:
  summary/syntenic_pair_plan.tsv
  summary/all_pair_alignment_stats.tsv
  summary/all_merged_4dtv_kaks.tsv
  summary/wgd_kaks_comparison_summary.tsv
  pairs/<comparison>.alignment_stats.tsv
  axt/<comparison>.axt
  kaks/<comparison>.KaKs_Calculator.YN.tsv
  kaks/<comparison>.merged_4dtv_kaks.tsv
  figures/ks_distribution_by_comparison.svg
  figures/4dtv_distribution_by_comparison.svg
  figures/kaks_ratio_distribution_by_comparison.svg

Interpretation cautions:
  Ka/Ks > 1 candidates require manual alignment inspection before biological
  interpretation. Very high Ks or saturated 4DTV values should be filtered or
  interpreted as saturated divergence.
EOF

echo "[INFO] WGD/4DTV/KaKs analysis finished: ${WGD_DIR}"
