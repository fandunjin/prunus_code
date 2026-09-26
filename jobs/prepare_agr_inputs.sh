#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/downstream_analysis/orthofinder/species/results}"
KARYOTYPE_ROOT="${KARYOTYPE_ROOT:-${ROOT_DIR}/downstream_analysis/karyotype}"
OUT_DIR="${OUT_DIR:-${KARYOTYPE_ROOT}/agr}"
PREPARE_SCRIPT="${ROOT_DIR}/scripts/prepare_agr_inputs.py"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
PYTHON_ENV="${PYTHON_ENV:-jcvi}"
FORCE="${FORCE:-0}"

# These are all available in the completed OrthoFinder/JCVI analyses and
# allow both the focal analysis and broader Prunus sensitivity analyses.
SPECIES="${SPECIES:-pzi_copy1 psp ppa ppe pdu pmu_refseq pcam pspe}"
read -r -a species_array <<< "$SPECIES"

mkdir -p "$OUT_DIR"

if [[ "$FORCE" != "1" \
      && -s "$OUT_DIR/orthogroups.csv" \
      && -s "$OUT_DIR/agr_gff.csv" \
      && -s "$OUT_DIR/agr_chr.csv" ]]; then
    echo "[INFO] AGR input tables already exist; skip. Set FORCE=1 to rebuild."
    echo "$OUT_DIR/orthogroups.csv"
    echo "$OUT_DIR/agr_gff.csv"
    echo "$OUT_DIR/agr_chr.csv"
    exit 0
fi

if [[ -z "${ORTHOGROUPS_TSV:-}" ]]; then
    latest_orthogroups=""
    while IFS= read -r candidate; do
        [[ -s "$candidate" ]] || continue
        if [[ -z "$latest_orthogroups" || "$candidate" -nt "$latest_orthogroups" ]]; then
            latest_orthogroups="$candidate"
        fi
    done < <(find "$OF_RESULTS_ROOT" -type f -path '*/Orthogroups/Orthogroups.tsv' -print)
    ORTHOGROUPS_TSV="$latest_orthogroups"
fi

[[ -n "${ORTHOGROUPS_TSV:-}" && -s "$ORTHOGROUPS_TSV" ]] || {
    echo "[ERROR] No completed OrthoFinder Orthogroups.tsv found under $OF_RESULTS_ROOT" >&2
    exit 1
}

if [[ -z "${JCVI_RUN_DIR:-}" ]]; then
    latest_jcvi=""
    for candidate in "$KARYOTYPE_ROOT"/run_*; do
        [[ -d "$candidate/inputs" ]] || continue
        candidate_complete=1
        for species in "${species_array[@]}"; do
            if [[ ! -d "$candidate/inputs/$species" ]]; then
                candidate_complete=0
                break
            fi
        done
        [[ "$candidate_complete" == "1" ]] || continue
        latest_jcvi="$candidate"
    done
    JCVI_RUN_DIR="$latest_jcvi"
fi

[[ -n "${JCVI_RUN_DIR:-}" && -d "$JCVI_RUN_DIR/inputs" ]] || {
    echo "[ERROR] No JCVI run with normalized inputs was found under $KARYOTYPE_ROOT" >&2
    exit 1
}

[[ -s "$PREPARE_SCRIPT" ]] || {
    echo "[ERROR] Missing converter: $PREPARE_SCRIPT" >&2
    exit 1
}

[[ -s "$CONDA_SH" ]] || {
    echo "[ERROR] Missing conda initialization script: $CONDA_SH" >&2
    exit 1
}

source "$CONDA_SH"
conda activate "$PYTHON_ENV"
command -v python >/dev/null 2>&1 || {
    echo "[ERROR] python is unavailable in $PYTHON_ENV" >&2
    exit 1
}

echo "=== Prepare AGR inputs ==="
echo "OrthoFinder : $ORTHOGROUPS_TSV"
echo "JCVI run    : $JCVI_RUN_DIR"
echo "Output      : $OUT_DIR"
echo "Species     : ${species_array[*]}"

python "$PREPARE_SCRIPT" \
    --orthogroups "$ORTHOGROUPS_TSV" \
    --jcvi-input-root "$JCVI_RUN_DIR/inputs" \
    --genome-root "$JCVI_RUN_DIR/stage/canonical" \
    --outdir "$OUT_DIR" \
    --min-match-rate 0.90 \
    --species "${species_array[@]}"

echo "=== AGR input QC ==="
cat "$OUT_DIR/species_summary.tsv"
echo "[OK] AGR inputs completed: $OUT_DIR"
