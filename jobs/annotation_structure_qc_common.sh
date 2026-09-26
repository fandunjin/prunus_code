#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# Shared runner for the PZI/PSP/PPA GFF3 structure QC jobs.

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}${CONDA_ROOT}/etc/profile.d/conda.sh}"
QC_SCRIPT="${QC_SCRIPT:-${ROOT_DIR}/scripts/gff3_structure_qc.py}"
THREADS="${THREADS:-8}"

die() { echo "ERROR: $*" >&2; exit 1; }
require_file() { [[ -s "$1" ]] || die "$2 does not exist or is empty: $1"; }

run_structure_qc() {
    : "${SPECIES:?SPECIES is required}"
    : "${GFF3:?GFF3 is required}"
    : "${OUT_DIR:?OUT_DIR is required}"
    require_file "$CONDA_SH" "Conda initialization script"
    require_file "$QC_SCRIPT" "GFF3 structure QC script"
    require_file "$GFF3" "${SPECIES} final annotation GFF3"
    mkdir -p "$OUT_DIR" "${ROOT_DIR}/logs"

    source "$CONDA_SH"
    conda activate "$PYTHON_ENV"
    command -v python3 >/dev/null 2>&1 || die "python3 is not available in ${PYTHON_ENV}"
    if [[ -n "${BAM_SPECS[*]:-}" ]]; then
        command -v samtools >/dev/null 2>&1 || die "samtools is required for splice support"
    fi

    args=(python3 "$QC_SCRIPT" --gff3 "$GFF3" --prefix "$SPECIES" --out-dir "$OUT_DIR" --threads "$THREADS")
    if [[ -n "${TE_GFF3:-}" ]]; then
        require_file "$TE_GFF3" "${SPECIES} EDTA TE GFF3"
        args+=(--te-gff3 "$TE_GFF3")
    fi
    if [[ -n "${BAM_SPECS[*]:-}" ]]; then
        for spec in "${BAM_SPECS[@]}"; do
            label="${spec%%=*}"
            path="${spec#*=}"
            [[ "$label" != "$spec" && -n "$path" ]] || die "BAM_SPECS entries must be LABEL=PATH"
            require_file "$path" "${SPECIES} splice BAM (${label})"
            [[ -s "${path}.bai" || -s "${path%.bam}.bai" ]] || die "BAM index missing for ${path}"
            args+=(--bam "$label=$path")
        done
    fi
    echo "=== ${SPECIES} GFF3 structure QC ==="
    echo "GFF3: $GFF3"
    echo "Output: $OUT_DIR"
    printf 'BAM: %s\n' "${BAM_SPECS[@]:-none}"
    "${args[@]}"
    require_file "$OUT_DIR/gene_structure_summary.tsv" "${SPECIES} structure summary"
    require_file "$OUT_DIR/gene_structure_metrics.tsv" "${SPECIES} structure metrics"
    require_file "$OUT_DIR/intron_details.tsv" "${SPECIES} intron details"
    {
        echo "python=$(python3 --version 2>&1)"
        echo "gff3_structure_qc=${QC_SCRIPT}"
        if command -v samtools >/dev/null 2>&1; then samtools --version | head -n 1; fi
    } > "$OUT_DIR/tool_versions.txt"
    echo "${SPECIES} structure QC completed: $OUT_DIR"
}
