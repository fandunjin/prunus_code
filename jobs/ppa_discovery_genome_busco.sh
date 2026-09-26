#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/data/ppa/busco_genome_discovery}"
RUN_NAME="${RUN_NAME:-ppa_discovery_genome_eudicotyledons_odb12}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
BUSCO_ENV="${BUSCO_ENV:-busco}"
BUSCO_LINEAGE="${BUSCO_LINEAGE:-${HOME}/data/busco_db/eudicotyledons_odb12/}"
THREADS="${THREADS:-24}"
FORCE="${FORCE:-0}"

die() { echo "ERROR: $*" >&2; exit 1; }
require_file() { [[ -s "$1" ]] || die "$2 does not exist or is empty: $1"; }
require_cmd() { command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"; }

case "$FORCE" in 0|1) ;; *) die "FORCE must be 0 or 1" ;; esac
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || die "THREADS must be a positive integer"
require_file "$PPA_GENOME" "PPA discovery genome FASTA"
require_file "$BUSCO_LINEAGE" "BUSCO lineage dataset"
require_file "$CONDA_SH" "Conda initialization script"
mkdir -p "$OUT_DIR" "${ROOT_DIR}/logs"

set +u
source "$CONDA_SH"
set -u
conda activate "$BUSCO_ENV"
require_cmd busco

RUN_DIR="${OUT_DIR}/${RUN_NAME}"
SUMMARY="$(find "$RUN_DIR" -maxdepth 3 -type f -name 'short_summary*.txt' -size +0c 2>/dev/null | head -n 1 || true)"
if [[ -n "$SUMMARY" && "$FORCE" == "0" ]]; then
    echo "Existing complete BUSCO result found; skipping: $SUMMARY"
    exit 0
fi

echo "=== BUSCO: PPA discovery primary genome ==="
echo "Genome  : $PPA_GENOME"
echo "Mode    : genome"
echo "Lineage : $BUSCO_LINEAGE"
echo "Threads : $THREADS"
echo "Output  : $RUN_DIR"
echo "Force   : $FORCE"

BUSCO_ARGS=(
    -i "$PPA_GENOME"
    -m genome
    -l "$BUSCO_LINEAGE"
    -c "$THREADS"
    --offline
    -o "$RUN_NAME"
    --out_path "$OUT_DIR"
)
if [[ "$FORCE" == "1" ]]; then
    BUSCO_ARGS+=(--force)
elif [[ -d "$RUN_DIR" ]]; then
    BUSCO_ARGS+=(--restart)
fi
busco "${BUSCO_ARGS[@]}"

SUMMARY="$(find "$RUN_DIR" -maxdepth 3 -type f -name 'short_summary*.txt' -size +0c | head -n 1 || true)"
[[ -n "$SUMMARY" ]] || die "BUSCO completed without a non-empty short summary in $RUN_DIR"
echo "BUSCO completed. Summary: $SUMMARY"
