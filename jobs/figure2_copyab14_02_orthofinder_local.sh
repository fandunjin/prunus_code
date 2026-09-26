#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
INPUT_DIR="${INPUT_DIR:-${ROOT_DIR}/0903/analysisphylogenetics/01_frozen_inputs_copyAB14/pep}"
OUT_PARENT="${OUT_PARENT:-${ROOT_DIR}/0903/analysisphylogenetics/02_orthofinder_copyAB14}"
THREADS="${THREADS:-12}"
ANALYSIS_THREADS="${ANALYSIS_THREADS:-4}"
SEARCH_PROGRAM="${SEARCH_PROGRAM:-diamond}"

die() { echo "[ERROR] $*" >&2; exit 2; }
[[ -d "${INPUT_DIR}" ]] || die "Missing input directory: ${INPUT_DIR}; run figure2_copyab14_01_prepare_local.sh first"
command -v orthofinder >/dev/null 2>&1 || die "orthofinder is unavailable on PATH"
command -v "${SEARCH_PROGRAM}" >/dev/null 2>&1 || die "${SEARCH_PROGRAM} is unavailable on PATH"
mkdir -p "${OUT_PARENT}"
[[ -z "${RUN_ID:-}" ]] && RUN_ID="$(date +%Y%m%d_%H%M%S)"
OUT_DIR="${OUT_PARENT}/Results_${RUN_ID}"
[[ ! -e "${OUT_DIR}" ]] || die "Output exists: ${OUT_DIR}"

orthofinder -f "${INPUT_DIR}" -S "${SEARCH_PROGRAM}" -t "${THREADS}" -a "${ANALYSIS_THREADS}" -o "${OUT_DIR}"
RESULT_DIR="$(find "${OUT_DIR}" -type f -path '*/Orthogroups/Orthogroups.tsv' -print -quit | sed 's#/Orthogroups/Orthogroups.tsv##')"
[[ -n "${RESULT_DIR}" ]] || die "OrthoFinder output was not found under ${OUT_DIR}"
printf '%s\n' "${RESULT_DIR}" > "${OUT_PARENT}/latest.txt"
echo "[OK] CopyAB14 OrthoFinder completed: ${RESULT_DIR}"
