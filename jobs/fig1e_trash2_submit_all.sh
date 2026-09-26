#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# fig1e_trash2_submit_all.sh
# Preflight + submit all 14 fig1e TRASH_2 jobs. Run on the HPC login node.
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/fig1e_satdna}"
INPUT_DIR="${INPUT_DIR:-${FIG_ROOT}/inputs/genomes}"

TARGETS=(pziA1 pziA2 pziB1 pziB2 pspA1 pspA2 pspB1 pspB2 ppaA ppaB ppe pcam gtr rch)

echo "=== input preflight (${INPUT_DIR}) ==="
missing=0
for t in "${TARGETS[@]}"; do
    if [[ -s "${INPUT_DIR}/${t}.fa" ]]; then
        printf '  %-6s OK   %s\n' "${t}" "$(du -h "${INPUT_DIR}/${t}.fa" | cut -f1)"
    else
        printf '  %-6s MISSING %s\n' "${t}" "${INPUT_DIR}/${t}.fa"
        missing=$((missing + 1))
    fi
done
[[ "${missing}" -eq 0 ]] || { echo "[ERROR] ${missing} input(s) missing; upload or symlink them first" >&2; exit 2; }

[[ -s "${FIG_ROOT}/tools/TRASH_2/src/TRASH.R" ]] || { echo "[ERROR] TRASH_2 not set up; run: bash jobs/fig1e_trash2_setup.sh" >&2; exit 2; }

cd "${ROOT_DIR}"
for t in "${TARGETS[@]}"; do
    export TARGET="${t}"
    echo ">>> submit ${t}"
    jsub < jobs/fig1e_trash2_run.sh
done
echo "=== 14 jobs submitted ==="
