#!/usr/bin/env bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
set -euo pipefail
ROOT_DIR="${PROJECT_ROOT}"
cd "${ROOT_DIR}"
export PATH="${CONDA_ROOT}/bin:${PATH}"
export ASTRAL_JAR="${ROOT_DIR}/software/astral_5.7.8/Astral/astral.5.7.8.jar"
export THREADS="${THREADS:-4}"
exec bash "${ROOT_DIR}/jobs/repartition_phylo_05_resume_local.sh"
