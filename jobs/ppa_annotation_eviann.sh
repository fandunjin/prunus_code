#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
[[ -s "${RUN_SCRIPT}" ]] || { echo "ERROR: PPA EviAnn run script missing: ${RUN_SCRIPT}" >&2; exit 1; }
exec bash "${RUN_SCRIPT}"
