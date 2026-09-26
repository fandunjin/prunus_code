#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
PROTEIN="${PROTEIN:-${DATA_DIR}/annotation/genes/psp32_eviann/postprocess/psp32.eviann.protein.fa}"; OUT_BASE="${DATA_DIR}/annotation/genes/psp32_eviann/busco"; LINEAGE="${BUSCO_LINEAGE:-${HOME}/data/busco_db/eudicotyledons_odb12/}"; THREADS="${THREADS:-24}"; RUN_NAME="psp32_eviann_protein_eudicots"; CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
[[ -s "$PROTEIN" ]] || { echo "Missing protein: $PROTEIN" >&2; exit 1; }; [[ -s "$LINEAGE" ]] || { echo "Missing BUSCO lineage: $LINEAGE" >&2; exit 1; }; mkdir -p "$OUT_BASE"; source "$CONDA_SH"; conda activate busco; command -v busco >/dev/null || { echo "busco not found" >&2; exit 1; }
cd "$OUT_BASE"; busco -i "$PROTEIN" -m proteins -l "$LINEAGE" -c "$THREADS" --offline -o "$RUN_NAME" --force
find "$OUT_BASE/$RUN_NAME" -name 'short_summary*.txt' -type f -print
