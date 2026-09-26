#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# The TE GFF3 is checked here as a workflow guard for the downstream EviAnn QC.
ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/data/psp}"
TE_GFF3="${TE_GFF3:-${DATA_DIR}/annotation/repeats/chr32_edta/psp_2n32.chromosomes.fa.mod.EDTA.TEanno.gff3}"
INPUT="${INPUT:-${DATA_DIR}/annotation/repeats/chr32_edta/psp_2n32.chromosomes.fa.mod.MAKER.masked}"; OUT_DIR="${OUT_DIR:-${DATA_DIR}/annotation/genes/psp32_helixer_edta_masked}"; OUT_GFF="$OUT_DIR/psp_2n32.chromosomes.edta_masked.helixer.gff3"; SIF="${HELIXER_SIF:-}"; SIF_DIR="${ROOT_DIR}/software/helixer_container"; MODEL="${ROOT_DIR}/annotation/genes/helixer_models/land_plant/land_plant_v0.3_a_0080.h5"; CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
[[ -s "$INPUT" ]] || { echo "Missing EDTA-masked PSP 32 FASTA: $INPUT; run jobs/psp_32_edta.sh first" >&2; exit 1; }; [[ -s "$TE_GFF3" ]] || { echo "Missing PSP 32 EDTA TE GFF3: $TE_GFF3; run jobs/psp_32_edta.sh first" >&2; exit 1; }; [[ -n "$SIF" ]] || for f in "$SIF_DIR/helixer-docker_latest.sif" "$SIF_DIR/helixer-docker_helixer_v0.3.6_cuda_12.2.2-cudnn8.sif"; do [[ -s "$f" ]] && SIF="$f" && break; done; [[ -s "$SIF" ]] || { echo "Set HELIXER_SIF to a Helixer SIF" >&2; exit 1; }; mkdir -p "$OUT_DIR" "$ROOT_DIR/annotation/tmp/singularity"; source "$CONDA_SH"; conda activate singularity; export SINGULARITY_TMPDIR="$ROOT_DIR/annotation/tmp/singularity"
INPUT_REL="${INPUT#${ROOT_DIR}/}"; OUT_REL="${OUT_GFF#${ROOT_DIR}/}"
[[ "$INPUT_REL" != "$INPUT" && "$OUT_REL" != "$OUT_GFF" ]] || { echo "INPUT/OUT_GFF must be under ROOT_DIR for Singularity binding" >&2; exit 1; }
if [[ ! -s "$OUT_GFF" ]]; then
  args=(Helixer.py --lineage land_plant --fasta-path "/work/$INPUT_REL" --species psp_2n32_edta_masked)
  [[ -s "$MODEL" ]] && args+=(--downloaded-model-path /work/annotation/genes/helixer_models)
  args+=(--gff-output-path "/work/$OUT_REL")
  singularity run --nv -B "$ROOT_DIR:/work" "$SIF" "${args[@]}"
fi
[[ -s "$OUT_GFF" ]] || { echo "Missing Helixer GFF3: $OUT_GFF" >&2; exit 1; }; echo "PSP chr32 Helixer completed: $OUT_GFF"
