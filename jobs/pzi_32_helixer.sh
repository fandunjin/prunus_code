#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
# accepted as Helixer inputs. Transcript evidence is supplied later to EviAnn.
ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/annotation/chr32_edta/pzi_2n32.chromosomes.fa.mod.MAKER.masked}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/annotation/genes/pzi32_helixer_edta_masked}"
OUT_GFF="${OUT_GFF:-${OUT_DIR}/pzi_2n32.chromosomes.edta_masked.helixer.gff3}"
CONDA_SH="${CONDA_SH:-${CONDA_ROOT}/etc/profile.d/conda.sh}"
SINGULARITY_ENV="${SINGULARITY_ENV:-singularity}"
HELIXER_CONTAINER_DIR="${HELIXER_CONTAINER_DIR:-${ROOT_DIR}/software/helixer_container}"
# Keep the same image used by the previously successful Helixer jobs.
HELIXER_SIF="${HELIXER_SIF:-${ROOT_DIR}/software/helixer_container/helixer-docker_helixer_v0.3.6_cuda_12.2.2-cudnn8.sif}"
MODEL_DIR="${MODEL_DIR:-${ROOT_DIR}/annotation/genes/helixer_models}"
MODEL_FILE="${MODEL_FILE:-${MODEL_DIR}/land_plant/land_plant_v0.3_a_0080.h5}"
LINEAGE="${LINEAGE:-land_plant}"
SPECIES="${SPECIES:-pzi_2n32_edta_masked}"
STATS_SCRIPT="${ROOT_DIR}/scripts/gff3_basic_stats.py"

die(){ echo "ERROR: $*" >&2; exit 1; }
[[ -s "$INPUT_FASTA" ]] || die "masked genome missing: $INPUT_FASTA"
[[ -s "$CONDA_SH" ]] || die "conda.sh missing: $CONDA_SH"
[[ -s "$HELIXER_SIF" ]] || die "Helixer image missing; set HELIXER_SIF=/path/to/helixer.sif"
mkdir -p "$OUT_DIR" "${OUT_DIR}/stats" "${ROOT_DIR}/logs" "$MODEL_DIR" "${ROOT_DIR}/annotation/tmp/singularity"
INPUT_REL="${INPUT_FASTA#${ROOT_DIR}/}"
OUT_REL="${OUT_GFF#${ROOT_DIR}/}"
[[ "${INPUT_REL}" != "${INPUT_FASTA}" ]] || die "INPUT_FASTA must be under ROOT_DIR for Singularity binding: ${INPUT_FASTA}"
[[ "${OUT_REL}" != "${OUT_GFF}" ]] || die "OUT_GFF must be under ROOT_DIR for Singularity binding: ${OUT_GFF}"
source "$CONDA_SH"
conda activate "$SINGULARITY_ENV"
export SINGULARITY_TMPDIR="${ROOT_DIR}/annotation/tmp/singularity"

if [[ ! -s "$OUT_GFF" ]]; then
  HELIXER_CMD=(
    singularity run --nv -B "$ROOT_DIR:/work" "$HELIXER_SIF" Helixer.py
    --lineage "$LINEAGE"
    --fasta-path "/work/${INPUT_REL}"
    --species "$SPECIES"
  )
  if [[ -s "$MODEL_FILE" ]]; then
    HELIXER_CMD+=(--downloaded-model-path "/work/annotation/genes/helixer_models")
  fi
  HELIXER_CMD+=(--gff-output-path "/work/${OUT_REL}")
  "${HELIXER_CMD[@]}"
else
  echo "Existing Helixer GFF3 found, skipping: $OUT_GFF"
fi
[[ -s "$OUT_GFF" ]] || die "Helixer GFF3 missing after run: $OUT_GFF"
if command -v python3 >/dev/null 2>&1 && [[ -s "$STATS_SCRIPT" ]]; then
  python3 "$STATS_SCRIPT" --gff3 "$OUT_GFF" --prefix "$SPECIES" --out-dir "$OUT_DIR/stats"
fi
echo "Helixer chr32 completed: $OUT_GFF"
