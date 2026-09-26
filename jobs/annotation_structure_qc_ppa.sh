#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
GFF3="${GFF3:-}"
for candidate in \
  "${ROOT_DIR}/data/ppa/annotation/eviann/postprocess/ppa.eviann.genes.gff3" \
  "${ROOT_DIR}/data/ppa/postprocess/ppa.eviann.genes.gff3" \
  "${ROOT_DIR}/data/ppa/annotation/genes/ppa_eviann/postprocess/ppa.eviann.genes.gff3"; do
  if [[ -z "$GFF3" && -s "$candidate" ]]; then GFF3="$candidate"; fi
done
[[ -n "$GFF3" ]] || GFF3="${ROOT_DIR}/data/ppa/annotation/eviann/postprocess/ppa.eviann.genes.gff3"
TE_GFF3="${TE_GFF3:-${PPA_EDTA_DIR}/ppa.fa.mod.EDTA.TEanno.gff3}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/data/ppa/annotation/eviann/postprocess/structure_qc}"
# PPA has no native PPA RNA-seq in this project.  Do not pass PSP/PZI BAMs as
# if they were PPA direct evidence; the transcript sequences are cross-species
# evidence used by EviAnn, not PPA splice support.
BAM_SPECS=()
source "${ROOT_DIR}/jobs/annotation_structure_qc_common.sh"
run_structure_qc
