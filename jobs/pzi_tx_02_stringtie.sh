#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
require_file "$GENOME_FA" "pzi32 genome FASTA"; activate_environment; require_cmd stringtie; require_cmd gffread
gtfs=()
for lib in "${LIB_IDS[@]}"; do
  bam="$ALIGN_DIR/${lib}.sorted.bam"; gtf="$STRINGTIE_DIR/${lib}.gtf"; require_file "$bam" "$lib BAM"
  if [[ -n "$GENOME_GFF" ]]; then stringtie "$bam" -p "$THREADS" -G "$GENOME_GFF" -o "$gtf" -l "$lib"; else stringtie "$bam" -p "$THREADS" -o "$gtf" -l "$lib"; fi
  gtfs+=("$gtf")
done
if [[ -n "$GENOME_GFF" ]]; then stringtie --merge -p "$THREADS" -G "$GENOME_GFF" -o "$STRINGTIE_DIR/pzi_illumina_stringtie_merged.gtf" "${gtfs[@]}"; else stringtie --merge -p "$THREADS" -o "$STRINGTIE_DIR/pzi_illumina_stringtie_merged.gtf" "${gtfs[@]}"; fi
gffread "$STRINGTIE_DIR/pzi_illumina_stringtie_merged.gtf" -g "$GENOME_FA" -w "$STRINGTIE_DIR/pzi_illumina_stringtie_merged.fa"
echo "StringTie completed: $STRINGTIE_DIR"
