#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
require_file "$GENOME_FA" "PSP 32-chromosome genome"; activate_environment
for x in fastp hisat2-build hisat2 samtools; do require_cmd "$x"; done
for i in 0 1 2; do
  lib="${LIB_IDS[$i]}"; require_file "$RAW_DIR/${R1S[$i]}" "$lib R1"; require_file "$RAW_DIR/${R2S[$i]}" "$lib R2"
  fastp -i "$RAW_DIR/${R1S[$i]}" -I "$RAW_DIR/${R2S[$i]}" -o "$CLEAN_DIR/${lib}_clean_R1.fq.gz" -O "$CLEAN_DIR/${lib}_clean_R2.fq.gz" -h "$CLEAN_DIR/${lib}_fastp.html" -j "$CLEAN_DIR/${lib}_fastp.json" -w "$THREADS" -q 5 -n 10
done
if ! index_is_complete "$IDX_PREFIX"; then hisat2-build -p "$THREADS" "$GENOME_FA" "$IDX_PREFIX"; fi
for lib in "${LIB_IDS[@]}"; do
  hisat2 -p "$ALIGN_THREADS" -x "$IDX_PREFIX" -1 "$CLEAN_DIR/${lib}_clean_R1.fq.gz" -2 "$CLEAN_DIR/${lib}_clean_R2.fq.gz" --dta 2> "$ALIGN_DIR/${lib}.hisat2.log" | samtools sort -@ "$SORT_THREADS" -o "$ALIGN_DIR/${lib}_sorted.bam" -
  samtools index -@ "$SORT_THREADS" "$ALIGN_DIR/${lib}_sorted.bam"; samtools flagstat "$ALIGN_DIR/${lib}_sorted.bam" > "$ALIGN_DIR/${lib}.flagstat.txt"
done
samtools merge -f -@ "$SAMTOOLS_THREADS" "$ALIGN_DIR/merged_aligned.bam" "$ALIGN_DIR/leaf_sorted.bam" "$ALIGN_DIR/root_sorted.bam" "$ALIGN_DIR/stem_sorted.bam"; samtools index -@ "$SORT_THREADS" "$ALIGN_DIR/merged_aligned.bam"
echo "PSP Illumina mapping completed: $ALIGN_DIR"
