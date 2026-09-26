#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
require_file "$GENOME_FA" "pzi32 genome FASTA"; activate_environment
for t in fastp hisat2-build hisat2 samtools; do require_cmd "$t"; done
for lib in "${LIB_IDS[@]}"; do
  require_file "$RAW_TRANS_DIR/$lib/${lib}_1.fq.gz" "$lib R1"; require_file "$RAW_TRANS_DIR/$lib/${lib}_2.fq.gz" "$lib R2"
  fastp -i "$RAW_TRANS_DIR/$lib/${lib}_1.fq.gz" -I "$RAW_TRANS_DIR/$lib/${lib}_2.fq.gz" -o "$CLEAN_DIR/${lib}_clean_R1.fq.gz" -O "$CLEAN_DIR/${lib}_clean_R2.fq.gz" -h "$CLEAN_DIR/$lib.fastp.html" -j "$CLEAN_DIR/$lib.fastp.json" -w "$THREADS" -q 5 -n 10
done
if ! index_is_complete "$IDX_PREFIX"; then hisat2-build -p "$THREADS" "$GENOME_FA" "$IDX_PREFIX"; fi
for lib in "${LIB_IDS[@]}"; do
  hisat2 -p "$ALIGN_THREADS" -x "$IDX_PREFIX" -1 "$CLEAN_DIR/${lib}_clean_R1.fq.gz" -2 "$CLEAN_DIR/${lib}_clean_R2.fq.gz" --dta --novel-splicesite-outfile "$ALIGN_DIR/${lib}.novel_splicesites.tsv" 2> "$ALIGN_STATS_DIR/${lib}.hisat2.log" | samtools sort -@ "$SORT_THREADS" -o "$ALIGN_DIR/${lib}.sorted.bam" -
  samtools index -@ "$SORT_THREADS" "$ALIGN_DIR/${lib}.sorted.bam"; samtools flagstat "$ALIGN_DIR/${lib}.sorted.bam" > "$ALIGN_STATS_DIR/${lib}.flagstat.txt"
done
samtools merge -f -@ "$SAMTOOLS_THREADS" -o "$ALIGN_DIR/all_illumina.sorted.bam" "$ALIGN_DIR/${LIB_IDS[0]}.sorted.bam" "$ALIGN_DIR/${LIB_IDS[1]}.sorted.bam" "$ALIGN_DIR/${LIB_IDS[2]}.sorted.bam"
samtools index -@ "$SORT_THREADS" "$ALIGN_DIR/all_illumina.sorted.bam"
printf 'library_id\ttissue\tr1\tr2\n' > "$LIBRARY_METADATA"
for i in "${!LIB_IDS[@]}"; do lib="${LIB_IDS[$i]}"; printf '%s\t%s\t%s\t%s\n' "$lib" "${TISSUES[$i]}" "$RAW_TRANS_DIR/$lib/${lib}_1.fq.gz" "$RAW_TRANS_DIR/$lib/${lib}_2.fq.gz" >> "$LIBRARY_METADATA"; done
echo "Illumina mapping completed: $ALIGN_DIR"
