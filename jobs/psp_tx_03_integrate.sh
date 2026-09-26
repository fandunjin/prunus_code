#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
require_file "$GENOME_FA" "PSP 32-chromosome genome"; activate_environment; require_cmd stringtie; require_cmd gffread
gtf="$STRINGTIE_DIR/psp_illumina_stringtie_merged.gtf"; require_file "$gtf" "PSP StringTie GTF"
cp -f "$gtf" "$INTEGRATED_DIR/psp_illumina_stringtie.gtf"; cp -f "$STRINGTIE_DIR/psp_illumina_stringtie_merged.fa" "$INTEGRATED_DIR/psp_illumina_stringtie.fa"
if [[ "$RUN_TRINITY" == 1 ]]; then
  require_cmd Trinity
  Trinity --seqType fq --left "$CLEAN_DIR/leaf_clean_R1.fq.gz,$CLEAN_DIR/root_clean_R1.fq.gz,$CLEAN_DIR/stem_clean_R1.fq.gz" --right "$CLEAN_DIR/leaf_clean_R2.fq.gz,$CLEAN_DIR/root_clean_R2.fq.gz,$CLEAN_DIR/stem_clean_R2.fq.gz" --CPU "$THREADS" --max_memory "$MEMORY" --output "$TRINITY_DIR/trinity_out"
  cat "$TRINITY_DIR/trinity_out/Trinity.fasta" >> "$INTEGRATED_DIR/psp_illumina_stringtie.fa"
  cp -f "$TRINITY_DIR/trinity_out/Trinity.fasta" "$EVIANN_EVIDENCE_DIR/trinity.supplementary.transcripts.fa"
fi
cp -f "$INTEGRATED_DIR/psp_illumina_stringtie.fa" "$EVIANN_EVIDENCE_DIR/psp_illumina.transcripts.fa"; cp -f "$INTEGRATED_DIR/psp_illumina_stringtie.gtf" "$EVIANN_EVIDENCE_DIR/psp_illumina.transcripts.gtf"
for lib in "${LIB_IDS[@]}"; do ln -sfn "$ALIGN_DIR/${lib}_sorted.bam" "$EVIANN_EVIDENCE_DIR/${lib}.sorted.bam"; done
ln -sfn "$ALIGN_DIR/merged_aligned.bam" "$EVIANN_EVIDENCE_DIR/merged_aligned.bam"
echo "PSP transcript evidence staged: $EVIANN_EVIDENCE_DIR"
