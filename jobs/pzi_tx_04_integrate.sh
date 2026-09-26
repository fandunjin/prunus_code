#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
require_file "$GENOME_FA" "pzi32 genome FASTA"; activate_environment; require_cmd stringtie; require_cmd gffread
illumina="$STRINGTIE_DIR/pzi_illumina_stringtie_merged.gtf"; pacbio="$PACBIO_WORK_DIR/pzi_isoseq.stringtie.gtf"; all="$INTEGRATED_DIR/pzi_illumina_pacbio_merged.gtf"
require_file "$illumina" "Illumina StringTie GTF"; require_file "$pacbio" "PacBio StringTie GTF"
if [[ -n "$GENOME_GFF" ]]; then stringtie --merge -p "$THREADS" -G "$GENOME_GFF" -o "$all" "$illumina" "$pacbio"; else stringtie --merge -p "$THREADS" -o "$all" "$illumina" "$pacbio"; fi
gffread "$all" -g "$GENOME_FA" -w "$INTEGRATED_DIR/pzi_illumina_pacbio_merged.fa"
if [[ "$RUN_TRINITY" == 1 ]]; then
  require_cmd Trinity
  Trinity --seqType fq --left "$CLEAN_DIR/${LIB_IDS[0]}_clean_R1.fq.gz,$CLEAN_DIR/${LIB_IDS[1]}_clean_R1.fq.gz,$CLEAN_DIR/${LIB_IDS[2]}_clean_R1.fq.gz" --right "$CLEAN_DIR/${LIB_IDS[0]}_clean_R2.fq.gz,$CLEAN_DIR/${LIB_IDS[1]}_clean_R2.fq.gz,$CLEAN_DIR/${LIB_IDS[2]}_clean_R2.fq.gz" --CPU "$THREADS" --max_memory "$MEMORY" --output "$TRINITY_DIR/trinity_out"
  cat "$TRINITY_DIR/trinity_out/Trinity.fasta" >> "$INTEGRATED_DIR/pzi_illumina_pacbio_merged.fa"
fi
cp -f "$INTEGRATED_DIR/pzi_illumina_pacbio_merged.fa" "$EVIANN_EVIDENCE_DIR/pzi_illumina_pacbio.transcripts.fa"; cp -f "$all" "$EVIANN_EVIDENCE_DIR/pzi_illumina_pacbio.transcripts.gtf"
for lib in "${LIB_IDS[@]}"; do ln -sfn "$ALIGN_DIR/${lib}.sorted.bam" "$EVIANN_EVIDENCE_DIR/${lib}.sorted.bam"; done
ln -sfn "$ALIGN_DIR/all_illumina.sorted.bam" "$EVIANN_EVIDENCE_DIR/all_illumina.sorted.bam"; ln -sfn "$PACBIO_WORK_DIR/pzi_isoseq.splice.sorted.bam" "$EVIANN_EVIDENCE_DIR/pzi_isoseq.splice.sorted.bam"
echo "Transcript evidence staged: $EVIANN_EVIDENCE_DIR"
