#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
require_file "$GENOME_FA" "pzi32 genome FASTA"; require_file "$PACBIO_FASTQ" "PacBio Iso-Seq FASTQ"; activate_environment
for t in minimap2 samtools stringtie gffread; do require_cmd "$t"; done
pb_bam="$PACBIO_WORK_DIR/pzi_isoseq.splice.sorted.bam"; pb_gtf="$PACBIO_WORK_DIR/pzi_isoseq.stringtie.gtf"
if [[ ! -s "$pb_bam" ]]; then minimap2 -ax splice:hq -uf --secondary=no -t "$PACBIO_ALIGN_THREADS" "$GENOME_FA" "$PACBIO_FASTQ" | samtools sort -@ "$SORT_THREADS" -o "$pb_bam" -; fi
require_file "$pb_bam" "PacBio splice BAM"; [[ -s "$pb_bam.bai" ]] || samtools index -@ "$SORT_THREADS" "$pb_bam"
if [[ -n "$GENOME_GFF" ]]; then stringtie -L "$pb_bam" -p "$THREADS" -G "$GENOME_GFF" -o "$pb_gtf" -l pziIsoSeq; else stringtie -L "$pb_bam" -p "$THREADS" -o "$pb_gtf" -l pziIsoSeq; fi
gffread "$pb_gtf" -g "$GENOME_FA" -w "$PACBIO_WORK_DIR/pzi_isoseq.stringtie.fa"; samtools flagstat "$pb_bam" > "$PACBIO_WORK_DIR/pzi_isoseq.flagstat.txt"
echo "PacBio modeling completed: $pb_gtf"
