#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
require_file "$GENOME_FA" "PSP 32-chromosome genome"; activate_environment; require_cmd stringtie; require_cmd gffread
gtfs=()
for lib in "${LIB_IDS[@]}"; do
  gtf="$STRINGTIE_DIR/${lib}.gtf"; stringtie "$ALIGN_DIR/${lib}_sorted.bam" -p "$THREADS" -o "$gtf" -l "psp_${lib}"; gtfs+=("$gtf")
done
stringtie --merge -p "$THREADS" -o "$STRINGTIE_DIR/psp_illumina_stringtie_merged.gtf" "${gtfs[@]}"
gffread "$STRINGTIE_DIR/psp_illumina_stringtie_merged.gtf" -g "$GENOME_FA" -w "$STRINGTIE_DIR/psp_illumina_stringtie_merged.fa"
echo "PSP StringTie completed: $STRINGTIE_DIR"
