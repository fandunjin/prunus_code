#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
[[ -s "$GFF3" ]] || { echo "Missing GFF3: $GFF3" >&2; exit 1; }; mkdir -p "$(dirname "$OUT")"; printf 'Set\tGene number\tRNA number\tCDS number\tExon\n' > "$OUT"; awk -F'\t' 'BEGIN{g=r=c=e=0} !/^#/ {if($3=="gene")g++; else if($3=="mRNA"||$3=="transcript")r++; else if($3=="CDS")c++; else if($3=="exon")e++} END{print "psp32\t"g"\t"r"\t"c"\t"e}' "$GFF3" >> "$OUT"; cat "$OUT"
