#!/bin/bash
# Adapt the scheduler header to your cluster (SLURM, PBS, JSUB, ...).
MAP="${MAP:-${DATA_DIR}/haploid/psp_A1A2B1B2.rename_map.tsv}"; IN_DIR="${IN_DIR:-${DATA_DIR}/annotation/genes/psp32_eviann/postprocess}"; PREFIX="${PREFIX:-psp32.eviann}"; OUT_DIR="${OUT_DIR:-${DATA_DIR}/annotation/genes/psp32_haplotypes}"
GFF="$IN_DIR/${PREFIX}.genes.gff3"; GENOME="${GENOME:-${DATA_DIR}/haphic/merge/04.build/final_assembly_qc_outputs/psp_2n32.chromosomes.fa}"
[[ -s "$MAP" && -s "$GFF" && -s "$GENOME" ]] || { echo "Missing map, PSP32 genome, or full PSP32 GFF3" >&2; exit 1; }; command -v gffread >/dev/null || { echo "gffread is required" >&2; exit 1; }; mkdir -p "$OUT_DIR"
declare -A H; while IFS=$'\t' read -r old new hap rest; do [[ "$old" == old_seqid || -z "$hap" ]] && continue; H["$old"]="$hap"; done < "$MAP"
for hap in hapA1 hapA2 hapB1 hapB2; do
  d="$OUT_DIR/$hap"; mkdir -p "$d"; out="$d/${PREFIX}.${hap}.genes.gff3"; : > "$out"
  for chr in "${!H[@]}"; do [[ "${H[$chr]}" == "$hap" ]] && awk -v c="$chr" 'BEGIN{OFS="\t"} /^#/||$1==c{print}' "$GFF" >> "$out"; done
  gffread "$out" -g "$GENOME" -x "$d/${PREFIX}.${hap}.cds.fa" -y "$d/${PREFIX}.${hap}.protein.fa" -w "$d/${PREFIX}.${hap}.mrna.fa" || { echo "gffread extraction failed for $hap" >&2; exit 1; }
done
printf 'Set\tGene number\tRNA number\tCDS number\tExon\n' > "$OUT_DIR/psp32_annotation_summary.tsv"
for item in psp32 hapA1 hapA2 hapB1 hapB2; do f="$GFF"; [[ "$item" != psp32 ]] && f="$OUT_DIR/$item/${PREFIX}.${item}.genes.gff3"; awk -F'\t' -v n="$item" 'BEGIN{g=r=c=e=0}!/^#/{if($3=="gene")g++;else if($3=="mRNA"||$3=="transcript")r++;else if($3=="CDS")c++;else if($3=="exon")e++}END{print n"\t"g"\t"r"\t"c"\t"e}' "$f" >> "$OUT_DIR/psp32_annotation_summary.tsv"; done
cat "$OUT_DIR/psp32_annotation_summary.tsv"
