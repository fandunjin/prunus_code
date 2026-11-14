# Key commands

Run from the repository root. Paths below are placeholders for input files and
output directories. Software versions and statistical settings follow the
manuscript Methods. Python scripts expose their full arguments with `--help`.

## Assembly and annotation

```bash
hifiasm --ont -o assembly -t 24 --h1 hic_R1.fq.gz --h2 hic_R2.fq.gz ont.fq.gz
minimap2 -x asm5 -t 16 hap1_contigs.fa curated_chromosomes.fa > hap1.paf
minimap2 -x asm5 -t 16 hap2_contigs.fa curated_chromosomes.fa > hap2.paf
EDTA.pl --genome haplotype.fa --cds peach.cds.fa --species others \
  --step all --anno 1 --sensitive 1 --force 1 --overwrite 0 --threads 24
busco -i haplotype.fa -m genome -l eudicotyledons_odb12 -o busco_genome
```

Hifiasm v0.25.0-r726, HapHiC v1.0.7, EDTA v2.2.2, Helixer v0.3.5, EviAnn v2.0.4
and BUSCO v5.8.3 were used. FCS screening and HapHiC scaffolding precede manual
Juicebox curation. Haplotype assignment uses the curated sequences, PAF files and
homolog/reference evidence:

```bash
python scripts/haplotype_provenance.py --final-fasta curated.fa \
  --hap1-paf hap1.paf --hap2-paf hap2.paf --homolog-pairs homologs.tsv \
  --out-prefix results/provenance --copy1-fasta results/copy1.fa \
  --copy2-fasta results/copy2.fa --min-mapq 10 --min-aln-len 10000 \
  --dominance-threshold 0.60
python scripts/haplotype_assignment.py --help
python scripts/assign_comparative_haplotypes.py --help
python scripts/partition_subgenomes.py --help
```

## Representative loci and gene families

```bash
python scripts/primary_locus_proteomes.py \
  --manifest inputs/primary_manifest.tsv --out results/primary
orthofinder -f results/primary/pep -S diamond -t 16 -a 4 -I 1.2 -og -o results/orthofinder
python scripts/locus_family_analysis.py --orthogroups results/orthofinder/Results_DATE/Orthogroups \
  --primary-inputs results/primary --annotations inputs/annotation/functional \
  --out results/families
python scripts/cafe_counts.py --counts results/orthofinder/Results_DATE/Orthogroups/Orthogroups.GeneCount.tsv \
  --out results/cafe
cafe5 -i results/cafe/full_counts.tsv -t inputs/seven_species.tre -k 3 -c 4 -o results/cafe/full
cafe5 -i results/cafe/diploid_counts.tsv -t inputs/four_diploids.tre -k 3 -c 4 -o results/cafe/diploid
```

Use completed OrthoFinder 3.1.5 results, not initial MCL clusters. For duplication
mapping, concatenate A/B representative proteins within each focal species into
14 species-level FASTA files and run OrthoFinder without `-og`:

```bash
orthofinder -f inputs/species14 -S diamond -t 8 -a 4 -I 1.2 -o results/species14
python scripts/summarize_copyab_duplications.py --result results/species14/Results_DATE \
  --outdir results/duplications --support 0.5
```

## Phylogenomics and dating

For each single-copy orthogroup:

```bash
mafft --localpair --maxiterate 1000 gene.fa > gene.aligned.fa
trimal -in gene.aligned.fa -out gene.trimmed.fa -automated1
```

Filter alignments, then apply the IQ-TREE command to each retained family:

```bash
python scripts/alignment_qc.py --alignment-dir inputs/trimmed \
  --species-file inputs/taxa.txt --out-dir results/alignment_qc \
  --min-length 100 --min-informative 10 --max-missing 0.20
python scripts/concatenate_alignments.py --alignment-dir results/alignment_qc/retained \
  --out-fasta results/concatenated.fa --out-partitions results/partitions.nex
iqtree3 -s gene.trimmed.fa -m MFP -B 2000 --alrt 1000 --bnni \
  -T 1 -seed 20260825 --prefix results/gene
python scripts/contract_gene_trees.py --tree-dir results/gene_trees \
  --threshold 10 --out results/contracted.tre
java -jar astral.jar -i results/contracted.tre -o results/astral.tre -t 2
java -jar astral.jar -i results/contracted.tre -b inputs/bootstrap_files.txt \
  -g -r 1000 -s 20260825 -o results/astral_bootstrap.tre -t 2
iqtree3 -s results/concatenated.fa -te inputs/dated_proxy.tre \
  -m LG+F+G4 --dating mcmctree --prefix results/dating
mcmctree inputs/mcmctree.ctl
```

MCMCtree uses `usedata=2`, `clock=2`, `ndata=1`, `BDparas=1 1 0.5 C`,
`rgene_gamma=2 2`, `sigma2_gamma=1 10`, burn-in 200,000, sampling frequency 50,
and 50,000 samples. The Prunus crown calibration is `>0.55<1.20`; the soft root
constraint is `<1.2`, in 100-Ma units. Original proxy trees, calibrated nodes and
alignments accompany the source data.

## Synteny and fractionation

```bash
python scripts/primary_synteny_inputs.py --catalog inputs/chromosome_catalog \
  --canonical inputs/canonical_annotations --primary-inputs results/primary \
  --out results/synteny_inputs
python scripts/build_fractionation_catalogs.py --inputs-root results/synteny_inputs \
  --outdir results/fractionation/01_catalog
lastdb -p -c results/GTR_pep results/fractionation/01_catalog/tracks/GTR.pep
lastal -p BL62 -P 2 -i256M -f BlastTab results/GTR_pep query_track.pep > query_track.last.tsv
```

Apply LAST separately to the eight query tracks and concatenate their outputs.
Quota Align preparation uses `blast_to_raw.py --no_strip_names --tandem_Nmax=10`,
`quota_raw_to_dag.py`, DAGchainer `-D 20 -g 10 -A 5`, and
`cluster_utils.py --format=dag --log_evalue`. Merge and depth-filter the blocks:

```bash
python quota_align.py --Dm=20 --merge inputs/initial.qa
python quota_align.py --Nm=40 --quota=8:1 --solver=GLPK inputs/initial.qa.merged
python scripts/parse_quota_anchors.py --quota-output inputs/initial.qa.merged.filtered \
  --catalog results/fractionation/01_catalog --query-bed inputs/query8.nolocaldups.bed \
  --target-bed inputs/GTR.nolocaldups.bed --outdir results/fractionation/04_quota_align/normalized
python scripts/retention_windows.py --anchors results/fractionation/04_quota_align/normalized/quota_filtered_anchors.tsv \
  --catalog results/fractionation/01_catalog --window-size 100 --outdir results/fractionation/05_fractbias
```

The logical query:GTR depth is 8:1; use `--quota=1:8` when the file columns put GTR
first. Preserve the retained BED order and merged QA files for regional analysis.
WGDI uses `wgdi -d`, `-icl`, `-bi`, `-c` with pair-specific configuration files for
ancestral-karyotype comparisons; accepted ancestral units are manually curated
from the resulting conserved blocks and reciprocal intervals.

## FixedAB POInT

```bash
python scripts/prepare_point_pillars.py --quota-anchors inputs/quota_filtered_anchors.tsv \
  --catalog results/fractionation/01_catalog --outdir results/point/pillars
python scripts/prepare_point_inputs.py --pillars results/point/pillars/gtr_centered_pillars.tsv \
  --catalog results/fractionation/01_catalog --outdir results/point/inputs
POInT_fixedAB -d:2 -g:PZI_POInT_geneorder.txt -g:PSP_POInT_geneorder.txt \
  -g:PPA_POInT_geneorder.txt -o:point_pillars_optimized.txt \
  -m:WGD_bias_fix_fixedAB_model.txt -t:species_tree.nex \
  -p:posterior_probs.txt -z -l:1e-7 -fixedAB
python scripts/prepare_point_pairs.py --inputs results/point/inputs \
  --optimized inputs/point_pillars_optimized.txt --outdir results/point/pairs
python scripts/prepare_point_singletons.py --inputs results/point/inputs \
  --optimized inputs/point_pillars_optimized.txt --outdir results/point/singletons
python scripts/compare_point_models.py --shared-fit inputs/point_model_fit.tsv \
  --pair-root results/point/pair_fits --independent-root results/point/singleton_fits \
  --pillars 12739 --optimized-pillars inputs/point_pillars_optimized.txt \
  --process-parameters 2 --out results/point/comparison.tsv
```

Run the corresponding POInT fits on each prepared pair/singleton input.
`summarize_point_fit.py` parses native outputs; `map_point_nodes.py --fixed-ab`
computes node retention/loss summaries. Singleton log likelihood is verified
against the exact three-state multinomial likelihood; model parameter counts are
7/7/6. `fractionation_checks.py` audits locus uniqueness, node conservation and
20,000 paired chromosome-bootstrap resamples.

## Regional enrichment and branch models

```bash
python scripts/pillar_retention.py --pillars inputs/gtr_centered_pillars.tsv --out analysis/figures
python scripts/regional_windows.py --source analysis/figures/figure4_panel_a_source_data.tsv \
  --density analysis/figures/figure4_panel_a_density_source.tsv --out analysis/GO/window_retention_profile.tsv
python scripts/regional_go.py --root analysis
python scripts/regional_go_meta.py --root analysis
python scripts/concordant_region_go.py --root analysis
python scripts/kegg_annotations.py --root analysis
python scripts/regional_kegg.py --root analysis
python scripts/regional_kegg_meta.py --root analysis
python scripts/regional_gene_summary.py --root analysis
python scripts/locus_pair_ks.py --canonical inputs/canonical_annotations \
  --pairs inputs/ks_plot_source.tsv --out results/ks
python scripts/recompute_branch_models.py --source inputs/branch_models \
  --out results/dnds --codeml /software/codeml --workers 8
python scripts/summarize_distinct_foregrounds.py --run results/dnds
```

GO/KEGG tests use the full eligible term universe before BH correction.
Cross-species support uses `min(1, 2 * second_smallest_p)` with missing species P=1.
Codeml uses both recorded starting omega values and checks nested likelihoods;
these finite optimization runs are part of the analysis. The 16 original
foreground labels are summarized as 11 distinct unrooted assignments.

## Satellite repeats

```bash
python scripts/satellite_family_audit.py --source inputs/satellites --out results/satellites
python scripts/satellite_unions.py --source inputs/satellites --census-dir results/satellites \
  --genomes inputs/genomes.tsv --centromeres inputs/centromere_proxies --out results/satellite_unions
python scripts/repeat_coverage.py --manifest inputs/repeat_targets.tsv \
  --census-dir results/satellites --out results/repeat_coverage
python scripts/satellite_alignment.py --source inputs/satellites \
  --census-dir results/satellites --out results/satellite_alignment
python scripts/satellite_increment.py --composition results/repeat_coverage/exclusive_composition.tsv \
  --summary results/repeat_coverage/exclusive_summary.tsv --out results/genome_increment
```

Family membership uses BLASTN identity >=70%, >=100 aligned columns and query
bases, E <=1e-5, and 163–170-bp monomers. The coverage sweep assigns tandem
monomers first, followed by EDTA score/span priority. Genome denominators use full
FASTA lengths. The five-species alignment uses observed representatives, circular
phase/strand normalization and MAFFT L-INS-i. SF13D contains the alignment and
conservation legend, with its former top/bottom explanatory lines removed.

## Corrected chromosome source tables

```bash
python scripts/pairwise_catalogs.py --catalog results/fractionation/01_catalog/tracks --out results/fig2
python scripts/chromosome_tracks.py --project-root data/project --out results/fig2
```

For each focal pair, run DIAMOND `blastp --sensitive --evalue 1e-5 --max-target-seqs 10`,
JCVI `blastfilter --cscore=.70`, and synteny scanning with minimum size 4 and
chaining distance 30. Place the resulting `pzi.psp.anchors`, `pzi.ppa.anchors` and
`psp.ppa.anchors` beside the pooled BED files:

```bash
python scripts/synteny_blocks.py --project-root data/project --out results/fig2
python scripts/plot_synteny.py --project-root data/project --out results/fig2
```

The tables contain 29 selected blocks per focal pair and 312 blocks across the
six adjacent-species comparisons. The final JCVI step updates display coordinates
in `riparian_layout.tsv`. These are source-table calculations; submitted main-figure
artwork is handled separately as described in [figures.md](figures.md).
