# Inputs

## Protein and synteny catalogs

`primary_locus_proteomes.py` reads a TSV with columns `track`, `frozen_pep`,
`id_map`, `restored_pep`, `cds`; see `examples/primary_manifest.tsv`.
The identifier map contains `stable_id`, `original_id`, `original_header`,
`protein_length`. Original headers identify the parent gene explicitly through
`gene=`, `locus=`, `gene-...|`, or the G. trifoliata transcript naming convention.
CDS lengths come from `cds_len=` or the supplied CDS FASTA. The longest valid CDS
wins, followed by protein length and transcript identifier; minimum CDS is 300 nt.
Identical proteins from different loci remain separate records.

Family tracks are `pzi_A`, `pzi_B`, `psp_A`, `psp_B`, `ppa_A`, `ppa_B`, `pcam`,
`pav`, `pspe`, `pmal`, `par`, `pmu`, `pdu`, `psal`, `ppe`, `gtr`, `rch`.
`primary_synteny_inputs.py` reads the existing chromosome-coordinate BED/protein/CDS
catalog, original comparator GFF3 files, and representative ID maps. It produces
locus-unique coordinate catalogs for the fractionation pipeline.

## Regional analysis root

The regional GO/KEGG commands accept `--root`, pointing to a directory with:

```
annotation/functional/{pzi,psp,ppa}/swissprot_*.tsv
figures/figure4_panel_a_source_data.tsv
figures/figure4_panel_a_density_source.tsv
GO/
downstream_analysis/figure2_polyploid/strict_wgd_mapping/
  04_fractbias_8track_reference_exact/
    01_catalog/gtr_gene_order.tsv
    02_last/{GTR,query8}.nolocaldups.bed
    04_quota_align/query8_vs_GTR.initial.qa.merged
    04_quota_align/normalized/quota_filtered_anchors.tsv
```

These are analysis-relative filenames retained for compatibility with the source
tables. Functional annotations are restricted to the selected representative loci.
GO rows contain `protein_id`, `go_id`, `go_aspect`, `go_description`, `go_evidence`.
KEGG mapping uses local REST response files in `annotation/functional/kegg_cache/`;
the script prints the required public URL if a file is missing.

## Satellite inputs

The satellite source directory contains:

```
analysis/published_monomers.fa
analysis/published_monomers_wang2022.fa
trash/<target>/<target>.fa_repeats_with_seq.csv
trash/<target>/<target>.fa_arrays.csv
```

TRASH targets are `pziA1`, `pziA2`, `pziB1`, `pziB2`, `pspA1`, `pspA2`, `pspB1`,
`pspB2`, `ppaA`, `ppaB`, `pcam`, `ppe`, `gtr`, `rch`. Original coordinates are
one-based inclusive; generated interval tables are zero-based half-open.

`repeat_coverage.py` reads `target`, `fasta`, `edta_gff`, `lineage_table` from a TSV.
See `examples/repeat_targets.tsv`. Its eight labels are `PZI-A1`, `PZI-B1`,
`PSP-A1`, `PSP-B1`, `PPA-A`, `PPA-B`, `PCAM`, `PPE`. Each FASTA has a `.fai` index.
REXdb lineage tables have `record_id`, `class`, `order`, `family`, `sublineage`,
`clade` columns. They are fixed inputs to the exclusive-coverage calculation.

`satellite_unions.py` additionally reads a TSV of `target` and `fai`, and a
predicted-centromere directory containing `<target>/pericentromere_regions.tsv`
with `chrom`, `start`, `end` columns. Genome lengths include ambiguous bases.

## Chromosome source-data root

`chromosome_tracks.py`, `synteny_blocks.py` and `plot_synteny.py` use `--project-root`
containing `0903/haploid_repartition`, `0903/annotation` and
`0903/native_multispecies_A1B1_pcam_pairwise`. The verified haplotype maps, current
FASTA/GFF/EDTA files and archived adjacent-species anchors are the required inputs.
`pairwise_catalogs.py` pools the final `PZI_A/PZI_B`, `PSP_A/PSP_B` and `PPA_A/PPA_B`
BED/protein catalogs without renumbering gene identifiers.
