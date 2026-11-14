# Prunus analysis code

Code for **Racemose Prunus subgenomes reveal ancestral karyotypes and satellite
expansion after shared allopolyploidy**.

The repository contains one current implementation per analysis. Inputs and outputs
are specified through command-line arguments. Start with [commands](docs/commands.md)
and [input formats](docs/inputs.md); each executable script also accepts `--help`.

| Analysis | Main scripts |
| --- | --- |
| Haplotype assignment | `haplotype_provenance.py`, `haplotype_assignment.py`, `assign_comparative_haplotypes.py`, `partition_subgenomes.py` |
| Locus representatives and families | `primary_locus_proteomes.py`, `primary_synteny_inputs.py`, `locus_family_analysis.py`, `cafe_counts.py` |
| Phylogenomics | `alignment_qc.py`, `concatenate_alignments.py`, `contract_gene_trees.py`, `summarize_copyab_duplications.py` |
| Synteny and fractionation | `build_fractionation_catalogs.py`, `quota_raw_to_dag.py`, `parse_quota_anchors.py`, `retention_windows.py`, `fractionation_checks.py` |
| FixedAB POInT | `prepare_point_*.py`, `map_point_pillars.py`, `summarize_point_fit.py`, `compare_point_models.py`, `map_point_nodes.py` |
| Regional enrichment | `regional_go.py`, `regional_go_meta.py`, `regional_kegg.py`, `regional_kegg_meta.py`, `concordant_region_go.py` |
| Ks and branch models | `locus_pair_ks.py`, `recompute_branch_models.py`, `summarize_distinct_foregrounds.py` |
| Satellite repeats | `satellite_family_audit.py`, `satellite_unions.py`, `repeat_coverage.py`, `satellite_alignment.py`, `satellite_increment.py` |
| Chromosome source data | `pairwise_catalogs.py`, `chromosome_tracks.py`, `synteny_blocks.py`, `plot_synteny.py` |

The current analysis uses one representative per parent gene locus. Reference
results include **53 retained families**, **12,739 fractionation pillars**, and
**4,003 focal-stem duplication records**, of which **3,147** have species-overlap
support >=0.5. PziCent167 accounts for **37.3% and 40.7%** of the corresponding A1
and B1 genome-size differences.

Python 3.10+ is required. Install Python packages with
`python -m pip install -r requirements.txt`. External programs and key parameters
are listed in [docs/commands.md](docs/commands.md). The modified POInT build is
available in [point_fix](https://github.com/fandunjin/point_fix).

Run checks with `python -m unittest discover -s tests -v`.
[Validation](docs/validation.md) records the comparisons against manuscript results.
Historical figure numbers can remain in input-table names; the current mapping is
in [docs/figures.md](docs/figures.md). Raw sequences and annotation files are inputs,
not bundled software.

License: [MIT](LICENSE).
