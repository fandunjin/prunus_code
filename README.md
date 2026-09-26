# prunus_code

Analysis code for the manuscript *Subgenome-resolved Prunus zippeliana reveals one
shared whole-genome duplication in the sampled racemose Prunus lineages*.

The repository contains the command-line pipeline used for genome assembly and
phasing, structural and functional annotation, comparative genomics,
phylogenomics and divergence-time estimation, subgenome-resolved duplication
mapping (POInT, including the modified build with fixed subgenome identity),
synonymous-site and branch-model selection analyses, CAFE5 family-size turnover
and the genome-wide retention scan, together with the analysis and plotting
scripts for the figures.

## Layout

- `jobs/` - pipeline steps as submitted to a compute cluster. Each script keeps
  its stage prefix (assembly, annotation, phylogenetics, duplication mapping,
  Ks, turnover, ...). Scheduler headers are generic; adapt them to your system.
- `scripts/` - helper Python scripts used by the pipeline (input preparation,
  POInT result parsing and summaries, figure generation).
- `scripts/figure_pipeline/` - analysis and plotting scripts for the main and
  supplementary figures (scan statistics, functional enrichment, dN/dS summaries).

## Requirements

Set `PROJECT_ROOT` to the root of your data tree (and `CONDA_ROOT` to your conda
installation). Software versions are listed in the Methods of the manuscript;
the principal tools are hifiasm, HapHiC, QUAST, BUSCO, EDTA, TRASH, Helixer,
EviAnn, OrthoFinder, MAFFT, trimAl, IQ-TREE, ASTRAL, MCMCtree (PAML), WGDI,
JCVI/MCscan, KaKs_Calculator, codeml (PAML), CAFE5 and the modified POInT build.

    export PROJECT_ROOT=/path/to/project
    export CONDA_ROOT=/path/to/conda
    bash jobs/<stage_script>.sh

## Modified POInT build

The modified POInT v1.63 build that holds subgenome identity fixed at the
physical A and B chromosome columns (the `-fixedAB` option) is openly available
at https://github.com/fandunjin/point_fix.

## Data

Raw sequencing data, assemblies and annotations generated in this study are
deposited in the Genome Sequence Archive of the China National Center for
Bioinformation (BioProject PRJCAxxxxx). Public data analysed here are listed in
Supplementary Table 8 of the manuscript (GSA CRA028383; NCBI accessions).

## License

MIT (see LICENSE).
