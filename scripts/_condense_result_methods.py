from pathlib import Path
import shutil

from docx import Document


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "result.docx"
BACKUP = ROOT / "result.before_methods_condense.docx"
TEMP = ROOT / "result.methods_condense.tmp.docx"


REPLACEMENTS = {
    6: (
        "The P. zippeliana assembly was generated from sample HD-14-2DNA using Oxford Nanopore Technologies (ONT) long reads and paired-end Hi-C reads."
    ),
    8: (
        "ONT and Hi-C reads were assembled with hifiasm v0.25.0-r726 in ONT mode (--ont)."
    ),
    9: (
        "Primary, hap1 and hap2 contigs were screened with FCS-Adaptor and FCS-GX v0.5.5 (Prunus taxon identifier 3754) in Apptainer v1.4.4. Hi-C reads were processed with BWA-MEM v0.7.18, samblaster v0.1.26 and samtools v1.20, and HapHiC v1.0.7-4-g90c792c produced 16- and 32-group assemblies. Contact maps were curated in Juicebox."
    ),
    11: (
        "The 32 chromosomes were aligned to FCS-cleaned hifiasm hap1 and hap2 contigs with minimap2 v2.30-r1287 (asm5). Alignments with MAPQ >= 10 and length >= 10 kb were merged, and a support-dominance threshold of 0.60 assigned two 16-chromosome phase sets (copy1 and copy2)."
    ),
    12: (
        "Four operational x = 8 haplotype groups were inferred by integrating self-synteny, hifiasm provenance, HapHiC-supported homolog pairs, chromosome-scale Prunus and Rosaceae references, and 500-kb segment support. Reference alignments used minimap2 asm10 with the same MAPQ and length thresholds."
    ),
    13: (
        "Quartet scores combined copy balance (0.16), corrected homolog support (0.14), P. padus pairing (0.14), x = 8 reference support (0.20), outgroup support (0.10), segment support (0.16), self-synteny (0.06) and length balance (0.04). Exact-cover optimization selected eight disjoint A1/A2/B1/B2 quartets. These are operational phase labels, not assignments to parental subgenomes."
    ),
    15: (
        "Assembly statistics were calculated with QUAST v5.3.0 and seqkit v2.10.0, and completeness with BUSCO v5.8.3 using eudicotyledons_odb12 (2,805 orthologues). Inspector v1.3.1 assessed ONT-read consistency; its switch and structural-error estimates were interpreted as cross-haplotype mapping diagnostics rather than allele-specific truth."
    ),
    17: (
        "Transposable elements were annotated per haplotype with EDTA v2.2.2 (--species others --step all --anno 1 --sensitive 1), using P. persica CDS (GCF_000346465.2) to limit coding-sequence misclassification. EDTA used RepeatMasker v4.2.2 and LTR_retriever v3.0.4; libraries were integrated with panEDTA (-f 3 -a 1 -o 0) for uniform reannotation."
    ),
    18: (
        "LTR Assembly Index values were calculated with LAI beta3.2 from EDTA-derived intact- and total-LTR annotations."
    ),
    19: (
        "Genes were predicted with Helixer v0.3.5 (land-plant model; GeenuFF v0.3.2 runtime) and refined with EviAnn v2.0.4 using CDS and protein evidence from P. persica, P. campanulata, P. mume and P. speciosa, plus Swiss-Prot. Models with CDS < 300 bp were removed; when EDTA annotations were supplied, gene-to-TE and CDS-to-TE overlap thresholds were 0.80 and 0.50. Primary isoforms were evaluated with BUSCO."
    ),
    21: (
        "Public P. spinulosa PacBio HiFi, Hi-C, Illumina and RNA-seq data were obtained from GSA BioProject CRA028383 (CRR2018654-CRR2018659). HiFi reads were assembled with hifiasm v0.25.0-r726 using Hi-C phasing. FCS-screened contigs (taxon identifier 170480) were scaffolded with the HapHiC workflow above into 16- and 32-group assemblies and reviewed in Juicebox."
    ),
    22: (
        "The 32 curated chromosomes were assigned 16:16 to hifiasm hap1/hap2 from minimap2 asm5 alignments (MAPQ >= 10; length >= 10 kb), then matched to the standardized P. zippeliana C01-C08_A1/A2/B1/B2 chromosomes using asm10. Identity-weighted, non-overlapping coverage defined a maximum-weight one-to-one assignment within each phase; strand-dominant chromosomes were reverse-complemented. The four x = 8 sets were annotated with EDTA, LAI and Helixer as above."
    ),
    24: (
        "The P. padus genome was obtained from NCBI GenBank (GCA_964059315.1; BioProject PRJEB75998; BioSample SAMEA110449347); chromosomes 1-16 were retained for correspondence analysis."
    ),
    25: (
        "P. padus genes were predicted with Helixer v0.3.5 and refined with EviAnn v2.0.4 using P. spinulosa copy1 CDS and protein evidence; CDS < 300 bp were removed."
    ),
    26: (
        "P. padus was aligned to P. zippeliana copy1 with minimap2 asm5. Primary alignments with MAPQ >= 20 and length >= 10 kb were merged, and greatest non-overlapping coverage defined a one-to-one correspondence among the 16 chromosomes."
    ),
    27: (
        "P. zippeliana C##_A1 and C##_B1 correspondence defined P. padus C##_A and C##_B labels after manual review; unplaced sequences were excluded."
    ),
    28: (
        "These A/B labels indicate structural correspondence, not parental-donor identity or a parentally phased P. padus assembly."
    ),
}


EXPECTED_HEADINGS = {
    4: "Methods",
    5: "Plant material and sequencing of Prunus zippeliana",
    7: "De novo assembly and chromosome-scale scaffolding of P. zippeliana",
    10: "Chromosome provenance and operational phase assignment in P. zippeliana",
    14: "Assembly and annotation quality assessment",
    16: "Repeat annotation, LTR quality and protein-coding gene prediction",
    20: "Public-data reassembly and P. zippeliana-constrained phasing of Prunus spinulosa",
    23: "P. zippeliana-guided chromosome correspondence and A/B partitioning of Prunus padus",
}


def replace_text_preserving_first_run(paragraph, text):
    if paragraph.runs:
        paragraph.runs[0].text = text
        for run in paragraph.runs[1:]:
            run.text = ""
    else:
        paragraph.add_run(text)


if not SOURCE.is_file() or SOURCE.stat().st_size == 0:
    raise SystemExit(f"Missing input document: {SOURCE}")

document = Document(SOURCE)
if len(document.paragraphs) != 31:
    raise SystemExit(f"Unexpected paragraph count: {len(document.paragraphs)}")
for index, expected in EXPECTED_HEADINGS.items():
    observed = document.paragraphs[index].text.strip()
    if observed != expected:
        raise SystemExit(f"Unexpected paragraph {index}: {observed!r}")

if not BACKUP.exists():
    shutil.copy2(SOURCE, BACKUP)

before_words = sum(len(document.paragraphs[i].text.split()) for i in range(5, 29))
for index, text in REPLACEMENTS.items():
    replace_text_preserving_first_run(document.paragraphs[index], text)
after_words = sum(len(document.paragraphs[i].text.split()) for i in range(5, 29))

document.save(TEMP)
if not TEMP.is_file() or TEMP.stat().st_size == 0:
    raise SystemExit("Updated document was not written")
TEMP.replace(SOURCE)

print(f"Updated: {SOURCE}")
print(f"Backup:  {BACKUP}")
print(f"Methods words: {before_words} -> {after_words}")
