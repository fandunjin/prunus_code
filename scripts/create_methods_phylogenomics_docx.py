from pathlib import Path

from docx import Document
from docx.enum.style import WD_STYLE_TYPE
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor


OUTPUT = Path(r"${PROJECT_ROOT}\article\main_text_methods_phylogenomics_orthofinder_to_mcmctree.docx")
BLUE = RGBColor(0x2E, 0x74, 0xB5)
BLACK = RGBColor(0x00, 0x00, 0x00)


def set_cell_border(*args, **kwargs):
    # Kept for compatibility with the existing manuscript style; no tables are used here.
    return None


def set_run(run, *, italic=False, bold=False, superscript=False, size=None, color=BLACK):
    run.font.name = "Calibri"
    run._element.rPr.rFonts.set(qn("w:eastAsia"), "Calibri")
    run.font.italic = italic
    run.font.bold = bold
    run.font.superscript = superscript
    run.font.color.rgb = color
    if size is not None:
        run.font.size = Pt(size)
    return run


def add_segments(doc, segments, style="Normal"):
    paragraph = doc.add_paragraph(style=style)
    for segment in segments:
        if isinstance(segment, str):
            set_run(paragraph.add_run(segment))
            continue
        text = segment[0]
        options = segment[1] if len(segment) > 1 else {}
        set_run(paragraph.add_run(text), **options)
    return paragraph


def add_citation(paragraph, citation):
    set_run(paragraph.add_run(citation), superscript=True, size=8)


def add_body(doc, segments):
    paragraph = add_segments(doc, segments)
    paragraph.paragraph_format.space_after = Pt(6)
    paragraph.paragraph_format.line_spacing = 1.25
    return paragraph


def add_heading(doc, text, level=1):
    paragraph = doc.add_paragraph(style=f"Heading {level}")
    run = paragraph.add_run(text)
    set_run(run, bold=True, size=16 if level == 1 else 13, color=BLUE)
    paragraph.paragraph_format.space_before = Pt(14 if level == 1 else 10)
    paragraph.paragraph_format.space_after = Pt(6)
    return paragraph


def add_reference(doc, number, segments):
    paragraph = doc.add_paragraph(style="Reference")
    paragraph.paragraph_format.left_indent = Inches(0.25)
    paragraph.paragraph_format.first_line_indent = Inches(-0.25)
    paragraph.paragraph_format.space_after = Pt(3)
    paragraph.paragraph_format.line_spacing = 1.1
    set_run(paragraph.add_run(f"{number}. "), bold=False, size=9.5)
    for segment in segments:
        if isinstance(segment, str):
            set_run(paragraph.add_run(segment), size=9.5)
        else:
            text = segment[0]
            options = dict(segment[1]) if len(segment) > 1 else {}
            options.setdefault("size", 9.5)
            set_run(paragraph.add_run(text), **options)
    return paragraph


def configure_styles(doc):
    section = doc.sections[0]
    section.top_margin = Inches(1)
    section.bottom_margin = Inches(1)
    section.left_margin = Inches(1)
    section.right_margin = Inches(1)
    section.header_distance = Inches(0.492)
    section.footer_distance = Inches(0.492)

    normal = doc.styles["Normal"]
    normal.font.name = "Calibri"
    normal._element.rPr.rFonts.set(qn("w:eastAsia"), "Calibri")
    normal.font.size = Pt(11)
    normal.font.color.rgb = BLACK
    normal.paragraph_format.space_after = Pt(6)
    normal.paragraph_format.line_spacing = 1.25

    for name, size in [("Heading 1", 16), ("Heading 2", 13)]:
        style = doc.styles[name]
        style.font.name = "Calibri"
        style._element.rPr.rFonts.set(qn("w:eastAsia"), "Calibri")
        style.font.size = Pt(size)
        style.font.bold = True
        style.font.color.rgb = BLUE

    if "Reference" not in [style.name for style in doc.styles]:
        reference = doc.styles.add_style("Reference", WD_STYLE_TYPE.PARAGRAPH)
    else:
        reference = doc.styles["Reference"]
    reference.font.name = "Calibri"
    reference._element.rPr.rFonts.set(qn("w:eastAsia"), "Calibri")
    reference.font.size = Pt(9.5)
    reference.font.color.rgb = BLACK


def main():
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    doc = Document()
    configure_styles(doc)
    doc.core_properties.title = "Methods: phylogenomic inference and divergence-time estimation"
    doc.core_properties.subject = "Prunus zippeliana phylogenomics"
    doc.core_properties.author = ""

    title = doc.add_paragraph()
    title.paragraph_format.space_after = Pt(10)
    title.paragraph_format.line_spacing = 1.0
    set_run(
        title.add_run("Methods: phylogenomic inference and divergence-time estimation"),
        bold=True,
        size=16,
        color=BLUE,
    )

    add_heading(doc, "Phylogenomic sampling and subgenome representation")
    p = add_body(doc, [
        "We assembled a 17-terminal protein data set for phylogenomic inference. The focal terminals comprised the A and B subgenomes of ",
        ("Prunus zippeliana", {"italic": True}),
        " (PZI), ",
        ("P. spinulosa", {"italic": True}),
        " (PSP) and ",
        ("P. padus", {"italic": True}),
        " (PPA). For PZI and PSP, the A1 and B1 haplotype assemblies were used to represent the A and B subgenome lineages, respectively; the repartitioned PPA-A and PPA-B assemblies were used analogously. These six terminals were analysed with nine ",
        ("Prunus", {"italic": True}),
        " comparators (",
        ("P. mahaleb", {"italic": True}),
        ", ",
        ("P. avium", {"italic": True}),
        ", ",
        ("P. campanulata", {"italic": True}),
        ", ",
        ("P. speciosa", {"italic": True}),
        ", ",
        ("P. persica", {"italic": True}),
        ", ",
        ("P. dulcis", {"italic": True}),
        ", ",
        ("P. salicina", {"italic": True}),
        ", ",
        ("P. mume", {"italic": True}),
        ", and ",
        ("P. armeniaca", {"italic": True}),
        "), the proximal outgroup ",
        ("Gillenia trifoliata", {"italic": True}),
        ", and the more distant rosaceous outgroup ",
        ("Rosa chinensis", {"italic": True}),
        ". A and B were treated as subgenome-resolved lineage terminals, rather than as independent species or biological replicates. Thus, each PZI and PSP subgenome was represented once in the main analysis, avoiding pseudoreplication of homologous haplotypes within an individual.",
    ])

    add_heading(doc, "Orthogroup inference and strict single-copy selection")
    p = add_body(doc, [
        "Primary protein sets were analysed with OrthoFinder v3.1.5 using DIAMOND similarity searches (",
    ])
    add_citation(p, "1")
    set_run(p.add_run("). OrthoFinder was run with `-S diamond -t 12 -a 4`; its internal FAMSA alignments and FastTree gene trees were used only for orthology inference. We retained an orthogroup only when it contained exactly one protein from each of the 17 terminals, including PZI-A, PZI-B, PSP-A, PSP-B, PPA-A, PPA-B and "))
    set_run(p.add_run("G. trifoliata"), italic=True)
    set_run(p.add_run(". This operationally strict data set excluded incomplete and multicopy orthogroups and yielded 1,234 complete single-copy orthogroups for the main analysis."))

    add_heading(doc, "Protein alignment, filtering and gene-tree inference")
    p = add_body(doc, [
        "Each selected protein family was aligned with MAFFT v7.526 using `--localpair --maxiterate 1000 --thread 1`",
    ])
    add_citation(p, "2")
    set_run(p.add_run(" and was trimmed with trimAl v1.5.rev1 using `-automated1`"))
    add_citation(p, "3")
    set_run(p.add_run(". A trimmed alignment was retained when its length was at least 100 amino-acid sites, it contained at least 10 parsimony-informative sites, and its missing-data fraction was at most 0.20. This filtering retained 1,105 alignments; 129 families were excluded, all because they contained fewer than 10 informative sites. The retained alignments were concatenated into a 459,656-site amino-acid supermatrix with one partition per orthogroup."))

    p = add_body(doc, [
        "Maximum-likelihood gene trees were inferred independently for the 1,105 retained families with IQ-TREE v3.1.2",
    ])
    add_citation(p, "4")
    set_run(p.add_run(". Substitution models were selected per alignment with ModelFinder (`-m MFP`)"))
    add_citation(p, "5")
    set_run(p.add_run(". Node support was evaluated with 2,000 ultrafast bootstrap replicates and 1,000 SH-aLRT replicates (`-B 2000 --alrt 1000 --bnni`), using one computational thread per family and a fixed random seed of 20260825. Ultrafast bootstrap and SH-aLRT support were interpreted as complementary branch-reliability measures"))
    add_citation(p, "6,7")
    set_run(p.add_run("."))

    add_heading(doc, "Species-tree inference, branch support and targeted topology frequencies")
    p = add_body(doc, [
        "Before species-tree inference, gene-tree branches with ultrafast bootstrap support below 10% were collapsed. The 1,105 contracted, taxonomically complete gene trees were summarized with ASTRAL v5.7.8 using `-t 2`",
    ])
    add_citation(p, "8")
    set_run(p.add_run(". This primary coalescent-based analysis provided the ASTRAL topology, local posterior probabilities and quartet frequencies (q1, q2 and q3); local posterior probabilities were evaluated from quartet frequencies as implemented in ASTRAL"))
    add_citation(p, "9")
    set_run(p.add_run(". The ASTRAL topology was inferred as an unrooted tree and was rooted after inference with "))
    set_run(p.add_run("R. chinensis"), italic=True)
    set_run(p.add_run(" for display."))

    p = add_body(doc, [
        "We assessed multilocus bootstrap support with 1,000 combined gene- and site-bootstrap replicates, resampling gene-wise ultrafast-bootstrap trees (`-b`, `-g -r 1000 -s 20260825`). As a corroborative analysis, we inferred a partitioned concatenation tree from the same supermatrix with IQ-TREE (`-p concat_partitions.nex -m MFP+MERGE -B 1000 --alrt 1000 --bnni`); this concatenation analysis was not used as the primary species-tree estimate. Gene concordance factors (gCF) were calculated from the contracted gene trees, and likelihood-based site concordance factors (sCF) were calculated from 1,000 sampled sites under LG+F+G4",
    ])
    add_citation(p, "10")
    set_run(p.add_run(". Final node labels report, in order, multilocus bootstrap support, ASTRAL local posterior probability, gCF and likelihood-sCF."))

    p = add_body(doc, [
        "To quantify targeted gene-tree concordance, we scored the three possible resolutions of five prespecified quartets: PZI-A/PSP-A/PPA-A/",
        ("G. trifoliata", {"italic": True}),
        "; PZI-B/PSP-B/PPA-B/",
        ("G. trifoliata", {"italic": True}),
        "; PZI-A/PSP-A/PZI-B/PSP-B; PZI-A/PPA-A/PZI-B/PPA-B; and PSP-A/PPA-A/PSP-B/PPA-B. A gene tree was counted as decisive only if it resolved the focal quartet. Frequencies and 95% gene-resampling intervals were calculated from 1,000 resamples of the decisive contracted gene trees (seed 20260825). These topology frequencies were used as descriptive measures of gene-tree concordance; they were not used alone to infer introgression, distinguish incomplete lineage sorting from reticulation, or date allopolyploidization.",
    ])

    add_heading(doc, "Divergence-time estimation")
    p = add_body(doc, [
        "Divergence times were estimated for a fully bifurcating 17-terminal proxy tree derived from the ASTRAL topology. The analysis used the same 1,105 retained single-copy protein families and 459,656-site concatenated alignment. We first generated approximate-likelihood gradients and Hessians with IQ-TREE v3.1.2 under LG+F+G4 using `--dating mcmctree`, an independent-rates clock and the following MCMC specification: 200,000 burn-in iterations, sampling every 50 iterations, and 50,000 posterior samples. Bayesian dating was then conducted in MCMCtree, implemented in PAML v4.10.10, using the approximate-likelihood option (`usedata = 2`) and an independent-rates clock (`clock = 2`)",
    ])
    add_citation(p, "11,12")
    set_run(p.add_run(". The birth-death prior was `BDparas = 1 1 0.5 C`, with `rgene_gamma = 2 2` and `sigma2_gamma = 1 10`; the nominal chain length was 2,700,000 iterations."))

    p = add_body(doc, [
        "We used the early Eocene fossil ",
        ("Prunus wutuensis", {"italic": True}),
        " (55 Ma) as a calibration on the sampled ",
        ("Prunus", {"italic": True}),
        " crown proxy (`>0.55<1.20`), following the taxonomic interpretation of the Wutu endocarps",
    ])
    add_citation(p, "13")
    set_run(p.add_run(". A broad root condition (`RootAge = <1.2`, 120 Ma) was used as a deep-root safeguard and not as an additional fossil calibration. The 35-Ma "))
    set_run(p.add_run("P. hirsutipetala"), italic=True)
    set_run(p.add_run(" calibration was not imposed on this 17-terminal tree because the corresponding monophyletic six-terminal racemose clade was not present in the ASTRAL topology. Because MCMCtree requires a bifurcating tree and the focal A and B copies were represented as lineage terminals, posterior ages describe splits in a bifurcating proxy tree; they were not interpreted as direct ages of hybridization or allopolyploidization."))

    add_heading(doc, "References")
    references = [
        [("Emms, D. M. & Kelly, S. OrthoFinder: phylogenetic orthology inference for comparative genomics. ", {}), ("Genome Biol.", {"italic": True}), (" ", {}), ("20", {"italic": True}), (", 238 (2019). https://doi.org/10.1186/s13059-019-1832-y.", {})],
        [("Katoh, K. & Standley, D. M. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. ", {}), ("Mol. Biol. Evol.", {"italic": True}), (" ", {}), ("30", {"italic": True}), (", 772-780 (2013). https://doi.org/10.1093/molbev/mst010.", {})],
        [("Capella-Gutierrez, S., Silla-Martinez, J. M. & Gabaldon, T. trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. ", {}), ("Bioinformatics", {"italic": True}), (" ", {}), ("25", {"italic": True}), (", 1972-1973 (2009). https://doi.org/10.1093/bioinformatics/btp348.", {})],
        [("Minh, B. Q. et al. IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. ", {}), ("Mol. Biol. Evol.", {"italic": True}), (" ", {}), ("37", {"italic": True}), (", 1530-1534 (2020). https://doi.org/10.1093/molbev/msaa015.", {})],
        [("Kalyaanamoorthy, S., Minh, B. Q., Wong, T. K. F., von Haeseler, A. & Jermiin, L. S. ModelFinder: fast model selection for accurate phylogenetic estimates. ", {}), ("Nat. Methods", {"italic": True}), (" ", {}), ("14", {"italic": True}), (", 587-589 (2017). https://doi.org/10.1038/nmeth.4285.", {})],
        [("Hoang, D. T., Chernomor, O., von Haeseler, A., Minh, B. Q. & Vinh, L. S. UFBoot2: improving the ultrafast bootstrap approximation. ", {}), ("Mol. Biol. Evol.", {"italic": True}), (" ", {}), ("35", {"italic": True}), (", 518-522 (2018). https://doi.org/10.1093/molbev/msx281.", {})],
        [("Guindon, S. et al. New algorithms and methods to estimate maximum-likelihood phylogenies: assessing the performance of PhyML 3.0. ", {}), ("Syst. Biol.", {"italic": True}), (" ", {}), ("59", {"italic": True}), (", 307-321 (2010). https://doi.org/10.1093/sysbio/syq010.", {})],
        [("Zhang, C., Rabiee, M., Sayyari, E. & Mirarab, S. ASTRAL-III: polynomial time species tree reconstruction from partially resolved gene trees. ", {}), ("BMC Bioinformatics", {"italic": True}), (" ", {}), ("19", {"italic": True}), (" (Suppl. 6), 153 (2018). https://doi.org/10.1186/s12859-018-2129-y.", {})],
        [("Sayyari, E. & Mirarab, S. Fast coalescent-based computation of local branch support from quartet frequencies. ", {}), ("Mol. Biol. Evol.", {"italic": True}), (" ", {}), ("33", {"italic": True}), (", 1654-1668 (2016). https://doi.org/10.1093/molbev/msw079.", {})],
        [("Minh, B. Q., Hahn, M. W. & Lanfear, R. New methods to calculate concordance factors for phylogenomic datasets. ", {}), ("Mol. Biol. Evol.", {"italic": True}), (" ", {}), ("37", {"italic": True}), (", 2727-2733 (2020). https://doi.org/10.1093/molbev/msaa106.", {})],
        [("dos Reis, M. & Yang, Z. Approximate likelihood calculation on a phylogeny for Bayesian estimation of divergence times. ", {}), ("Mol. Biol. Evol.", {"italic": True}), (" ", {}), ("28", {"italic": True}), (", 2161-2172 (2011). https://doi.org/10.1093/molbev/msr045.", {})],
        [("Yang, Z. PAML 4: phylogenetic analysis by maximum likelihood. ", {}), ("Mol. Biol. Evol.", {"italic": True}), (" ", {}), ("24", {"italic": True}), (", 1586-1591 (2007). https://doi.org/10.1093/molbev/msm088.", {})],
        [("Li, Y. et al. Endocarps of ", {}), ("Prunus", {"italic": True}), (" (Rosaceae: Prunoideae) from the early Eocene of Wutu, Shandong Province, China. ", {}), ("Taxon", {"italic": True}), (" ", {}), ("60", {"italic": True}), (", 555-564 (2011).", {})],
    ]
    for index, reference in enumerate(references, start=1):
        add_reference(doc, index, reference)

    # Avoid an accidental blank trailing page in Word/LibreOffice.
    for paragraph in doc.paragraphs:
        paragraph.paragraph_format.widow_control = True

    doc.save(OUTPUT)
    print(OUTPUT)


if __name__ == "__main__":
    main()
