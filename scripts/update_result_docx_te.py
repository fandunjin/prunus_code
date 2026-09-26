from pathlib import Path

from docx import Document
from docx.oxml import OxmlElement
from docx.text.paragraph import Paragraph


DOCX = Path(r"${PROJECT_ROOT}\article\result.docx")

EN_HEADING = "Transposable-element expansion and chromosome-resolved distribution"
ZH_HEADING = "转座子扩张与染色体分辨率分布"

EN_PARAGRAPHS = [
    (
        "Using the revised, independently generated EDTA annotations and the union of overlapping intervals, "
        "transposable elements (TEs) occupied 1,264.4 Mb (70.02%) of the pooled PZI chromosome set "
        "(1,805.8 Mb), compared with 472.6 Mb (44.99%) in PSP (1,050.4 Mb) and 173.1 Mb (39.12%) in PPA "
        "(442.5 Mb) (Fig. 1). Within PZI, TE fractions were 70.62%, 71.55%, 69.13% and 68.66% in hapA1, "
        "hapA2, hapB1 and hapB2, respectively. The elevated repeat fraction was therefore shared by both "
        "operational subgenomes and all four phased haplotypes, rather than being restricted to a single "
        "haplotype or chromosome set."
    ),
    (
        "Category-resolved composition showed that the PZI excess was not attributable to enrichment of one "
        "canonical TE class alone. Other TE accounted for 36.87% of PZI sequence, compared with 16.61% in PSP "
        "and 12.62% in PPA; TIR/DNA elements contributed 16.27%, 11.03% and 16.62%, respectively. Gypsy "
        "elements represented 11.77% of PZI, 13.31% of PSP and 4.64% of PPA, whereas Copia elements represented "
        "5.10%, 4.03% and 5.22%. Thus, relative to PSP, PZI contained 20.26 percentage points more Other TE "
        "and 5.24 points more TIR/DNA, while Gypsy was not enriched. In this analysis, Other TE includes "
        "unclassified LTR and other EDTA classes; these class-level estimates therefore establish the composition "
        "of the expansion but do not by themselves assign it to particular TE families."
    ),
    (
        "The 1-Mb chromosome-resolved profiles showed that high TE density was distributed across most PZI "
        "chromosomes, with local maxima varying among chromosomes and between the A and B tracks (Fig. 1). "
        "Projection of the comparative non-aligned candidate intervals onto the same C1A-C8A/C1B-C8B coordinate "
        "system identified discrete expansion windows. A subset of these intervals was supported by ONT reads "
        "spanning both boundaries and passing the depth filters, providing orthogonal evidence that part of the "
        "excess sequence represents genuine genomic sequence rather than assembly inflation. These candidate and "
        "ONT-support tracks were kept separate from the TE union-length estimates."
    ),
    (
        "Family-level attribution and insertion-age comparisons were withheld from the independent-library "
        "summary. The revised independent EDTA results are retained as a sensitivity analysis, whereas definitive "
        "cross-species family comparisons will use the frozen Prunus pan-TE library remasked onto the revised "
        "chromosome coordinates."
    ),
]

ZH_PARAGRAPHS = [
    (
        "基于重新分配后独立运行的 EDTA 注释，并对重叠区间取并集以避免嵌套注释重复计数，转座子（TE）在 PZI 四套染色体组中共占据 "
        "1,264.4 Mb（70.02%），而 PSP 和 PPA 分别占据 472.6 Mb（44.99%）和 173.1 Mb（39.12%）（图 1）。在 PZI 内部，hapA1、"
        "hapA2、hapB1 和 hapB2 的 TE 比例分别为 70.62%、71.55%、69.13% 和 68.66%。因此，较高的重复序列比例同时存在于两个操作性亚基因组及四套分相单倍型中，而不是局限于某一套单倍型或某一组染色体。"
    ),
    (
        "类别分辨率分析表明，PZI 的额外序列并非由单一经典 TE 类别的富集所解释。Other TE 在 PZI 中占 36.87%，在 PSP 和 PPA 中分别占 16.61% 和 12.62%；TIR/DNA 分别占 16.27%、11.03% 和 16.62%。Gypsy 在 PZI、PSP 和 PPA 中分别占 11.77%、13.31% 和 4.64%，Copia 分别占 5.10%、4.03% 和 5.22%。与 PSP 相比，PZI 的 Other TE 比例高 20.26 个百分点，TIR/DNA 高 5.24 个百分点，而 Gypsy 并未富集。需要指出的是，本分析中的 Other TE 包含未分类 LTR 及其他 EDTA 类别，因此这些结果可以确定扩张的类别组成，但不能单独据此将扩张归因于特定 TE 家族。"
    ),
    (
        "基于 1 Mb 窗口的染色体分辨率图谱显示，高 TE 密度分布于 PZI 的大多数染色体，并且在不同染色体以及 A、B 两条轨道之间呈现局部峰值差异（图 1）。将比较基因组中的未比对候选区间投影到统一的 C1A-C8A/C1B-C8B 坐标后，可以识别出离散的扩张窗口。其中一部分区间同时得到跨越两侧边界的 ONT reads 支持，并通过深度筛选，为部分额外序列是真实基因组序列而非组装膨胀提供了独立证据。候选区间和 ONT 支持轨道未被计入 TE 并集长度。"
    ),
    (
        "本研究未依据独立 library 汇总结果直接进行 TE 家族归因或插入年龄比较。重新生成的独立 EDTA 结果作为敏感性分析保留；最终的跨物种家族比较将采用冻结的 Prunus pan-TE library，并在重新修订的染色体坐标上完成统一 mask。"
    ),
]


def insert_paragraph_after(paragraph: Paragraph, text: str, style: str) -> Paragraph:
    new_p = OxmlElement("w:p")
    paragraph._p.addnext(new_p)
    new_para = Paragraph(new_p, paragraph._parent)
    new_para.style = style
    new_para.add_run(text)
    return new_para


def replace_paragraph(paragraph: Paragraph, text: str) -> None:
    for child in list(paragraph._p):
        if child.tag.endswith("}pPr"):
            continue
        paragraph._p.remove(child)
    paragraph.add_run(text)


def main() -> None:
    doc = Document(str(DOCX))
    headings = [p for p in doc.paragraphs if p.style.name == "Heading 2"]
    te_index, te_heading = next((i, p) for i, p in enumerate(doc.paragraphs)
                                if p.style.name == "Heading 2" and
                                ("Transposable-element" in p.text or "Transposable" in p.text))
    te_body = doc.paragraphs[te_index + 1: te_index + 5]
    if len(te_body) != 4:
        raise RuntimeError(f"Expected four English TE paragraphs, found {len(te_body)}")

    replace_paragraph(te_heading, EN_HEADING)
    for paragraph, text in zip(te_body, EN_PARAGRAPHS):
        replace_paragraph(paragraph, text)

    anchor = te_body[-1]
    zh_heading = insert_paragraph_after(anchor, ZH_HEADING, "Heading 2")
    for text in ZH_PARAGRAPHS:
        zh_heading = insert_paragraph_after(zh_heading, text, "Body Text")

    doc.save(str(DOCX))


if __name__ == "__main__":
    main()
