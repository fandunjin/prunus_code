from pathlib import Path
from docx import Document

DOCX = Path(r"${PROJECT_ROOT}\article\result.docx")
ZH_HEADING = "转座子扩张与染色体分辨率分布"
ZH_PARAGRAPHS = [
    "基于重新分配后独立运行的 EDTA 注释，并对重叠区间取并集以避免嵌套注释重复计数，转座子（TE）在 PZI 四套染色体组中共占据 1,264.4 Mb（70.02%），而 PSP 和 PPA 分别占据 472.6 Mb（44.99%）和 173.1 Mb（39.12%）（图 1）。在 PZI 内部，hapA1、hapA2、hapB1 和 hapB2 的 TE 比例分别为 70.62%、71.55%、69.13% 和 68.66%。因此，较高的重复序列比例同时存在于两个操作性亚基因组及四套分相单倍型中，而不是局限于某一套单倍型或某一组染色体。",
    "类别分辨率分析表明，PZI 的额外序列并非由单一经典 TE 类别的富集所解释。Other TE 在 PZI 中占 36.87%，在 PSP 和 PPA 中分别占 16.61% 和 12.62%；TIR/DNA 分别占 16.27%、11.03% 和 16.62%。Gypsy 在 PZI、PSP 和 PPA 中分别占 11.77%、13.31% 和 4.64%，Copia 分别占 5.10%、4.03% 和 5.22%。与 PSP 相比，PZI 的 Other TE 比例高 20.26 个百分点，TIR/DNA 高 5.24 个百分点，而 Gypsy 并未富集。需要指出的是，本分析中的 Other TE 包含未分类 LTR 及其他 EDTA 类别，因此这些结果可以确定扩张的类别组成，但不能单独据此将扩张归因于特定 TE 家族。",
    "基于 1 Mb 窗口的染色体分辨率图谱显示，高 TE 密度分布于 PZI 的大多数染色体，并且在不同染色体以及 A、B 两条轨道之间呈现局部峰值差异（图 1）。将比较基因组中的未比对候选区间投影到统一的 C1A-C8A/C1B-C8B 坐标后，可以识别出离散的扩张窗口。其中一部分区间同时得到跨越两侧边界的 ONT reads 支持，并通过深度筛选，为部分额外序列是真实基因组序列而非组装膨胀提供了独立证据。候选区间和 ONT 支持轨道未被计入 TE 并集长度。",
    "本研究未依据独立 library 汇总结果直接进行 TE 家族归因或插入年龄比较。重新生成的独立 EDTA 结果作为敏感性分析保留；最终的跨物种家族比较将采用冻结的 Prunus pan-TE library，并在重新修订的染色体坐标上完成统一 mask。",
]

def set_east_asia(run, font="Noto Sans SC"):
    run.font.name = font
    rpr = run._element.get_or_add_rPr()
    rfonts = rpr.rFonts
    if rfonts is None:
        from docx.oxml import OxmlElement
        rfonts = OxmlElement("w:rFonts")
        rpr.insert(0, rfonts)
    rfonts.set("{http://schemas.openxmlformats.org/wordprocessingml/2006/main}eastAsia", font)
    rfonts.set("{http://schemas.openxmlformats.org/wordprocessingml/2006/main}hAnsi", font)
    rfonts.set("{http://schemas.openxmlformats.org/wordprocessingml/2006/main}ascii", font)

def replace(p, text, style=None, font=None, bold=None):
    if style:
        p.style = style
    for child in list(p._p):
        if child.tag.endswith("}pPr"):
            continue
        p._p.remove(child)
    run = p.add_run(text)
    if bold is not None:
        run.bold = bold
    if font:
        set_east_asia(run, font)

doc = Document(str(DOCX))
paras = doc.paragraphs
# Replace the corrupted Chinese heading and body paragraphs in-place.
replace(paras[20], ZH_HEADING, style="Heading 2", font="Noto Sans SC", bold=True)
for p, txt in zip(paras[21:25], ZH_PARAGRAPHS):
    replace(p, txt, style="Body Text", font="Noto Serif SC")
doc.save(str(DOCX))
