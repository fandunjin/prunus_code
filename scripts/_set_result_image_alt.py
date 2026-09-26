from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

from lxml import etree


DOCX = Path(__file__).resolve().parents[1] / "result.docx"
TMP = DOCX.with_suffix(".alt-text.tmp.docx")
WP_NS = "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"
ALT_TEXT = (
    "Chromosome-scale assembly and haplotype-phasing workflow for "
    "Prunus zippeliana, from ONT and Hi-C data to four phased chromosome sets."
)


with ZipFile(DOCX, "r") as source, ZipFile(TMP, "w", ZIP_DEFLATED) as target:
    for item in source.infolist():
        data = source.read(item.filename)
        if item.filename == "word/document.xml":
            root = etree.fromstring(data)
            nodes = root.xpath(
                './/wp:docPr[@id="3"]', namespaces={"wp": WP_NS}
            )
            if len(nodes) != 1:
                raise RuntimeError(f"Expected one image with wp:docPr id=3, found {len(nodes)}")
            nodes[0].set("descr", ALT_TEXT)
            data = etree.tostring(
                root, xml_declaration=True, encoding="UTF-8", standalone=True
            )
        target.writestr(item, data)

TMP.replace(DOCX)
