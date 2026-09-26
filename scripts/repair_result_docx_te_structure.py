from pathlib import Path

from docx import Document
from docx.oxml import OxmlElement
from docx.text.paragraph import Paragraph


DOCX = Path(r"${PROJECT_ROOT}\article\result.docx")


def replace_paragraph(paragraph: Paragraph, text: str, style: str | None = None) -> None:
    if style:
        paragraph.style = style
    for child in list(paragraph._p):
        if child.tag.endswith("}pPr"):
            continue
        paragraph._p.remove(child)
    paragraph.add_run(text)


def insert_before(paragraph: Paragraph, text: str, style: str) -> Paragraph:
    new_p = OxmlElement("w:p")
    paragraph._p.addprevious(new_p)
    new_para = Paragraph(new_p, paragraph._parent)
    new_para.style = style
    new_para.add_run(text)
    return new_para


def main() -> None:
    doc = Document(str(DOCX))
    # Current layout after the first insertion:
    # 15 TE heading, 16-18 English TE body, 19 English caveat carrying the
    # accidental Heading 2 style, 20-24 Chinese TE section, 25-27 karyotype body.
    if len(doc.paragraphs) < 28:
        raise RuntimeError("Unexpected paragraph count; refusing structural repair")

    doc.paragraphs[19].style = "Body Text"
    doc.paragraphs[25].style = "Body Text"
    insert_before(doc.paragraphs[25], "An inferred nine-unit ancestral karyotype", "Heading 2")
    doc.save(str(DOCX))


if __name__ == "__main__":
    main()
