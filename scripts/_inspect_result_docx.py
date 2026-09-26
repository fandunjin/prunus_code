from pathlib import Path

from docx import Document


path = Path(__file__).resolve().parents[1] / "result.docx"
document = Document(path)
terms = (
    "method", "材料", "方法", "hifiasm", "haphic", "minimap", "quast",
    "seqkit", "busco", "inspector", "edta", "repeatmasker",
    "ltr_retriever", "lai", "helixer", "eviann", "fcs", "samtools",
    "samblaster", "bwa", "singularity", "apptainer", "pzi", "psp", "ppa",
)

print(f"paragraphs={len(document.paragraphs)} tables={len(document.tables)}")
for index, paragraph in enumerate(document.paragraphs):
    text = paragraph.text.strip()
    if text and any(term in text.lower() for term in terms):
        print(f"P{index}\t{paragraph.style.name}\t{text}")

for table_index, table in enumerate(document.tables):
    for row_index, row in enumerate(table.rows):
        text = " | ".join(cell.text.replace("\n", " / ") for cell in row.cells)
        if any(term in text.lower() for term in terms):
            print(f"T{table_index}R{row_index}\t{text}")
