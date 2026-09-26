from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "result.before_methods_version_update.docx"
OUTPUT = ROOT / "result.docx"
TEMP = ROOT / "result.methods_version_update.ooxml.tmp.docx"


def replace_once(data: bytes, old: str, new: str) -> bytes:
    old_bytes = old.encode("utf-8")
    count = data.count(old_bytes)
    if count != 1:
        raise RuntimeError(f"Expected one occurrence of {old!r}; observed {count}")
    return data.replace(old_bytes, new.encode("utf-8"), 1)


if not SOURCE.is_file() or SOURCE.stat().st_size == 0:
    raise SystemExit(f"Missing source document: {SOURCE}")

with ZipFile(SOURCE) as source_zip:
    document_xml = source_zip.read("word/document.xml")

    # Normalize punctuation that rendered inconsistently in the local shell.
    artifact_counts = {
        "\u2013": document_xml.count("\u2013".encode("utf-8")),
        "\u2265": document_xml.count("\u2265".encode("utf-8")),
        "\u00d7": document_xml.count("\u00d7".encode("utf-8")),
    }
    if not all(count > 0 for count in artifact_counts.values()):
        raise RuntimeError(f"Expected punctuation was not found: {artifact_counts}")
    document_xml = document_xml.replace("\u2013".encode("utf-8"), b"-")
    document_xml = document_xml.replace("\u2265".encode("utf-8"), b">=")
    document_xml = document_xml.replace("\u00d7".encode("utf-8"), b" x ")

    replacements = [
        (
            "assembled with hifiasm in ONT mode",
            "assembled with hifiasm v0.25.0-r726 in ONT mode",
        ),
        (
            "with NCBI FCS-Adaptor and FCS-GX.",
            "with NCBI FCS-Adaptor and FCS-GX (both v0.5.5), executed as "
            "Apptainer v1.4.4 containers (MD5 "
            "27fb2938ed60835a9b41d1ce9533dbf6 and "
            "449690356545bd4ddfc6e9c6153e7de8, respectively).",
        ),
        (
            "filtered BAM files were generated with ",
            "filtered BAM files were generated with samtools v1.20 using ",
        ),
        (
            "HapHiC (v1.0.7) was then used",
            "HapHiC v1.0.7-4-g90c792c (commit "
            "90c792c45d76a499bee186273585d4242b245343) was then used",
        ),
        ("Minimap2 was run", "Minimap2 v2.30-r1287 was run"),
        (
            "with QUAST and seqkit.",
            "with QUAST v5.3.0 and seqkit v2.10.0.",
        ),
        (
            "evaluated with Inspector using",
            "evaluated with Inspector v1.3.1 using",
        ),
        (
            "annotated with EDTA v2.2.2 using ",
            "annotated with EDTA v2.2.2, including RepeatMasker v4.2.2 and "
            "LTR_retriever v3.0.4, using ",
        ),
        (
            "with panEDTA (",
            "with the panEDTA script distributed with EDTA v2.2.2 (",
        ),
        (
            "). The shared library",
            "; script MD5 63c0ec43e011b8d62729836c68ef5e83). The shared library",
        ),
        (
            "Protein-coding genes were predicted with Helixer v0.3.5 using",
            "Protein-coding genes were predicted with Helixer v0.3.5 and "
            "GeenuFF v0.3.2, executed in Singularity v3.8.6, using",
        ),
        (
            "EviAnn refinement.",
            "EviAnn v2.0.4 refinement (Bioconda package release 2.0.5).",
        ),
        (
            "converted to FASTQ with samtools and assembled",
            "converted to FASTQ with samtools v1.20 and assembled",
        ),
        (
            ", GFF3, CDS and protein packages were exported.",
            ", and an approved review status before the GFF3, CDS and protein "
            "packages were exported.",
        ),
    ]

    for old, new in replacements:
        document_xml = replace_once(document_xml, old, new)

    with ZipFile(TEMP, "w", compression=ZIP_DEFLATED) as output_zip:
        for entry in source_zip.infolist():
            payload = document_xml if entry.filename == "word/document.xml" else source_zip.read(entry.filename)
            output_zip.writestr(entry, payload)

TEMP.replace(OUTPUT)
print(f"Updated without rebuilding package: {OUTPUT}")
print(f"Punctuation replacements: {artifact_counts}")
