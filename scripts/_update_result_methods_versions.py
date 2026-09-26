from pathlib import Path
import shutil

from docx import Document


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "result.docx"
BACKUP = ROOT / "result.before_methods_version_update.docx"
TEMP = ROOT / "result.methods_version_update.tmp.docx"


def replace_in_runs(paragraph, old, new, expected=1):
    count = 0
    for run in paragraph.runs:
        if old in run.text:
            occurrences = run.text.count(old)
            run.text = run.text.replace(old, new)
            count += occurrences
    if count != expected:
        raise RuntimeError(
            f"Expected {expected} replacement(s) in paragraph {paragraph.text!r}: "
            f"{old!r}; observed {count}"
        )


if not SOURCE.is_file() or SOURCE.stat().st_size == 0:
    raise SystemExit(f"Missing input document: {SOURCE}")

if not BACKUP.exists():
    shutil.copy2(SOURCE, BACKUP)

document = Document(SOURCE)
paragraphs = document.paragraphs
if len(paragraphs) < 29 or paragraphs[4].text.strip() != "Methods":
    raise SystemExit("Unexpected result.docx structure; Methods section was not found")

# Repair encoding artifacts already present in the document while preserving runs.
artifact_counts = {"\u2013": 0, "\u2265": 0, "\u00d7": 0}
for paragraph in paragraphs:
    for run in paragraph.runs:
        for token in artifact_counts:
            artifact_counts[token] += run.text.count(token)
        run.text = run.text.replace("\u2013", "-")
        run.text = run.text.replace("\u2265", ">=")
        run.text = run.text.replace("\u00d7", " x ")

if not all(count > 0 for count in artifact_counts.values()):
    raise RuntimeError(f"Expected encoding artifacts were not all present: {artifact_counts}")

# PZI assembly and scaffolding software versions.
replace_in_runs(
    paragraphs[8],
    "assembled with hifiasm in ONT mode",
    "assembled with hifiasm v0.25.0-r726 in ONT mode",
)
replace_in_runs(
    paragraphs[9],
    "with NCBI FCS-Adaptor and FCS-GX.",
    "with NCBI FCS-Adaptor and FCS-GX (both v0.5.5), executed as Apptainer v1.4.4 containers "
    "(MD5 27fb2938ed60835a9b41d1ce9533dbf6 and 449690356545bd4ddfc6e9c6153e7de8, respectively).",
)
replace_in_runs(
    paragraphs[9],
    "filtered BAM files were generated with ",
    "filtered BAM files were generated with samtools v1.20 using ",
)
replace_in_runs(
    paragraphs[9],
    "HapHiC (v1.0.7) was then used",
    "HapHiC v1.0.7-4-g90c792c (commit 90c792c45d76a499bee186273585d4242b245343) was then used",
)

# Alignment and quality-control software versions.
replace_in_runs(
    paragraphs[11],
    "Minimap2 was run",
    "Minimap2 v2.30-r1287 was run",
)
replace_in_runs(
    paragraphs[15],
    "with QUAST and seqkit.",
    "with QUAST v5.3.0 and seqkit v2.10.0.",
)
replace_in_runs(
    paragraphs[15],
    "evaluated with Inspector using",
    "evaluated with Inspector v1.3.1 using",
)

# Repeat and gene-annotation software versions.
replace_in_runs(
    paragraphs[17],
    "annotated with EDTA v2.2.2 using ",
    "annotated with EDTA v2.2.2, including RepeatMasker v4.2.2 and LTR_retriever v3.0.4, using ",
)
replace_in_runs(
    paragraphs[17],
    "with panEDTA (",
    "with the panEDTA script distributed with EDTA v2.2.2 (",
)
replace_in_runs(
    paragraphs[17],
    "). The shared library",
    "; script MD5 63c0ec43e011b8d62729836c68ef5e83). The shared library",
)
replace_in_runs(
    paragraphs[19],
    "Protein-coding genes were predicted with Helixer v0.3.5 using",
    "Protein-coding genes were predicted with Helixer v0.3.5 and GeenuFF v0.3.2, "
    "executed in Singularity v3.8.6, using",
)
replace_in_runs(
    paragraphs[19],
    "EviAnn refinement.",
    "EviAnn v2.0.4 refinement (Bioconda package release 2.0.5).",
)

# PSP reassembly used the same recorded samtools installation.
replace_in_runs(
    paragraphs[21],
    "converted to FASTQ with samtools and assembled",
    "converted to FASTQ with samtools v1.20 and assembled",
)

# Restore the omitted review gate in the PPA export sentence.
replace_in_runs(
    paragraphs[27],
    ", GFF3, CDS and protein packages were exported.",
    ", and an approved review status before the GFF3, CDS and protein packages were exported.",
)

document.save(TEMP)
if not TEMP.is_file() or TEMP.stat().st_size == 0:
    raise SystemExit("Updated document was not written")
TEMP.replace(SOURCE)

print(f"Updated: {SOURCE}")
print(f"Backup:  {BACKUP}")
print(f"Encoding repairs: {artifact_counts}")
