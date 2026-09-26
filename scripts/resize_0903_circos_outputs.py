#!/usr/bin/env python3
"""Add 50-Mb ticks and normalize 0903 CircOS outputs to 45 mm square.

The native CircOS outputs are retained as the source artwork.  This post-
processor only adds missing coordinate ticks and changes the output canvas;
the karyotype, tracks, links and colours are not recalculated.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from PIL import Image
from pypdf import PdfReader, PdfWriter


TARGET_MM = 45.0
TARGET_PDF_PT = TARGET_MM / 25.4 * 72.0
PNG_DPI = 600
SPECIES = ("pzi", "psp", "ppa")


def parse_karyotype(path: Path) -> list[tuple[str, int]]:
    rows: list[tuple[str, int]] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            fields = raw.split()
            if len(fields) >= 6 and fields[0] == "chr":
                rows.append((fields[2], int(fields[5])))
    if not rows:
        raise ValueError(f"No chromosomes found in {path}")
    return rows


def polar_angle(x: float, y: float, cx: float, cy: float) -> float:
    return math.atan2(y - cy, x - cx)


def canvas_dimensions(svg: str) -> tuple[float, float]:
    match = re.search(
        r'<svg\s+[^>]*width="([0-9.]+)px"[^>]*height="([0-9.]+)px"', svg
    )
    if match:
        return float(match.group(1)), float(match.group(2))
    match = re.search(r'<svg\s+[^>]*viewBox="0 0 ([0-9.]+) ([0-9.]+)"', svg)
    if match:
        return float(match.group(1)), float(match.group(2))
    raise ValueError("SVG dimensions are missing or not numeric")


def ideogram_spans(svg: str, chromosome_count: int) -> list[tuple[float, float]]:
    """Extract start/end angles from the filled ideogram paths."""
    match = re.search(r'<g id="ideograms">(.*?)</g>', svg, flags=re.S)
    if not match:
        raise ValueError("SVG does not contain the CircOS ideograms group")
    group = match.group(1)
    path_re = re.compile(r'<path d="([^"]+)"\s+style="([^"]+)"\s*/>')
    coordinate_re = re.compile(
        r'^M([0-9.+-]+),([0-9.+-]+)\s+A[^L]+?L([0-9.+-]+),([0-9.+-]+)\s+A'
    )
    paths: list[tuple[str, str, str, str]] = []
    for path_data, style in path_re.findall(group):
        if "fill:rgb" not in style:
            continue
        coordinates = coordinate_re.search(path_data)
        if coordinates:
            paths.append(coordinates.groups())
    if len(paths) != chromosome_count:
        raise ValueError(
            f"Expected {chromosome_count} filled ideograms, found {len(paths)}"
        )
    # CircOS centres the image and uses the outer radius of the ideogram.  The
    # centre can be inferred from the symmetry of the first path coordinates.
    coords = [(float(a), float(b), float(c), float(d)) for a, b, c, d in paths]
    cx = max(max(a, c) for a, _b, c, _d in coords) / 2.0 + min(min(a, c) for a, _b, c, _d in coords) / 2.0
    cy = max(max(b, d) for _a, b, _c, d in coords) / 2.0 + min(min(b, d) for _a, b, _c, d in coords) / 2.0
    # In practice the centre is the midpoint of the SVG canvas.  This avoids
    # propagating rounding noise from the individual path coordinates.
    width, height = canvas_dimensions(svg)
    cx, cy = width / 2.0, height / 2.0
    spans: list[tuple[float, float]] = []
    for x1, y1, x2, y2 in coords:
        start = polar_angle(x1, y1, cx, cy)
        end = polar_angle(x2, y2, cx, cy)
        while end <= start:
            end += 2.0 * math.pi
        spans.append((start, end))
    return spans


def tick_group(svg: str, karyotype: list[tuple[str, int]]) -> str:
    spans = ideogram_spans(svg, len(karyotype))
    width, height = canvas_dimensions(svg)
    cx, cy = width / 2.0, height / 2.0
    # Filled ideogram paths use an outer radius of 0.9r.  Recover it from the
    # first path and use a small outward tick and label offset.
    first = re.search(
        r'<path d="M[0-9.+-]+,[0-9.+-]+\s+A[^L]+?L([0-9.+-]+),([0-9.+-]+)',
        svg,
    )
    if not first:
        raise ValueError("SVG has no ideogram path")
    outer = math.hypot(float(first.group(1)) - cx, float(first.group(2)) - cy)
    tick_outer = outer
    tick_inner = outer + 6.0
    label_radius = outer + 10.0
    elements = ['<g id="local_50mb_ticks" fill="none" stroke="rgb(150,150,150)" stroke-width="1.0">']
    for (seqid, length), (start, end) in zip(karyotype, spans):
        for mb in range(50_000_000, length, 50_000_000):
            angle = start + (end - start) * (mb / length)
            cos_a, sin_a = math.cos(angle), math.sin(angle)
            x1, y1 = cx + tick_outer * cos_a, cy + tick_outer * sin_a
            x2, y2 = cx + tick_inner * cos_a, cy + tick_inner * sin_a
            lx, ly = cx + label_radius * cos_a, cy + label_radius * sin_a
            rotation = math.degrees(angle) % 360.0
            elements.append(
                f'<line x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" y2="{y2:.3f}"/>'
            )
            elements.append(
                f'<text x="{lx:.3f}" y="{ly:.3f}" fill="rgb(150,150,150)" stroke="none" '
                f'font-size="11.7px" font-family="CMUBright-Roman" text-anchor="middle" '
                f'transform="rotate({rotation:.2f},{lx:.3f},{ly:.3f})">{mb // 1_000_000}</text>'
            )
    elements.append('</g>')
    return "".join(elements)


def update_svg(path: Path, karyotype_path: Path) -> None:
    svg = path.read_text(encoding="utf-8")
    native_width, native_height = canvas_dimensions(svg)
    karyotype = parse_karyotype(karyotype_path)
    svg = re.sub(r'<g id="local_50mb_ticks">.*?</g>', "", svg, flags=re.S)
    # CircOS itself writes these labels when the configuration has
    # ``spacing = 50u``.  Add the group only for the older outputs that had
    # origin ticks alone, avoiding duplicate labels on later reruns.
    if not re.search(r'>50</text>', svg):
        group = tick_group(svg, karyotype)
        ideograms_end = svg.find('</g>', svg.find('<g id="ideograms">'))
        if ideograms_end < 0:
            raise ValueError("Cannot locate the end of the ideograms group")
        svg = svg[:ideograms_end] + group + svg[ideograms_end:]
    svg = re.sub(
        r'<svg\s+width="[0-9.]+px"\s+height="[0-9.]+px"',
        f'<svg width="45mm" height="45mm" viewBox="0 0 {native_width:.2f} {native_height:.2f}"',
        svg,
        count=1,
    )
    path.write_text(svg, encoding="utf-8")


def update_png(path: Path) -> None:
    image = Image.open(path).convert("RGB")
    pixels = round(TARGET_MM / 25.4 * PNG_DPI)
    resized = image.resize((pixels, pixels), Image.Resampling.LANCZOS)
    with tempfile.NamedTemporaryFile(prefix="circos_", suffix=".png", dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
    try:
        resized.save(temporary, format="PNG", dpi=(PNG_DPI, PNG_DPI), optimize=True)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def update_pdf(path: Path, svg_path: Path, soffice: str | None) -> None:
    source = path
    if soffice:
        with tempfile.TemporaryDirectory(prefix="circos_svg_pdf_", dir=path.parent) as directory:
            output_dir = Path(directory)
            subprocess.run(
                [soffice, "--headless", "--convert-to", "pdf", "--outdir", str(output_dir), str(svg_path)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
            converted = output_dir / f"{svg_path.stem}.pdf"
            if not converted.is_file() or converted.stat().st_size == 0:
                raise RuntimeError(f"LibreOffice did not create a PDF for {svg_path}")
            _scale_pdf(converted, path)
            return
    _scale_pdf(source, path)


def _scale_pdf(source: Path, target: Path) -> None:
    reader = PdfReader(str(source))
    writer = PdfWriter()
    for page in reader.pages:
        width = float(page.mediabox.width)
        if width <= 0:
            raise ValueError(f"Invalid PDF page width in {source}")
        page.scale_by(TARGET_PDF_PT / width)
        writer.add_page(page)
    with tempfile.NamedTemporaryFile(prefix="circos_", suffix=".pdf", dir=target.parent, delete=False) as handle:
        temporary = Path(handle.name)
    try:
        with temporary.open("wb") as handle:
            writer.write(handle)
        temporary.replace(target)
    finally:
        temporary.unlink(missing_ok=True)


def process_species(root: Path, species: str, soffice: str | None, suffix: str) -> None:
    species_root = root / species
    figures = species_root / "figures"
    stem = f"{species}_ab_circos{suffix}"
    update_svg(figures / f"{stem}.svg", species_root / "data" / "karyotype.txt")
    update_png(figures / f"{stem}.png")
    update_pdf(figures / f"{stem}.pdf", figures / f"{stem}.svg", soffice)
    print(f"[OK] {species}: 50-Mb ticks and 45-mm outputs")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True, help="0903/haploid_repartition/circos root")
    parser.add_argument("--species", nargs="*", choices=SPECIES, default=SPECIES)
    parser.add_argument(
        "--soffice", type=str,
        default=os.environ.get("SOFFICE_BIN") or shutil.which("soffice"),
        help="LibreOffice executable used to rebuild PDF from the updated SVG",
    )
    parser.add_argument(
        "--suffix", type=str, default="",
        help="Filename suffix after <species>_ab_circos, for example _ncrna",
    )
    args = parser.parse_args()
    for species in args.species:
        process_species(args.root.resolve(), species, args.soffice, args.suffix)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
