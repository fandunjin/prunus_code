#!/usr/bin/env python3
"""Add an A-E radial track guide to the intentionally blank top Circos gap.

Circos has no native free-positioned axis annotation.  This script appends a
small SVG overlay after Circos finishes, leaving all data tracks and links
untouched.  It is idempotent so a resumed job never duplicates the guide.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


START = "<!-- PZI_TRACK_GUIDE_START -->"
END = "<!-- PZI_TRACK_GUIDE_END -->"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--svg", type=Path, required=True)
    parser.add_argument("--layout", type=Path, required=True)
    parser.add_argument("--ideogram-radius-factor", type=float, default=0.940)
    parser.add_argument("--ideogram-thickness-p", type=float, default=22.0)
    return parser.parse_args()


def svg_size(svg: str) -> tuple[float, float]:
    viewbox = re.search(r'viewBox="\s*[-+0-9.eE]+\s+[-+0-9.eE]+\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*"', svg)
    if viewbox:
        return float(viewbox.group(1)), float(viewbox.group(2))
    width = re.search(r'\bwidth="([0-9.]+)', svg)
    height = re.search(r'\bheight="([0-9.]+)', svg)
    if width and height:
        return float(width.group(1)), float(height.group(1))
    raise ValueError("Cannot determine the SVG canvas dimensions")


def radial_fraction(value: str) -> float:
    match = re.fullmatch(r"([0-9.]+)r", value.strip())
    if match is None:
        raise ValueError(f"Invalid radial coordinate: {value!r}")
    return float(match.group(1))


def axis_value(value: float) -> str:
    if abs(value) >= 10:
        return f"{value:.0f}"
    # Keep small fractions readable while avoiding unnecessary trailing zeros.
    return f"{value:.3f}".rstrip("0").rstrip(".")


def load_layout(path: Path) -> list[tuple[str, float, float, float, float]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty track layout: {path}")
    rows: list[tuple[str, float, float, float, float]] = []
    with path.open("r", encoding="ascii", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"panel_label", "r0", "r1", "min", "max"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"Track layout lacks {sorted(required)}: {path}")
        for row in reader:
            rows.append(
                (
                    row["panel_label"],
                    radial_fraction(row["r0"]),
                    radial_fraction(row["r1"]),
                    float(row["min"]),
                    float(row["max"]),
                )
            )
    if not rows:
        raise ValueError(f"No track rows parsed from {path}")
    return rows


def make_guide(
    width: float,
    height: float,
    tracks: list[tuple[str, float, float, float, float]],
    ideogram_radius_factor: float = 0.940,
    ideogram_thickness_p: float = 22.0,
) -> str:
    center_x = width / 2.0
    center_y = height / 2.0
    radius = min(width, height) / 2.0
    # Circos unit_parse() resolves values below 1r against the ideogram's
    # inner radius, not its outer radius.  Subtract the ideogram thickness so
    # the appended A-E axes use exactly the same radial reference as plots.
    track_reference_radius = radius * ideogram_radius_factor - ideogram_thickness_p
    if track_reference_radius <= 0:
        raise ValueError("Ideogram inner radius must be positive")
    pieces = [START, '<g id="pzi_track_guide" fill="none" stroke="#6d7379" stroke-width="1.5">']
    for label, r0, r1, minimum, maximum in tracks:
        y_outer = center_y - track_reference_radius * r1
        y_inner = center_y - track_reference_radius * r0
        y_mid = (y_outer + y_inner) / 2.0
        pieces.append(
            f'<line x1="{center_x:.1f}" y1="{y_outer:.1f}" x2="{center_x:.1f}" y2="{y_inner:.1f}"/>'
        )
        # Five evenly spaced ticks make the compact radial axes interpretable
        # without changing any of the underlying Circos tracks.  The top of
        # each axis represents the track maximum and the bottom its minimum.
        span = maximum - minimum
        tick_fractions = (0.0, 0.25, 0.5, 0.75, 1.0)
        for fraction in tick_fractions:
            y = y_outer + fraction * (y_inner - y_outer)
            value = maximum - fraction * span
            pieces.append(
                f'<line x1="{center_x - 8:.1f}" y1="{y:.1f}" '
                f'x2="{center_x + 8:.1f}" y2="{y:.1f}"/>'
            )
            pieces.append(
                f'<text x="{center_x - 14:.1f}" y="{y + 4:.1f}" fill="#596067" stroke="none" '
                'text-anchor="end" font-family="Arial,Helvetica,sans-serif" font-size="9">'
                f'{axis_value(value)}</text>'
            )
        pieces.append(
            f'<text x="{center_x + 18:.1f}" y="{y_mid + 8:.1f}" fill="#3d4349" stroke="none" '
            'font-family="Arial,Helvetica,sans-serif" font-size="17">'
            f'{label}</text>'
        )
    pieces.append("</g>")
    pieces.append(END)
    return "\n".join(pieces)


def main() -> None:
    args = parse_args()
    path = args.svg.resolve()
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty SVG: {path}")
    svg = path.read_text(encoding="utf-8")
    svg = re.sub(re.escape(START) + r".*?" + re.escape(END) + r"\s*", "", svg, flags=re.DOTALL)
    width, height = svg_size(svg)
    tracks = load_layout(args.layout.resolve())
    closing = "</svg>"
    if closing not in svg:
        raise ValueError(f"Not a complete SVG file: {path}")
    svg = svg.rsplit(closing, 1)[0] + make_guide(
        width,
        height,
        tracks,
        args.ideogram_radius_factor,
        args.ideogram_thickness_p,
    ) + "\n" + closing
    path.write_text(svg, encoding="utf-8")
    print(f"[DONE] Added A-E track guide to {path}")


if __name__ == "__main__":
    main()
