#!/usr/bin/env python3
"""Add a compact central legend for the 0903 CircOS data tracks.

CircOS draws the tracks themselves but has no portable, version-independent
way to place a free legend in the central blank area.  This small post-process
adds editable SVG text and coloured swatches without changing any CircOS
geometry or data files.  ``--include-ncrna`` adds the optional ncRNA-density
track used by the PZI enhanced figure.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


BASE_TRACKS = (
    ("gene density", None),
    ("total TE", "#708734"),
    ("LTR / Gypsy", "#c45b52"),
    ("LTR / Copia", "#4c78a8"),
    ("TIR / DNA", "#7e58a0"),
    ("Helitron", "#289691"),
    ("GC", "#2a586c"),
)
NCRNA_TRACK = ("ncRNA density", "#a6365f")


def tracks_for(include_ncrna: bool):
    return BASE_TRACKS + ((NCRNA_TRACK,) if include_ncrna else ())


def markers_for(include_ncrna: bool):
    markers = ("i", "ii", "iii", "iv", "v", "vi", "vii")
    return markers + (("viii",) if include_ncrna else ())


def svg_dimension(svg: str, name: str) -> float:
    match = re.search(rf'\b{name}="([0-9.]+)', svg)
    if not match:
        raise ValueError(f"SVG is missing numeric {name} dimension")
    return float(match.group(1))


def legend_markup(width: float, height: float, include_ncrna: bool = False) -> str:
    tracks = tracks_for(include_ncrna)
    markers = markers_for(include_ncrna)
    cx, cy = width / 2.0, height / 2.0
    # Keep the original proportions for the large CircOS canvas, while
    # scaling the legend down for the 45-mm export canvas.  The lower bound
    # leaves enough room for seven labels at approximately 4.5-pt text.
    scale = min(1.0, max(0.30, width / 540.0))
    row_height = 43.0 * scale
    box_width = 360.0 * scale
    box_height = (row_height * len(tracks) + 42.0 * scale)
    x = cx - box_width / 2.0
    y = cy - box_height / 2.0
    swatch_height = 24.0 * scale
    swatch_width = 58.0 * scale
    marker_size = 17.0 * scale
    label_size = 19.0 * scale
    hint_size = 12.0 * scale
    parts = [
        '<g id="central_track_legend" font-family="Arial,Helvetica,sans-serif">',
        '<defs><linearGradient id="density_gradient" x1="0" y1="0" x2="1" y2="0">'
        '<stop offset="0%" stop-color="#ffffcc"/><stop offset="20%" stop-color="#ffec8b"/>'
        '<stop offset="40%" stop-color="#c6dbee"/><stop offset="60%" stop-color="#6a97c6"/>'
        '<stop offset="80%" stop-color="#4364aa"/><stop offset="100%" stop-color="#302174"/>'
        '</linearGradient></defs>',
        f'<rect x="{x:.1f}" y="{y:.1f}" width="{box_width:.1f}" height="{box_height:.1f}" '
        f'rx="{10.0 * scale:.1f}" fill="#ffffff" fill-opacity="0.94" stroke="#d7dce0" '
        f'stroke-width="{2.0 * scale:.1f}"/>',
    ]
    first_y = y + 34.0 * scale
    for index, (label, colour) in enumerate(tracks):
        yy = first_y + index * row_height
        parts.append(
            f'<text x="{x + 30 * scale:.1f}" y="{yy + 5 * scale:.1f}" fill="#4e555b" stroke="none" '
            f'font-size="{marker_size:.1f}px" text-anchor="end">'
            f'{markers[index]}</text>'
        )
        if index == 0:
            parts.append(
                f'<rect x="{x + 52 * scale:.1f}" y="{yy - 12 * scale:.1f}" '
                f'width="{swatch_width:.1f}" height="{swatch_height:.1f}" '
                f'fill="url(#density_gradient)" stroke="#b8c1c8" stroke-width="{scale:.1f}"/>'
            )
        else:
            parts.append(
                f'<rect x="{x + 52 * scale:.1f}" y="{yy - 12 * scale:.1f}" '
                f'width="{24 * scale:.1f}" height="{24 * scale:.1f}" fill="{colour}" stroke="none"/>'
            )
        parts.append(
            f'<text x="{x + (126 if index == 0 else 94) * scale:.1f}" y="{yy + 5 * scale:.1f}" '
            f'fill="#2f3438" stroke="none" font-size="{label_size:.1f}px">{label}</text>'
        )
        if index == 0:
            parts.append(
                f'<text x="{x + 52 * scale:.1f}" y="{yy + 28 * scale:.1f}" fill="#6b747b" '
                f'stroke="none" font-size="{hint_size:.1f}px">low</text>'
            )
            parts.append(
                f'<text x="{x + 110 * scale:.1f}" y="{yy + 28 * scale:.1f}" fill="#6b747b" '
                f'stroke="none" font-size="{hint_size:.1f}px" text-anchor="end">high</text>'
            )
    parts.append('</g>')
    return "".join(parts)


def draw_png(path: Path, include_ncrna: bool = False) -> None:
    """Draw the same legend on the native CircOS PNG without re-rasterizing SVG."""
    image = Image.open(path).convert("RGBA")
    width, height = image.size
    draw = ImageDraw.Draw(image, "RGBA")
    cx, cy = width / 2.0, height / 2.0
    scale = min(1.0, max(0.30, width / 540.0))
    row_height = 43.0 * scale
    box_width = 360.0 * scale
    tracks = tracks_for(include_ncrna)
    box_height = row_height * len(tracks) + 42.0 * scale
    x = cx - box_width / 2.0
    y = cy - box_height / 2.0
    draw.rounded_rectangle(
        (x, y, x + box_width, y + box_height),
        radius=max(2, round(10.0 * scale)),
        fill=(255, 255, 255, 242),
        outline=(215, 220, 224, 255),
        width=max(1, round(2.0 * scale)),
    )
    font_path = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
    font_size = max(5, round(19.0 * scale))
    marker_font = ImageFont.truetype(font_path, max(5, round(17.0 * scale)))
    label_font = ImageFont.truetype(font_path, font_size)
    first_y = y + 34.0 * scale
    for index, (label, colour) in enumerate(tracks):
        yy = first_y + index * row_height
        draw.text((x + 30.0 * scale, yy), markers_for(include_ncrna)[index], fill=(78, 85, 91, 255), font=marker_font, anchor="ra")
        if index == 0:
            gradient_left, gradient_right = x + 52.0 * scale, x + 110.0 * scale
            gradient_top, gradient_bottom = yy - 12.0 * scale, yy + 12.0 * scale
            stops = ((255, 255, 204), (255, 236, 139), (198, 219, 239), (106, 151, 198), (67, 100, 170), (48, 33, 116))
            for xx in range(int(gradient_left), int(gradient_right) + 1):
                frac = (xx - gradient_left) / max(1, gradient_right - gradient_left)
                position = frac * (len(stops) - 1)
                lo = min(len(stops) - 2, int(position))
                hi = lo + 1
                local = position - lo
                colour_rgb = tuple(round(stops[lo][j] * (1 - local) + stops[hi][j] * local) for j in range(3))
                draw.line((xx, gradient_top, xx, gradient_bottom), fill=colour_rgb + (255,))
            draw.rectangle((gradient_left, gradient_top, gradient_right, gradient_bottom), outline=(184, 193, 200, 255), width=1)
            draw.text((gradient_left, yy + 28.0 * scale), "low", fill=(107, 116, 123, 255), font=marker_font, anchor="la")
            draw.text((gradient_right, yy + 28.0 * scale), "high", fill=(107, 116, 123, 255), font=marker_font, anchor="ra")
            draw.text((x + 126.0 * scale, yy), label, fill=(47, 52, 56, 255), font=label_font, anchor="la")
        else:
            colour_hex = colour.lstrip("#")
            swatch = tuple(int(colour_hex[i : i + 2], 16) for i in (0, 2, 4))
            draw.rectangle((x + 52.0 * scale, yy - 12.0 * scale, x + 76.0 * scale, yy + 12.0 * scale), fill=swatch + (255,))
            draw.text((x + 94.0 * scale, yy), label, fill=(47, 52, 56, 255), font=label_font, anchor="la")
    image.convert("RGB").save(path, format="PNG", optimize=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("svg", type=Path)
    parser.add_argument("--png", type=Path, help="Native CircOS PNG to annotate with the same legend")
    parser.add_argument("--include-ncrna", action="store_true", help="Add the ncRNA-density track to the legend")
    args = parser.parse_args()
    path = args.svg.resolve()
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(path)
    svg = path.read_text(encoding="utf-8")
    if "id=\"central_track_legend\"" in svg:
        print(f"[SKIP] Legend already present: {path}")
        return 0
    if not re.search(r"</svg>\s*$", svg):
        raise ValueError(f"Not a complete SVG file: {path}")
    width = svg_dimension(svg, "width")
    height = svg_dimension(svg, "height")
    svg = re.sub(r"</svg>\s*$", legend_markup(width, height, args.include_ncrna) + "</svg>\n", svg)
    path.write_text(svg, encoding="utf-8")
    if args.png:
        png = args.png.resolve()
        if not png.is_file() or png.stat().st_size == 0:
            raise FileNotFoundError(png)
        draw_png(png, args.include_ncrna)
    print(f"[DONE] Added central track legend: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
