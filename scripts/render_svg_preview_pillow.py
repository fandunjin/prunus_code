#!/usr/bin/env python3
"""Render SVGs emitted by plot_fig2_gtr_target_dotplots.py to PNG previews."""

from __future__ import annotations

import argparse
import re
import xml.etree.ElementTree as ET
from pathlib import Path

from PIL import Image, ImageColor, ImageDraw, ImageFont


STYLE = {
    "title": {"size": 17, "bold": True, "fill": "#000000"},
    "label": {"size": 11, "bold": False, "fill": "#222222"},
    "tick": {"size": 8, "bold": False, "fill": "#444444"},
    "meta": {"size": 9, "bold": False, "fill": "#555555"},
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", help="SVG files or directories to render")
    parser.add_argument("--scale", type=float, default=1.0)
    return parser.parse_args()


def number(value, default=0.0):
    if value is None:
        return default
    match = re.match(r"[-+0-9.eE]+", value)
    return float(match.group(0)) if match else default


def rgb(value, default="#000000"):
    return ImageColor.getrgb(value or default)


def style_fill(element, default=None):
    style = element.attrib.get("style", "")
    for item in style.split(";"):
        key, separator, value = item.partition(":")
        if separator and key.strip() == "fill":
            return value.strip()
    return element.attrib.get("fill", default)


def font_path(bold):
    windows_fonts = Path("C:/Windows/Fonts")
    candidates = ["arialbd.ttf", "Arial.ttf"] if bold else ["arial.ttf", "Arial.ttf"]
    for name in candidates:
        candidate = windows_fonts / name
        if candidate.exists():
            return str(candidate)
    return None


def load_font(size, bold, scale):
    path = font_path(bold)
    if path:
        return ImageFont.truetype(path, max(1, round(size * scale)))
    return ImageFont.load_default()


def draw_text(draw, element, scale):
    class_name = element.attrib.get("class", "")
    style = STYLE.get(class_name, {"size": 10, "bold": False, "fill": "#000000"})
    font = load_font(style["size"], style["bold"], scale)
    label = "".join(element.itertext())
    x = number(element.attrib.get("x")) * scale
    y = number(element.attrib.get("y")) * scale
    bbox = draw.textbbox((0, 0), label, font=font)
    text_width = bbox[2] - bbox[0]
    anchor = element.attrib.get("text-anchor", "start")
    if anchor == "middle":
        x -= text_width / 2.0
    elif anchor == "end":
        x -= text_width
    # SVG y is the baseline; Pillow y is the top of the text box.
    y -= (bbox[3] - bbox[1])
    draw.text((x, y), label, font=font, fill=rgb(style_fill(element), style["fill"]))


def render(svg_path, scale):
    root = ET.parse(svg_path).getroot()
    width = round(number(root.attrib.get("width")) * scale)
    height = round(number(root.attrib.get("height")) * scale)
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image, "RGBA")

    for element in root:
        tag = element.tag.rsplit("}", 1)[-1]
        attrs = element.attrib
        if tag in {"style"}:
            continue
        if tag == "rect":
            x = number(attrs.get("x")) * scale
            y = number(attrs.get("y")) * scale
            w = number(attrs.get("width"), width / scale) * scale
            h = number(attrs.get("height"), height / scale) * scale
            if attrs.get("width", "").endswith("%"):
                w = width
            if attrs.get("height", "").endswith("%"):
                h = height
            fill = "#FCFCFC" if attrs.get("class") == "frame" else attrs.get("fill", "white")
            outline = "#333333" if attrs.get("class") == "frame" else attrs.get("stroke")
            fill_color = None if fill == "none" else rgb(fill)
            outline_color = None if not outline or outline == "none" else rgb(outline)
            draw.rectangle((x, y, x + w, y + h), fill=fill_color, outline=outline_color,
                           width=max(1, round(number(attrs.get("stroke-width"), 0.8) * scale)))
        elif tag == "line":
            color = "#CBD1D8" if attrs.get("class") == "grid" else attrs.get("stroke", "#000000")
            draw.line(
                (
                    number(attrs.get("x1")) * scale,
                    number(attrs.get("y1")) * scale,
                    number(attrs.get("x2")) * scale,
                    number(attrs.get("y2")) * scale,
                ),
                fill=rgb(color),
                width=max(1, round(number(attrs.get("stroke-width"), 0.55) * scale)),
            )
        elif tag == "circle":
            x = number(attrs.get("cx")) * scale
            y = number(attrs.get("cy")) * scale
            radius = max(0.5, number(attrs.get("r"), 0.75) * scale)
            opacity = max(0.0, min(1.0, number(attrs.get("opacity"), 1.0)))
            color = rgb(attrs.get("fill", "#000000")) + (round(255 * opacity),)
            draw.ellipse((x - radius, y - radius, x + radius, y + radius), fill=color)
        elif tag in {"polyline", "polygon"}:
            points = []
            for point in attrs.get("points", "").split():
                x_value, separator, y_value = point.partition(",")
                if separator:
                    points.append((number(x_value) * scale, number(y_value) * scale))
            if not points:
                continue
            opacity = max(0.0, min(1.0, number(attrs.get("opacity"), number(attrs.get("fill-opacity"), 1.0))))
            fill_value = attrs.get("fill", "none")
            stroke_value = attrs.get("stroke", "none")
            if tag == "polygon" and fill_value != "none":
                draw.polygon(points, fill=rgb(fill_value) + (round(255 * opacity),))
            if stroke_value != "none":
                draw.line(points + ([points[0]] if tag == "polygon" else []), fill=rgb(stroke_value) + (round(255 * opacity),),
                          width=max(1, round(number(attrs.get("stroke-width"), 1.0) * scale)), joint="curve")
        elif tag == "text":
            draw_text(draw, element, scale)

    output_path = svg_path.with_suffix(".png")
    image.save(output_path, dpi=(150, 150), optimize=True)
    print("[OK]", output_path)


def main():
    args = parse_args()
    if args.scale <= 0:
        raise SystemExit("ERROR: --scale must be greater than zero")
    svg_paths = []
    for value in args.inputs:
        path = Path(value)
        if path.is_dir():
            svg_paths.extend(sorted(path.rglob("*.svg")))
        else:
            svg_paths.append(path)
    if not svg_paths:
        raise SystemExit("ERROR: no SVG files found")
    for svg_path in svg_paths:
        render(svg_path, args.scale)


if __name__ == "__main__":
    main()
