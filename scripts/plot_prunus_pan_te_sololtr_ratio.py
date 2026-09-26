#!/usr/bin/env python3
"""Plot solo-LTR turnover and intact LTR-RT superfamily composition."""

import argparse
import csv
import html
from collections import defaultdict
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
from reportlab.lib.colors import HexColor, white
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[1]
GROUPS = [
    ("PZI-A", ["pziA1", "pziA2"], "#C44E52"),
    ("PZI-B", ["pziB1", "pziB2"], "#C44E52"),
    ("PSP-A", ["pspA1", "pspA2"], "#2A9D8F"),
    ("PSP-B", ["pspB1", "pspB2"], "#2A9D8F"),
    ("PPA-A", ["ppaA"], "#4C78A8"),
    ("PPA-B", ["ppaB"], "#4C78A8"),
    ("PCAM", ["pcam"], "#767676"),
    ("PPE", ["ppe"], "#B279A2"),
]
TARGETS = [target for _, targets, _ in GROUPS for target in targets]
FAMILIES = ("Copia", "Gypsy", "unknown")
FAMILY_LABELS = {"Copia": "Copia", "Gypsy": "Gypsy", "unknown": "Unknown"}
FAMILY_COLORS = {"Copia": "#4C78A8", "Gypsy": "#E69F00", "unknown": "#A0A0A0"}


def resolve(path):
    path = Path(path)
    if not path.is_absolute():
        path = next((candidate for candidate in (path, ROOT / path) if candidate.is_file()), ROOT / path)
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def read_ratios(path):
    ratios = {}
    with resolve(path).open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            intact = int(row["quality_filtered_intact_records"])
            if intact <= 0:
                raise ValueError(f"non-positive intact count for {row['target']}")
            ratios[row["target"]] = int(row["final_solo_ltr_records"]) / intact
    missing = [target for target in TARGETS if target not in ratios]
    if missing:
        raise ValueError("missing ratio targets: " + ", ".join(missing))
    return ratios


def read_compositions(path):
    counts = defaultdict(lambda: defaultdict(int))
    with resolve(path).open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            family = row["superfamily"].strip()
            family = family if family in ("Copia", "Gypsy") else "unknown"
            counts[row["target"]][family] += int(row["intact_count"])
    compositions = {}
    for target in TARGETS:
        total = sum(counts[target].get(family, 0) for family in FAMILIES)
        if total <= 0:
            raise ValueError(f"no intact LTR-RT counts for {target}")
        compositions[target] = {family: counts[target].get(family, 0) / total for family in FAMILIES}
    return compositions


def group_data(ratios, compositions):
    ratio_data = []
    comp_data = []
    for label, targets, color in GROUPS:
        ratio_data.append((label, [ratios[target] for target in targets], color))
        comp_data.append((label, {family: [compositions[target][family] for target in targets] for family in FAMILIES}))
    return ratio_data, comp_data


def write_source_table(path, comp_data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["group", "Copia_fraction", "Gypsy_fraction", "Unknown_fraction"])
        for label, values in comp_data:
            writer.writerow([label] + [f"{sum(values[family]) / len(values[family]):.8f}" for family in FAMILIES])


def geometry():
    width, height = 2.15 * 72, 4.15 * 72
    return {
        "width": width, "height": height, "left": width * 0.28, "right": width * 0.96,
        "ratio_top": height * 0.04, "ratio_bottom": height * 0.43,
        "comp_top": height * 0.61, "comp_bottom": height * 0.86,
    }


def x_positions(left, right, count):
    return [left + (index + 0.55) * (right - left) / (count + 0.1) for index in range(count)]


def font(size, bold=False):
    candidates = [
        r"C:\Windows\Fonts\arialbd.ttf" if bold else r"C:\Windows\Fonts\arial.ttf",
        r"C:\Windows\Fonts\segoeui.ttf",
        r"C:\Windows\Fonts\DejaVuSans.ttf",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, max(1, round(size)))
        except OSError:
            pass
    return ImageFont.load_default()


def y_ratio(value, top, bottom):
    return bottom - value / 3.65 * (bottom - top)


def y_comp(value, top, bottom):
    return bottom - value * (bottom - top)


def svg_text(x, y, value, size=7, anchor="start", rotate=None):
    transform = f' transform="rotate({rotate:g} {x:g} {y:g})"' if rotate is not None else ""
    return (f'<text x="{x:.3f}" y="{y:.3f}" font-family="Arial, Helvetica, sans-serif" '
            f'font-size="{size:g}pt" fill="#222222" text-anchor="{anchor}"{transform}>'
            f'{html.escape(str(value))}</text>')


def svg_figure(ratio_data, comp_data):
    g = geometry(); width, height = g["width"], g["height"]
    left, right = g["left"], g["right"]
    rt, rb, ct, cb = g["ratio_top"], g["ratio_bottom"], g["comp_top"], g["comp_bottom"]
    xs = x_positions(left, right, len(GROUPS)); bar = (right - left) / (len(GROUPS) + 0.1) * 0.62
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.1f}pt" height="{height:.1f}pt" viewBox="0 0 {width:.1f} {height:.1f}">', '<rect width="100%" height="100%" fill="white"/>']
    out.append(svg_text(left - 24, rt + 4, "a", size=7, anchor="middle"))
    for tick in (0, .5, 1, 1.5, 2, 2.5, 3, 3.5):
        y = y_ratio(tick, rt, rb); out.append(f'<line x1="{left:.3f}" y1="{y:.3f}" x2="{right:.3f}" y2="{y:.3f}" stroke="#e8e8e8" stroke-width="0.45"/>'); out.append(svg_text(left - 4, y + 2.4, tick, size=7, anchor="end"))
    y = y_ratio(1, rt, rb); out.append(f'<line x1="{left:.3f}" y1="{y:.3f}" x2="{right:.3f}" y2="{y:.3f}" stroke="#a7a7a7" stroke-width="0.65" stroke-dasharray="2,2"/>')
    out += [f'<line x1="{left:.3f}" y1="{rt:.3f}" x2="{left:.3f}" y2="{rb:.3f}" stroke="#4d4d4d" stroke-width="0.7"/>', f'<line x1="{left:.3f}" y1="{rb:.3f}" x2="{right:.3f}" y2="{rb:.3f}" stroke="#4d4d4d" stroke-width="0.7"/>', svg_text(18, (rt + rb) / 2, "solo-LTR / intact LTR-RT", size=7, anchor="middle", rotate=-90)]
    for i, (_, values, color) in enumerate(ratio_data):
        x = xs[i]; mean = sum(values) / len(values); low, high = y_ratio(min(values), rt, rb), y_ratio(max(values), rt, rb); ym = y_ratio(mean, rt, rb)
        out.append(f'<rect x="{x - bar/2:.3f}" y="{ym:.3f}" width="{bar:.3f}" height="{rb - ym:.3f}" fill="{color}" stroke="#272727" stroke-width="0.45"/>'); out.append(f'<line x1="{x:.3f}" y1="{high:.3f}" x2="{x:.3f}" y2="{low:.3f}" stroke="#272727" stroke-width="0.65"/>')
        cap = bar * .55
        for yy in (low, high): out.append(f'<line x1="{x-cap/2:.3f}" y1="{yy:.3f}" x2="{x+cap/2:.3f}" y2="{yy:.3f}" stroke="#272727" stroke-width="0.65"/>')
    out.append(svg_text(left - 24, ct + 4, "b", size=7, anchor="middle"))
    for tick in (0, .25, .5, .75, 1):
        y = y_comp(tick, ct, cb); out.append(f'<line x1="{left:.3f}" y1="{y:.3f}" x2="{right:.3f}" y2="{y:.3f}" stroke="#e8e8e8" stroke-width="0.45"/>'); out.append(svg_text(left - 4, y + 2.4, round(tick * 100), size=6.4, anchor="end"))
    out += [f'<line x1="{left:.3f}" y1="{ct:.3f}" x2="{left:.3f}" y2="{cb:.3f}" stroke="#4d4d4d" stroke-width="0.7"/>', f'<line x1="{left:.3f}" y1="{cb:.3f}" x2="{right:.3f}" y2="{cb:.3f}" stroke="#4d4d4d" stroke-width="0.7"/>', svg_text(18, (ct + cb) / 2, "intact LTR-RT (%)", size=6.7, anchor="middle", rotate=-90)]
    for i, (_, values) in enumerate(comp_data):
        x = xs[i]; bottom = 0
        for family in FAMILIES:
            value = sum(values[family]) / len(values[family]); top = bottom + value; y1, y0 = y_comp(top, ct, cb), y_comp(bottom, ct, cb); out.append(f'<rect x="{x-bar/2:.3f}" y="{y1:.3f}" width="{bar:.3f}" height="{y0-y1:.3f}" fill="{FAMILY_COLORS[family]}" stroke="white" stroke-width="0.35"/>'); bottom = top
    lx, ly = left + 2, ct - 10
    for family in FAMILIES:
        out.append(f'<rect x="{lx:.3f}" y="{ly-5:.3f}" width="5.5" height="5.5" fill="{FAMILY_COLORS[family]}"/>'); out.append(svg_text(lx + 7, ly, FAMILY_LABELS[family], size=5.7)); lx += 34 if family != "unknown" else 0
    for i, (label, _, _) in enumerate(GROUPS): out.append(svg_text(xs[i], cb + 16, label, size=5.8, anchor="end", rotate=60))
    out.append("</svg>")
    return "\n".join(out)


def draw_raster(path, ratio_data, comp_data, dpi, tiff=False):
    scale = dpi / 72; gp = geometry(); width, height = round(gp["width"] * scale), round(gp["height"] * scale); g = {key: value * scale for key, value in gp.items()}; left, right = g["left"], g["right"]; rt, rb, ct, cb = g["ratio_top"], g["ratio_bottom"], g["comp_top"], g["comp_bottom"]; xs = x_positions(left, right, len(GROUPS)); bar = (right-left)/(len(GROUPS)+.1)*.62
    image = Image.new("RGB", (width, height), "white"); draw = ImageDraw.Draw(image); text_font = font(7 * scale); small_font = font(5.8 * scale)
    draw.text((left - 24 * scale, rt + 4 * scale), "a", fill="#222222", font=text_font, anchor="mm")
    for tick in (0, .5, 1, 1.5, 2, 2.5, 3, 3.5):
        y = y_ratio(tick, rt, rb); draw.line((left, y, right, y), fill="#e8e8e8", width=max(1, round(.45*scale))); draw.text((left-4*scale, y), f"{tick:g}", fill="#222222", font=text_font, anchor="ra")
    y = y_ratio(1, rt, rb); dash = max(2, round(2*scale))
    for start in range(round(left), round(right), dash*2): draw.line((start, y, min(start+dash, right), y), fill="#a7a7a7", width=max(1, round(.65*scale)))
    draw.line((left, rt, left, rb), fill="#4d4d4d", width=max(1, round(.7*scale))); draw.line((left, rb, right, rb), fill="#4d4d4d", width=max(1, round(.7*scale)))
    layer = Image.new("RGBA", (round(180*scale), round(15*scale)), (255,255,255,0)); ImageDraw.Draw(layer).text((0,0), "solo-LTR / intact LTR-RT", fill="#222222", font=text_font); layer = layer.rotate(90, expand=True); image.paste(layer, (round(7*scale), round((rt+rb)/2-layer.height/2)), layer)
    for i, (_, values, color) in enumerate(ratio_data):
        x=xs[i]; mean=sum(values)/len(values); low,high=y_ratio(min(values),rt,rb),y_ratio(max(values),rt,rb); ym=y_ratio(mean,rt,rb); draw.rectangle((x-bar/2,ym,x+bar/2,rb),fill=color,outline="#272727",width=max(1,round(.45*scale))); draw.line((x,high,x,low),fill="#272727",width=max(1,round(.65*scale))); cap=bar*.55
        for yy in (low,high): draw.line((x-cap/2,yy,x+cap/2,yy),fill="#272727",width=max(1,round(.65*scale)))
    draw.text((left - 24 * scale, ct + 4 * scale), "b", fill="#222222", font=text_font, anchor="mm")
    for tick in (0,.25,.5,.75,1):
        y=y_comp(tick,ct,cb); draw.line((left,y,right,y),fill="#e8e8e8",width=max(1,round(.45*scale))); draw.text((left-4*scale,y),str(round(tick*100)),fill="#222222",font=small_font,anchor="ra")
    draw.line((left,ct,left,cb),fill="#4d4d4d",width=max(1,round(.7*scale))); draw.line((left,cb,right,cb),fill="#4d4d4d",width=max(1,round(.7*scale)))
    layer=Image.new("RGBA",(round(130*scale),round(15*scale)),(255,255,255,0)); ImageDraw.Draw(layer).text((0,0),"intact LTR-RT (%)",fill="#222222",font=small_font); layer=layer.rotate(90,expand=True); image.paste(layer,(round(8*scale),round((ct+cb)/2-layer.height/2)),layer)
    for i,(_,values) in enumerate(comp_data):
        x=xs[i]; bottom=0
        for family in FAMILIES:
            value=sum(values[family])/len(values[family]); top=bottom+value; draw.rectangle((x-bar/2,y_comp(top,ct,cb),x+bar/2,y_comp(bottom,ct,cb)),fill=FAMILY_COLORS[family],outline="white",width=max(1,round(.35*scale))); bottom=top
    lx,ly=left+2*scale,ct-10*scale
    for family in FAMILIES:
        draw.rectangle((lx,ly-5*scale,lx+5.5*scale,ly+.5*scale),fill=FAMILY_COLORS[family]); draw.text((lx+7*scale,ly),FAMILY_LABELS[family],fill="#222222",font=small_font,anchor="lm"); lx += 34*scale if family != "unknown" else 0
    for i,(label,_,_) in enumerate(GROUPS):
        bbox=draw.textbbox((0,0),label,font=small_font); layer=Image.new("RGBA",(bbox[2]-bbox[0]+round(8*scale),bbox[3]-bbox[1]+round(8*scale)),(255,255,255,0)); ImageDraw.Draw(layer).text((4*scale,4*scale),label,fill="#222222",font=small_font); layer=layer.rotate(60,expand=True,resample=Image.Resampling.BICUBIC); image.paste(layer,(min(max(0,round(xs[i]-layer.width*.2)),width-layer.width-1),round(cb+4*scale)),layer)
    kwargs={"dpi":(dpi,dpi)}
    if tiff: kwargs.update(format="TIFF", compression="tiff_lzw")
    image.save(path,**kwargs)


def draw_pdf(path, ratio_data, comp_data):
    g=geometry(); width,height=g["width"],g["height"]; left,right=g["left"],g["right"]; rt,rb,ct,cb=g["ratio_top"],g["ratio_bottom"],g["comp_top"],g["comp_bottom"]; xs=x_positions(left,right,len(GROUPS)); bar=(right-left)/(len(GROUPS)+.1)*.62; pdf=canvas.Canvas(str(path),pagesize=(width,height),pageCompression=1); pdf.setFillColor(white); pdf.rect(0,0,width,height,stroke=0,fill=1)
    def py(y): return height-y
    pdf.setFillColor(HexColor("#222222")); pdf.setFont("Helvetica",7); pdf.drawCentredString(left-24,py(rt+4),"a")
    for tick in (0,.5,1,1.5,2,2.5,3,3.5):
        y=y_ratio(tick,rt,rb); pdf.setStrokeColor(HexColor("#e8e8e8")); pdf.setLineWidth(.45); pdf.line(left,py(y),right,py(y)); pdf.setFillColor(HexColor("#222222")); pdf.drawRightString(left-4,py(y)-2.3,f"{tick:g}")
    y=y_ratio(1,rt,rb); pdf.setStrokeColor(HexColor("#a7a7a7")); pdf.setDash(2,2); pdf.setLineWidth(.65); pdf.line(left,py(y),right,py(y)); pdf.setDash(); pdf.setStrokeColor(HexColor("#4d4d4d")); pdf.setLineWidth(.7); pdf.line(left,py(rt),left,py(rb)); pdf.line(left,py(rb),right,py(rb)); pdf.saveState(); pdf.translate(18,py((rt+rb)/2)); pdf.rotate(90); pdf.drawCentredString(0,0,"solo-LTR / intact LTR-RT"); pdf.restoreState()
    for i,(_,values,color) in enumerate(ratio_data):
        x=xs[i]; mean=sum(values)/len(values); low,high=y_ratio(min(values),rt,rb),y_ratio(max(values),rt,rb); ym=y_ratio(mean,rt,rb); pdf.setFillColor(HexColor(color)); pdf.setStrokeColor(HexColor("#272727")); pdf.setLineWidth(.45); pdf.rect(x-bar/2,py(rb),bar,py(ym)-py(rb),stroke=1,fill=1); pdf.setLineWidth(.65); pdf.line(x,py(high),x,py(low)); cap=bar*.55; pdf.line(x-cap/2,py(low),x+cap/2,py(low)); pdf.line(x-cap/2,py(high),x+cap/2,py(high))
    pdf.setFillColor(HexColor("#222222")); pdf.setFont("Helvetica",7); pdf.drawCentredString(left-24,py(ct+4),"b")
    for tick in (0,.25,.5,.75,1):
        y=y_comp(tick,ct,cb); pdf.setStrokeColor(HexColor("#e8e8e8")); pdf.setLineWidth(.45); pdf.line(left,py(y),right,py(y)); pdf.setFillColor(HexColor("#222222")); pdf.setFont("Helvetica",6.4); pdf.drawRightString(left-4,py(y)-2.1,str(round(tick*100)))
    pdf.setStrokeColor(HexColor("#4d4d4d")); pdf.setLineWidth(.7); pdf.line(left,py(ct),left,py(cb)); pdf.line(left,py(cb),right,py(cb)); pdf.saveState(); pdf.translate(18,py((ct+cb)/2)); pdf.rotate(90); pdf.drawCentredString(0,0,"intact LTR-RT (%)"); pdf.restoreState()
    for i,(_,values) in enumerate(comp_data):
        x=xs[i]; bottom=0
        for family in FAMILIES:
            value=sum(values[family])/len(values[family]); top=bottom+value; pdf.setFillColor(HexColor(FAMILY_COLORS[family])); pdf.setStrokeColor(white); pdf.setLineWidth(.35); pdf.rect(x-bar/2,py(y_comp(bottom,ct,cb)),bar,py(y_comp(top,ct,cb))-py(y_comp(bottom,ct,cb)),stroke=1,fill=1); bottom=top
    lx,ly=left+2,ct-10; pdf.setFont("Helvetica",5.7)
    for family in FAMILIES:
        pdf.setFillColor(HexColor(FAMILY_COLORS[family])); pdf.rect(lx,py(ly),5.5,5.5,stroke=0,fill=1); pdf.setFillColor(HexColor("#222222")); pdf.drawString(lx+7,py(ly)+.3,FAMILY_LABELS[family]); lx += 34 if family != "unknown" else 0
    pdf.setFillColor(HexColor("#222222")); pdf.setFont("Helvetica",5.8)
    for i,(label,_,_) in enumerate(GROUPS): pdf.saveState(); pdf.translate(xs[i]+2,py(cb+16)); pdf.rotate(60); pdf.drawRightString(0,0,label); pdf.restoreState()
    pdf.showPage(); pdf.save()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",default="analysis/prunus_pan_te/followup_sololtr_v2/sololtr_summary/sololtr_qc.tsv")
    parser.add_argument("--family-input",default="analysis/prunus_pan_te/followup_sololtr_v2/sololtr_summary/sololtr_by_target_lineage.tsv")
    parser.add_argument("--output-dir",default=ROOT/"analysis/prunus_pan_te/figures_v3")
    parser.add_argument("--stem",default="supp_te_sololtr_ratio_family")
    args=parser.parse_args(); output=Path(args.output_dir); output=output if output.is_absolute() else ROOT/output; output.mkdir(parents=True,exist_ok=True)
    ratios=read_ratios(args.input); compositions=read_compositions(args.family_input); ratio_data,comp_data=group_data(ratios,compositions)
    write_source_table(output/f"{args.stem}_source.tsv",comp_data)
    (output/f"{args.stem}.svg").write_text(svg_figure(ratio_data,comp_data),encoding="utf-8")
    draw_pdf(output/f"{args.stem}.pdf",ratio_data,comp_data); draw_raster(output/f"{args.stem}.png",ratio_data,comp_data,300); draw_raster(output/f"{args.stem}.tiff",ratio_data,comp_data,600,True)
    print(f"Wrote {args.stem} outputs to {output}")


if __name__ == "__main__":
    main()
