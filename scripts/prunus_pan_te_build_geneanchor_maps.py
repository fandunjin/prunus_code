#!/usr/bin/env python3
"""Build auditable GFF-gene to OrthoFinder orthogroup maps.

The map is intentionally separate from the TE-site analysis.  A gene-anchor
run is only meaningful after this file reports a high mapping fraction for
the selected samples.
"""
import argparse, csv, gzip, re
from collections import defaultdict
from pathlib import Path

def opener(path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, encoding="utf-8")
def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as h: return list(csv.DictReader(h, delimiter="\t"))
def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as h:
        w=csv.DictWriter(h, fieldnames=fields, delimiter="\t", extrasaction="ignore"); w.writeheader(); w.writerows(rows)
def parse_attrs(raw):
    d={}
    for item in raw.split(";"):
        if "=" in item:
            k,v=item.split("=",1); d[k]=v
    return d
def parse_genes(path):
    rows=[]
    with opener(path) as h:
        for line in h:
            if not line.strip() or line[0] == "#": continue
            f=line.rstrip("\n").split("\t")
            if len(f)<9 or f[2].lower() not in ("gene", "protein_coding_gene"): continue
            a=parse_attrs(f[8])
            if a.get("gene_biotype") and a["gene_biotype"] not in ("protein_coding", "protein_coding_gene"): continue
            gid=a.get("ID", "")
            if not gid: continue
            aliases={gid, a.get("geneID", ""), a.get("original_id", ""), a.get("gene", "")}
            aliases={x for x in aliases if x}
            rows.append({"gff_gene_id":gid,"seqid":f[0],"start":f[3],"end":f[4],"strand":f[6],"aliases":aliases})
    return rows
def parse_idmap(path):
    aliases={}
    for r in read_tsv(path):
        stable=r.get("stable_id", ""); original=r.get("original_id", "")
        if not stable or not original: continue
        keys={original, original.rsplit("-mRNA",1)[0], re.sub(r"\.\d+$", "", original)}
        if "|" in original: keys.add(original.split("|",1)[0])
        header=r.get("original_header", "")
        for m in re.finditer(r"(?:gene|orig_id|locus)=([^ ]+)", header): keys.add(m.group(1))
        for key in keys: aliases.setdefault(key, stable)
    return aliases
def parse_ogs(path):
    result=defaultdict(set)
    for r in read_tsv(path):
        og=r["Orthogroup"]
        for col,val in r.items():
            if col != "Orthogroup" and val:
                for gene in val.split(", "): result[gene].add(og)
    return result
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--sample", action="append", required=True, help="name,gff,idmap")
    p.add_argument("--orthogroups", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args=p.parse_args(); ogs=parse_ogs(args.orthogroups)
    manifest=[]
    for spec in args.sample:
        f=spec.split(",")
        if len(f)!=3: raise SystemExit("--sample requires name,gff,idmap")
        name,gff,idmap= f[0],Path(f[1]),Path(f[2])
        if not gff.is_file() or not idmap.is_file(): raise SystemExit(f"ERROR: missing input for {name}")
        id_alias=parse_idmap(idmap); rows=[]; status=defaultdict(int)
        for g in parse_genes(gff):
            stable=""
            for alias in g["aliases"]:
                if alias in id_alias: stable=id_alias[alias]; break
            groups=sorted(ogs.get(stable, set())) if stable else []
            state="mapped" if groups else ("stable_id_only" if stable else "unmapped")
            status[state]+=1
            rows.append({"sample":name,"gff_gene_id":g["gff_gene_id"],"seqid":g["seqid"],"start":g["start"],"end":g["end"],"strand":g["strand"],"stable_id":stable,"orthogroups":",".join(groups),"mapping_status":state})
        out=args.output_dir/f"{name}.gene_to_orthogroup.tsv"; write_tsv(out,["sample","gff_gene_id","seqid","start","end","strand","stable_id","orthogroups","mapping_status"],rows)
        total=len(rows); mapped=status["mapped"]
        manifest.append({"sample":name,"genes":total,"mapped":mapped,"stable_id_only":status["stable_id_only"],"unmapped":status["unmapped"],"mapping_fraction":f"{mapped/total:.6f}" if total else "0"})
    write_tsv(args.output_dir/"geneanchor_mapping_qc.tsv",["sample","genes","mapped","stable_id_only","unmapped","mapping_fraction"],manifest)
    print("\n".join("\t".join(str(x[k]) for k in ("sample","genes","mapped","stable_id_only","unmapped","mapping_fraction")) for x in manifest))
if __name__=="__main__": main()
