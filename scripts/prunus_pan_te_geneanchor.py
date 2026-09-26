#!/usr/bin/env python3
"""Test TE presence between conserved upstream/downstream gene anchors.

This is deliberately conservative: missing anchors, discordant order and
multiple target placements are reported as unresolved rather than interpreted
as absence.  Coordinates are 1-based inclusive in the output.
"""
import argparse, csv, gzip, re
from bisect import bisect_left, bisect_right
from collections import defaultdict
from pathlib import Path

def op(path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, encoding="utf-8")
def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as h: return list(csv.DictReader(h, delimiter="\t"))
def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as h:
        w=csv.DictWriter(h, fieldnames=fields, delimiter="\t", extrasaction="ignore"); w.writeheader(); w.writerows(rows)
def attrs(s):
    d={}
    for x in s.split(";"):
        if "=" in x:
            k,v=x.split("=",1); d[k]=v
    return d
def parse_gff(path):
    out=defaultdict(list)
    with op(path) as h:
        for line in h:
            if not line.strip() or line[0]=='#': continue
            f=line.rstrip().split("\t")
            if len(f)<9 or f[2].lower() not in ("gene","protein_coding_gene"): continue
            try: start,end=int(f[3]),int(f[4])
            except ValueError: continue
            a=attrs(f[8])
            if a.get("gene_biotype") and a["gene_biotype"] not in ("protein_coding", "protein_coding_gene"):
                continue
            gid=a.get("ID","")
            if not gid: continue
            out[f[0]].append({"id":gid,"aliases":set(x for x in (a.get("geneID"), a.get("original_id"), a.get("gene")) if x),"start":start,"end":end,"strand":f[6],"seqid":f[0]})
    for seq in out: out[seq].sort(key=lambda x:x["start"])
    return out
def aliases(gff):
    result={}
    for seq, genes in gff.items():
        for g in genes:
            gid=g["id"]; result[gid]=g
            for alias in g.get("aliases", set()): result[alias]=g
            result[gid.split(".1",1)[0]]=g
            result[gid.replace("_edta_masked","")]=g
    return result
def parse_idmap(path):
    result={}
    if not path or not Path(path).is_file(): return result
    for r in read_tsv(Path(path)):
        stable=r.get("stable_id",""); original=r.get("original_id","")
        if stable and original:
            result[original]=stable; result[original.rsplit("-mRNA",1)[0]]=stable
            if re.search(r"\.\d+$", original): result[re.sub(r"\.\d+$", "", original)]=stable
            if "|" in original: result[original.split("|", 1)[0]]=stable
            header=r.get("original_header","")
            m=re.search(r"gene=([^ ]+)", header)
            if m: result[m.group(1)]=stable
    return result
def parse_orthogroups(path):
    rows=read_tsv(path); result=defaultdict(set)
    for r in rows:
        og=r["Orthogroup"]
        for col,val in r.items():
            if col=="Orthogroup" or not val: continue
            for gene in val.split(", "): result[gene].add(og)
    return result
def parse_rmout(path, families):
    out=defaultdict(list)
    if not Path(path).is_file(): return out
    with op(path) as h:
        for line in h:
            f=line.split()
            if len(f)<11 or not f[0].isdigit(): continue
            fam=f[9].split("#",1)[0]
            if fam not in families: continue
            try: start,end=int(f[5]),int(f[6])
            except ValueError: continue
            if end>=start: out[fam].append((f[4],start,end))
    return out
def parse_gene_map(path):
    result={}
    if not path.is_file(): return result
    for r in read_tsv(path):
        groups={x for x in r.get("orthogroups", "").split(",") if x}
        result[r.get("gff_gene_id", "")]=groups
    return result
def nearest(genes, start, end):
    ups=[g for g in genes if g["end"] < start]
    dns=[g for g in genes if g["start"] > end]
    return (ups[-1] if ups else None, dns[0] if dns else None)
def build_target_index(gff, ogmap, rmout):
    """Pre-index orthogroup anchors and TE intervals once per target."""
    by_seq_og=defaultdict(lambda: defaultdict(list))
    for seq, genes in gff.items():
        for g in genes:
            for og in ogmap.get(g["id"], set()):
                by_seq_og[seq][og].append(g)
    family_by_seq=defaultdict(lambda: defaultdict(list))
    for fam, loci in rmout.items():
        for seq, start, end in loci:
            family_by_seq[fam][seq].append((start, end))
    interval_index={}
    for fam, seqs in family_by_seq.items():
        interval_index[fam]={}
        for seq, intervals in seqs.items():
            intervals=sorted(intervals)
            interval_index[fam][seq]=(intervals, [x[0] for x in intervals])
    return by_seq_og, interval_index
def target_interval(index, up_ogs, dn_ogs, family, max_gap):
    by_seq_og, interval_index=index
    hits=[]; anchor_pairs=0
    for seq, byog in by_seq_og.items():
        for uo in up_ogs:
            for do in dn_ogs:
                for u in byog.get(uo,[]):
                    for d in byog.get(do,[]):
                        if u["start"] < d["start"] and d["start"]-u["end"] <= max_gap:
                            anchor_pairs += 1
                            if seq not in interval_index.get(family, {}):
                                continue
                            intervals, starts=interval_index[family][seq]
                            left=bisect_right(starts, u["end"])
                            right=bisect_left(starts, d["start"])
                            for a,b in intervals[left:right]:
                                if b < d["start"]: hits.append((u,d,a,b))
    return hits, anchor_pairs
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--query", action="append", required=True, help="name,gff,rmout,idmap,orthofinder_column")
    p.add_argument("--targets", type=Path, required=True, help="TSV name,gff,rmout,idmap,group")
    p.add_argument("--orthogroups", type=Path, required=True)
    p.add_argument("--families", type=Path, required=True, help="curated_te_driver_families.tsv")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--gene-map-dir", type=Path, default=None)
    p.add_argument("--max-gap", type=int, default=1000000)
    args=p.parse_args(); drivers=read_tsv(args.families); families={r["pan_family"] for r in drivers if r.get("manual_decision")=="include"}
    if not families: families={r["pan_family"] for r in drivers}
    ogs=parse_orthogroups(args.orthogroups)
    configs=[]
    for q in args.query:
        f=q.split(",");
        if len(f)!=5: raise SystemExit("--query requires name,gff,rmout,idmap,orthofinder_column")
        configs.append({"name":f[0],"gff":Path(f[1]),"rmout":Path(f[2]),"idmap":Path(f[3]),"column":f[4]})
    targets=read_tsv(args.targets)
    parsed={}
    for c in configs+[{"name":r["name"],"gff":Path(r["gff"]),"rmout":Path(r["rmout"]),"idmap":Path(r.get("idmap", "")),"column":r.get("column","")} for r in targets]:
        if not c["gff"].is_file(): raise SystemExit(f"ERROR: missing GFF {c['gff']}")
        if not c["rmout"].is_file(): raise SystemExit(f"ERROR: missing RepeatMasker .out {c['rmout']}")
        if c["idmap"] and not c["idmap"].is_file(): raise SystemExit(f"ERROR: missing ID map {c['idmap']}")
        parsed[c["name"]]=(parse_gff(c["gff"]), parse_idmap(c["idmap"]), parse_rmout(c["rmout"], families))
    ogmap=ogs
    sample_ogmaps={}
    target_indexes={}
    all_configs=configs + [{"name":r["name"]} for r in targets]
    for c in all_configs:
        name=c["name"]; g,ids,_=parsed[name]
        map_path=(args.gene_map_dir / f"{name}.gene_to_orthogroup.tsv") if args.gene_map_dir else None
        tm=parse_gene_map(map_path) if map_path else {}
        if not tm:
            for seq, genes in g.items():
                for gene in genes:
                    mapped=set()
                    for alias in aliases({seq: [gene]}):
                        if alias in ids: mapped.update(ogs.get(ids[alias],set()))
                    tm[gene["id"]]=mapped
        sample_ogmaps[name]=tm
        if name in {r["name"] for r in targets}:
            target_indexes[name]=build_target_index(g, tm, parsed[name][2])
    rows=[]; anchor_rows=[]
    for q in configs:
        qg, qids, qrm=parsed[q["name"]]; q_ogmap=sample_ogmaps[q["name"]]
        for fam, loci in qrm.items():
            print(f"[geneanchor] query={q['name']} family={fam} loci={len(loci)}", flush=True)
            for seq,start,end in loci:
                genes=qg.get(seq,[]); up,dn=nearest(genes,start,end)
                up_ogs=q_ogmap.get(up["id"],set()) if up else set(); dn_ogs=q_ogmap.get(dn["id"],set()) if dn else set()
                base={"query":q["name"],"pan_family":fam,"query_seqid":seq,"query_start":start,"query_end":end,"upstream_gene":up["id"] if up else "","downstream_gene":dn["id"] if dn else "","upstream_orthogroups":",".join(sorted(up_ogs)),"downstream_orthogroups":",".join(sorted(dn_ogs))}
                for t in targets:
                    hits, anchor_pairs=target_interval(target_indexes[t["name"]],up_ogs,dn_ogs,fam,args.max_gap)
                    status="shared" if hits else ("absent" if anchor_pairs else "anchor_missing")
                    r=dict(base); r.update({"target":t["name"],"target_group":t.get("group", ""),"target_hit_count":len(hits),"target_anchor_pair_count":anchor_pairs,"target_status":status,"target_hit_coordinates":";".join(f"{u['seqid']}:{a}-{b}" for u,d,a,b in hits[:10])}); rows.append(r)
                anchor_rows.append(base)
    fields=["query","pan_family","query_seqid","query_start","query_end","upstream_gene","downstream_gene","upstream_orthogroups","downstream_orthogroups","target","target_group","target_hit_count","target_anchor_pair_count","target_status","target_hit_coordinates"]
    write_tsv(args.output_dir/"geneanchor_locus_evidence.tsv",fields,rows)
    write_tsv(args.output_dir/"geneanchor_anchor_catalog.tsv",list(anchor_rows[0]) if anchor_rows else ["query","pan_family"],anchor_rows)
    summary=[]
    by=defaultdict(list)
    for r in rows: by[(r["query"],r["pan_family"],r["target"])].append(r)
    for (q,f,t),rs in by.items():
        counts=defaultdict(int)
        for r in rs: counts[r["target_status"]]+=1
        summary.append({"query":q,"pan_family":f,"target":t,"loci":len(rs),"shared":counts["shared"],"absent":counts["absent"],"anchor_missing":counts["anchor_missing"]})
    write_tsv(args.output_dir/"geneanchor_family_summary.tsv",["query","pan_family","target","loci","shared","absent","anchor_missing"],summary)
    calls=[]
    for q in configs:
        qname=q["name"]
        qgroup="A" if "A" in qname else "B"
        counterpart="pziA2" if qname=="pziA1" else "pziB2" if qname=="pziB1" else ""
        psp_names=["pspA1","pspA2"] if qgroup=="A" else ["pspB1","pspB2"]
        for fam in sorted(families):
            locus_count=sum(1 for r in rows if r["query"]==qname and r["pan_family"]==fam and r["target"]==qname)
            if not locus_count: continue
            states={t: [r["target_status"] for r in rows if r["query"]==qname and r["pan_family"]==fam and r["target"]==t] for t in [counterpart]+psp_names if t}
            def all_state(t, state): return bool(states.get(t)) and all(x==state for x in states[t])
            pzi_fixed=all_state(counterpart,"shared") if counterpart else False
            psp_absent=all(all_state(t,"absent") for t in psp_names)
            psp_resolved=all(bool(states.get(t)) and all(x in ("absent", "shared") for x in states[t]) for t in psp_names)
            representative_psp=[t for t in psp_names if states.get(t)]
            if not counterpart and len(representative_psp)==1 and all_state(representative_psp[0],"absent"):
                call="PZI_candidate_representative_only"
            elif pzi_fixed and psp_absent and psp_resolved: call="fixed_PZI_specific"
            elif any(any(x=="shared" for x in states.get(t,[])) for t in psp_names): call="shared_with_PSP"
            elif not psp_resolved or not all_state(counterpart,"shared"): call="unresolved_or_polymorphic"
            else: call="PZI_enriched"
            calls.append({"query":qname,"subgenome":qgroup,"pan_family":fam,"query_locus_count":locus_count,"pzi_counterpart_status":";".join(states.get(counterpart,[])),"psp_status":";".join(f"{t}:{','.join(states.get(t,[]))}" for t in psp_names),"species_call":call})
    write_tsv(args.output_dir/"geneanchor_species_calls.tsv",["query","subgenome","pan_family","query_locus_count","pzi_counterpart_status","psp_status","species_call"],calls)
    qc=[]
    for c in configs + [{"name":r["name"]} for r in targets]:
        g, ids, rm = parsed[c["name"]]
        gene_count=sum(len(v) for v in g.values())
        mapped=sum(1 for seq, genes in g.items() for gene in genes if sample_ogmaps[c["name"]].get(gene["id"]))
        qc.append({"query":c["name"],"family_count":len(families),"locus_count":sum(len(v) for v in rm.values()),"protein_coding_gene_count":gene_count,"genes_with_orthogroup":mapped,"orthogroup_mapping_fraction":f"{mapped/gene_count:.6f}" if gene_count else "0"})
    write_tsv(args.output_dir/"geneanchor_qc.tsv",["query","family_count","locus_count","protein_coding_gene_count","genes_with_orthogroup","orthogroup_mapping_fraction"],qc)
if __name__=="__main__": main()
