#!/usr/bin/env python3
"""Full-precision family scan from one-locus OrthoFinder output.

Counts and denominators use assigned parent loci. All 17 tracks are unique;
overlapping biological groupings are used only to form contrasts. Each
selected protein inherits its own existing Swiss-Prot annotation, and TE
exclusion searches every member's complete description. GO tests use all
background terms represented in >=3 tested, GO-annotated families; the
foreground count threshold is a reporting filter, not the BH universe.
"""
import argparse, json, math, re
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.stats import binomtest, fisher_exact
from primary_locus_proteomes import rows, write, sha

RAC = ["pzi_A", "pzi_B", "psp_A", "psp_B", "ppa_A", "ppa_B"]
COR = ["pcam", "pav", "pspe", "pmal"]
SOL = ["par", "pmu", "pdu", "psal", "ppe"]
OUT = ["gtr", "rch"]
TRACKS = RAC + COR + SOL + OUT
GROUPS = {
    "rac": RAC,
    "racA": RAC[::2],
    "racB": RAC[1::2],
    "pzi": RAC[:2],
    "psp": RAC[2:4],
    "ppa": RAC[4:6],
    "cor": COR,
    "sol": SOL,
    "out": OUT,
}
TE = re.compile(
    r"transposon|retrotransposon|retrovirus|reverse transcriptase|integrase|transposase|gag.polyprotein|copia|gypsy|pol polyprotein",
    re.I,
)


def bh(p):
    p = np.asarray(p, float)
    if len(p) == 0:
        return []
    assert np.isfinite(p).all() and ((p >= 0) & (p <= 1)).all()
    order = np.argsort(p, kind="stable")
    q = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[
        ::-1
    ]
    ans = np.empty(len(p))
    ans[order] = np.minimum(q, 1)
    return ans.tolist()


def contrast(obs, n, base, m):
    assert 0 <= obs <= n and 0 <= base <= m and n > 0 and m > 0
    expected = n * base / m
    return (
        expected,
        (obs / expected if expected else math.nan),
        binomtest(obs, n, base / m, alternative="greater").pvalue,
        binomtest(obs, n, base / m, alternative="less").pvalue,
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--orthogroups", type=Path, required=True)
    ap.add_argument("--primary-inputs", type=Path, required=True)
    ap.add_argument("--annotations", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    if a.out.exists():
        raise SystemExit("Use a new output directory")
    a.out.mkdir(parents=True)
    maps = {
        tr: {
            r["stable_id"]: r
            for r in rows(a.primary_inputs / f"id_maps/{tr}.id_map.tsv")
        }
        for tr in TRACKS
    }
    for tr, m in maps.items():
        assert len(m) == len({r["gene_locus"] for r in m.values()}), tr
    ann = {}
    godesc = {}
    for sp in dict.fromkeys([x.split("_")[0] for x in TRACKS]):
        p = a.annotations / sp / "swissprot_functional_annotation.tsv"
        ann[sp] = {r["protein_id"]: r for r in rows(p)}
        gp = a.annotations / sp / "swissprot_go.tsv"
        for r in rows(gp):
            godesc.setdefault(r["go_id"], r["go_description"])
    counts = rows(a.orthogroups / "Orthogroups.GeneCount.tsv")
    members = rows(a.orthogroups / "Orthogroups.tsv")
    members = {r["Orthogroup"]: r for r in members}
    tot = {tr: sum(int(r[tr]) for r in counts) for tr in TRACKS}
    totals = {g: sum(tot[t] for t in ts) for g, ts in GROUPS.items()}
    result = []
    allann = []
    seen = set()
    for row in counts:
        raw_og = row["Orthogroup"]
        og = "LGF" + raw_og.removeprefix("OG")
        d = {tr: int(row[tr]) for tr in TRACKS}
        mr = members[raw_og]
        gos = set()
        symbols = Counter()
        rsymbols = Counter()
        descs = []
        rdescs = []
        te = False
        src = set()
        for tr in TRACKS:
            ids = [x.strip() for x in mr[tr].split(",") if x.strip()]
            assert len(ids) == d[tr], (og, tr)
            for sid in ids:
                assert sid in maps[tr] and (tr, sid) not in seen, (og, tr, sid)
                seen.add((tr, sid))
                m = maps[tr][sid]
                sp = tr.split("_")[0]
                key = m["original_id"] if tr in RAC else sid
                ar = ann[sp].get(key)
                if ar is None:
                    continue
                src.add(tr)
                desc = ar.get("protein_description", "")
                gn = re.split(r"[;\s]+", ar.get("gene_names", "").strip())[0]
                gn = "" if gn.startswith("{") else gn
                if gn:
                    symbols[gn] += 1
                if desc:
                    descs.append(desc.split(";")[0])
                if tr in RAC:
                    if gn:
                        rsymbols[gn] += 1
                    if desc:
                        rdescs.append(desc.split(";")[0])
                te = te or bool(TE.search(gn + " " + desc))
                gos.update(re.findall(r"GO:\d{7}", ar.get("go_terms", "")))
        symbol = (
            (rsymbols or symbols).most_common(1)[0][0] if (rsymbols or symbols) else ""
        )
        description = " | ".join(dict.fromkeys(rdescs or descs))[:250]
        allann.append(
            dict(
                og=og,
                raw_orthofinder_id=raw_og,
                symbol=symbol,
                all_symbols=";".join(sorted(symbols)),
                description=description,
                go_terms=";".join(sorted(gos)),
                te=te,
                source_tracks=";".join(sorted(src)),
            )
        )
        if (
            sum(v > 0 for v in d.values()) < 2
            or sum(d.values()) > 150
            or sum(d[t] for t in COR + SOL + OUT) < 3
        ):
            continue
        r = {
            "og": og,
            "raw_orthofinder_id": raw_og,
            **{g: sum(d[t] for t in ts) for g, ts in GROUPS.items()},
        }
        for suffix, base, base_total in [
            ("cs", r["cor"] + r["sol"], totals["cor"] + totals["sol"]),
            ("cor", r["cor"], totals["cor"]),
            ("sol", r["sol"], totals["sol"]),
            ("out", r["out"], totals["out"]),
        ]:
            exp, en, p, pl = contrast(r["rac"], totals["rac"], base, base_total)
            r.update({f"enr_{suffix}": en, f"p_{suffix}": p, f"p_{suffix}_less": pl})
            if suffix == "cs":
                r["expected_copies"] = exp
        sig = 0
        for sp in ["pzi", "psp", "ppa"]:
            _, en, p, _ = contrast(
                r[sp], totals[sp], r["cor"] + r["sol"], totals["cor"] + totals["sol"]
            )
            r[f"{sp}_p"] = p
            r[f"{sp}_enrichment"] = en
            sig += p < 0.05 and en > 1
        r["n_species_up"] = sig
        for g in ["cor", "sol"]:
            base = sum(r[x] for x in ["rac", "cor", "sol", "out"] if x != g)
            m = sum(totals[x] for x in ["rac", "cor", "sol", "out"] if x != g)
            _, en, p, _ = contrast(r[g], totals[g], base, m)
            r[f"enr_{g}_rest"] = en
            r[f"p_{g}_rest"] = p
        r.update(
            symbol=symbol,
            all_symbols=";".join(sorted(symbols)),
            description=description,
            te=te,
            symbol17=symbol,
            description17=description,
            te17=te,
            go17=";".join(sorted(gos)),
            src_tracks17=";".join(sorted(src)),
        )
        result.append(r)
    assert len(seen) == sum(tot.values()), "Membership/count mismatch"
    for pc, qc in [
        ("p_cs", "q_cs"),
        ("p_cs_less", "q_cs_less"),
        ("p_cor_rest", "q_cor_rest"),
        ("p_sol_rest", "q_sol_rest"),
    ]:
        for r, q in zip(result, bh([r[pc] for r in result])):
            r[qc] = q
    sets = {"racemose": set(), "corymbose": set(), "solitary": set()}
    for r in result:
        plot = r["rac"] > 0 and r["expected_copies"] > 0
        r.update(
            plot_eligible=plot,
            higher=plot and r["enr_cs"] > 1 and r["q_cs"] < 0.05,
            lower=plot and r["enr_cs"] < 1 and r["q_cs_less"] < 0.05,
        )
        r["retained"] = (
            plot
            and r["q_cs"] < 0.05
            and r["enr_cs"] > 1.3
            and r["n_species_up"] >= 2
            and not r["te17"]
        )
        if r["retained"]:
            sets["racemose"].add(r["og"])
        if r["q_cor_rest"] < 0.05 and r["enr_cor_rest"] > 1.3:
            sets["corymbose"].add(r["og"])
        if r["q_sol_rest"] < 0.05 and r["enr_sol_rest"] > 1.3:
            sets["solitary"].add(r["og"])
    sets["cor_sol_shared"] = sets["corymbose"] & sets["solitary"]
    write(a.out / "clade_family_scan.tsv", result, list(result[0]))
    matrix = [
        {
            "Orthogroup": "LGF" + r["Orthogroup"].removeprefix("OG"),
            "raw_orthofinder_id": r["Orthogroup"],
            **{tr: int(r[tr]) for tr in TRACKS},
        }
        for r in counts
    ]
    write(a.out / "all_track_family_counts.tsv", matrix, list(matrix[0]))
    write(a.out / "og_annotation_all17.tsv", allann, list(allann[0]))
    if sets["racemose"]:
        write(
            a.out / "racemose_retained.tsv",
            [r for r in result if r["retained"]],
            list(result[0]),
        )
    famgo = {r["og"]: set(r["go17"].split(";")) for r in result if r["go17"]}
    bg = Counter(g for gs in famgo.values() for g in gs)
    universe = {g: n for g, n in bg.items() if n >= 3}
    gotables = {}
    for name, ids in sets.items():
        fg = ids & set(famgo)
        fc = Counter(g for og in fg for g in famgo[og])
        gr = []
        M = len(famgo)
        N = len(fg)
        for g, K in sorted(universe.items()):
            k = fc[g]
            _, p = fisher_exact(
                [[k, N - k], [K - k, M - N - K + k]], alternative="greater"
            )
            gr.append(
                dict(
                    go_id=g,
                    go_description=godesc.get(g, ""),
                    expanded=k,
                    bg=K,
                    foreground_families=N,
                    background_families=M,
                    p=p,
                )
            )
        for r, q in zip(gr, bh([r["p"] for r in gr])):
            r["q"] = q
            r["reportable"] = r["expanded"] >= 3
        gr.sort(key=lambda r: (r["q"], r["p"], r["go_id"]))
        write(a.out / f"go_{name}.tsv", gr, list(gr[0]))
        gotables[name] = {
            "foreground_with_GO": N,
            "background_with_GO": M,
            "tested_GO_terms": len(gr),
            "significant_reportable_terms": sum(
                r["q"] < 0.05 and r["reportable"] for r in gr
            ),
        }
    stat = dict(
        orthogroups=len(counts),
        assigned_loci_by_track=tot,
        assigned_loci_by_group=totals,
        tested_families=len(result),
        plot_eligible=sum(r["plot_eligible"] for r in result),
        higher=sum(r["higher"] for r in result),
        lower=sum(r["lower"] for r in result),
        retained=len(sets["racemose"]),
        retained_A=sum(r["racA"] for r in result if r["retained"]),
        retained_B=sum(r["racB"] for r in result if r["retained"]),
        clade_set_sizes={k: len(v) for k, v in sets.items()},
        go=gotables,
        source_counts_sha256=sha(a.orthogroups / "Orthogroups.GeneCount.tsv"),
    )
    (a.out / "summary.json").write_text(json.dumps(stat, indent=2))
    print(json.dumps(stat, indent=2))


if __name__ == "__main__":
    main()
