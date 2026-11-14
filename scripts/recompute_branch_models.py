#!/usr/bin/env python3
"""Re-fit audited codon alignments with free, gene-specific branch lengths.

Retains original alignments, topology and foreground definitions, but never
fixes arbitrary starting branch lengths. Every fit has an isolated directory.
The same branch-length policy is used for the null, focal and all alternative
foreground models. Reported alignment length is in codons, not nucleotides.
"""
from pathlib import Path
import argparse, csv, re, json, shutil, subprocess, concurrent.futures, time
import numpy as np
from scipy.stats import chi2, wilcoxon
from locus_family_analysis import bh
from primary_locus_proteomes import sha, write

PARAMS = """noisy = 0
verbose = 0
runmode = 0
seqtype = 1
CodonFreq = 2
clock = 0
aaDist = 0
model = {model}
NSsites = 0
icode = 0
Mgene = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = {omega}
fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 8
getSE = 0
RateAncestor = 0
Small_Diff = .0000005
cleandata = 0
fix_blength = 0
method = 0
"""


def read(p):
    with p.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def parse(p, model):
    t = p.read_text()
    m = re.search(r"lnL\(ntime.*?\):\s*(-?[\d.]+)", t)
    if not m:
        raise ValueError("No likelihood: " + str(p))
    if model == 0:
        om = re.search(r"omega \(dN/dS\) =\s*([\d.]+)", t)
        if not om:
            raise ValueError("No omega: " + str(p))
        bg = fg = float(om[1])
    else:
        om = re.search(r"w \(dN/dS\) for branches:\s*([\d.]+)\s+([\d.]+)", t)
        if not om:
            raise ValueError("No branch omegas: " + str(p))
        bg, fg = map(float, om.groups())
    return {"lnL": float(m[1]), "omega_background": bg, "omega_foreground": fg}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--codeml", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--resume", action="store_true")
    a = ap.parse_args()
    if a.out.exists() and not a.resume:
        raise SystemExit("Output exists; use --resume for identical inputs")
    a.out.mkdir(parents=True, exist_ok=a.resume)
    sets = []
    tasks = []
    manifest = []
    for file, folder, key, group in [
        ("branch_dnds_results.tsv", "runs", "set", None),
        ("branch_dnds_background.tsv", "runs_bg", "og", "Background A"),
        ("branch_dnds_background_B.tsv", "runs_bg_B", "og", "Background B"),
    ]:
        for r in read(a.source / "results" / file):
            name = r[key]
            g = group or (
                "Candidates" if r["candidate"] == "True" else "Family controls"
            )
            source = a.source / folder / name
            d = a.out / g.replace(" ", "_") / name
            d.mkdir(parents=True, exist_ok=a.resume)
            aln = source / "aln.phy"
            header = aln.read_text().splitlines()[0].split()
            assert int(header[0]) == 6 and int(header[1]) % 3 == 0
            shutil.copyfile(aln, d / "aln.phy")
            sets.append(
                dict(
                    group=g,
                    name=name,
                    directory=d,
                    codons=int(header[1]) // 3,
                    family=r.get("family", "Single-copy orthogroup"),
                    subfamily=r.get("subfamily", ""),
                )
            )
            manifest.append(
                dict(
                    group=g,
                    gene_set=name,
                    alignment_sha256=sha(aln),
                    nucleotides=int(header[1]),
                    codons=int(header[1]) // 3,
                )
            )
            for tag in ["null", "branch"] + (
                [f"perm{i}" for i in range(16)] if group else []
            ):
                tree = source / f"tree_{tag}.nwk"
                assert tree.exists(), tree
                td = d / tag
                td.mkdir(exist_ok=a.resume)
                # Codeml needs an explicit tree count. Suppress the degree-two
                # root for a non-clock, reversible branch model.
                raw = tree.read_text().strip()
                m = re.fullmatch(
                    r"\(gtr:[0-9.eE+-]+(?P<label> #1)?, \((?P<rest>.*)\):[0-9.eE+-]+(?: #1)?\);",
                    raw,
                )
                if not m:
                    raise ValueError("Unexpected source tree: " + raw)
                unrooted = "(" + m["rest"] + ",gtr:0.5" + (m["label"] or "") + ");"
                (td / "tree.nwk").write_text("6 1\n" + unrooted + "\n")
                tasks.append(
                    (g, name, tag, td, d / "aln.phy", 0 if tag == "null" else 2)
                )
    write(a.out / "input_manifest.tsv", manifest, list(manifest[0]))

    def fit(task):
        g, name, tag, d, aln, model = task

        def once(omega, attempt):
            ctl = d / f"fit{attempt}.ctl"
            outfile = d / f"mlc{attempt}"
            if a.resume and outfile.exists():
                try:
                    return parse(outfile, model)
                except ValueError:
                    pass
            ctl.write_text(
                f"seqfile = {aln}\ntreefile = {d}/tree.nwk\noutfile = {outfile}\n"
                + PARAMS.format(model=model, omega=omega)
            )
            with (d / f"fit{attempt}.log").open("w") as f:
                p = subprocess.run(
                    [str(a.codeml), str(ctl)], cwd=d, stdout=f, stderr=subprocess.STDOUT
                )
            if p.returncode:
                raise RuntimeError(f"codeml exit {p.returncode}: {d}")
            return parse(outfile, model)

        value = once(0.4, 1)
        if tag in ["null", "branch"]:
            value = max([value, once(1.2, 2)], key=lambda x: x["lnL"])
        return (g, name, tag), value

    fits = {}
    restart_count = 0
    start = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
        for i, (key, value) in enumerate(pool.map(fit, tasks), 1):
            fits[key] = value
            if i % 50 == 0:
                (a.out / "progress.json").write_text(
                    json.dumps(
                        {
                            "finished": i,
                            "total": len(tasks),
                            "elapsed_seconds": round(time.time() - start),
                        }
                    )
                )
                print(i, "/", len(tasks), flush=True)
    out = []
    perms = []
    for s in sets:
        key = (s["group"], s["name"])
        null = fits[(*key, "null")]
        alt = fits[(*key, "branch")]
        if alt["lnL"] < null["lnL"] - 0.0001:
            raise ValueError(f"Nested likelihood inversion: {key}")
        lrt = max(0, 2 * (alt["lnL"] - null["lnL"]))
        bg = alt["omega_background"]
        fg = alt["omega_foreground"]
        if bg <= 0:
            raise ValueError(f"Non-positive omega denominator: {key}")
        out.append(
            dict(
                Group=s["group"],
                Gene_set=s["name"],
                Family=s["family"],
                Subfamily=s["subfamily"],
                Codon_length=s["codons"],
                Omega_null=null["omega_background"],
                Omega_foreground=fg,
                Omega_background=bg,
                Omega_ratio=fg / bg,
                Log2_omega_ratio=float(np.log2(fg / bg)),
                LogL_null=null["lnL"],
                LogL_branch=alt["lnL"],
                LRT=lrt,
                P=float(chi2.sf(lrt, 1)),
            )
        )
        if s["group"].startswith("Background"):
            for i in range(16):
                f = fits[(*key, f"perm{i}")]
                if f["lnL"] < null["lnL"] - 0.0001:
                    td = s["directory"] / f"perm{i}"
                    alternatives = [f]
                    for j, omega in enumerate([null["omega_background"], 1.2, 2.0], 3):
                        ctl = td / f"fit{j}.ctl"
                        outfile = td / f"mlc{j}"
                        ctl.write_text(
                            f"seqfile = {s['directory']}/aln.phy\ntreefile = {td}/tree.nwk\noutfile = {outfile}\n"
                            + PARAMS.format(model=2, omega=omega)
                        )
                        with (td / f"fit{j}.log").open("w") as log:
                            pr = subprocess.run(
                                [str(a.codeml), str(ctl)],
                                cwd=td,
                                stdout=log,
                                stderr=subprocess.STDOUT,
                            )
                        if pr.returncode:
                            raise RuntimeError("codeml restart failed")
                        alternatives.append(parse(outfile, 2))
                    f = max(alternatives, key=lambda x: x["lnL"])
                    restart_count += 1
                    if f["lnL"] < null["lnL"] - 0.0001:
                        raise ValueError(
                            f"Persistent nested likelihood inversion: {key}, {i}"
                        )
                perms.append(
                    dict(
                        Subgenome=s["group"][-1],
                        Gene_set=s["name"],
                        Configuration=i,
                        Omega_foreground=f["omega_foreground"],
                        Omega_background=f["omega_background"],
                        Omega_ratio=f["omega_foreground"] / f["omega_background"],
                        LogL=f["lnL"],
                    )
                )
    for r, q in zip(out, bh([r["P"] for r in out])):
        r["Q"] = q
    write(a.out / "selection_gene_sets.tsv", out, list(out[0]))
    write(a.out / "alternative_foregrounds.tsv", perms, list(perms[0]))
    summary = {
        "restarted_alternative_foregrounds": restart_count,
        "sets": len(out),
        "fits": len(fits),
        "branch_lengths": "estimated separately for every gene set under every model",
        "groups": {},
    }
    for g in dict.fromkeys(r["Group"] for r in out):
        rs = [r for r in out if r["Group"] == g]
        summary["groups"][g] = {
            "n": len(rs),
            "significant_q05": sum(r["Q"] < 0.05 for r in rs),
            "median_ratio": float(np.median([r["Omega_ratio"] for r in rs])),
            "mean_ratio": float(np.mean([r["Omega_ratio"] for r in rs])),
        }
    by = {
        g: {r["Gene_set"]: r for r in out if r["Group"] == g}
        for g in ["Background A", "Background B"]
    }
    common = sorted(set(by["Background A"]) & set(by["Background B"]))
    stat = wilcoxon(
        [by["Background A"][k]["Log2_omega_ratio"] for k in common],
        [by["Background B"][k]["Log2_omega_ratio"] for k in common],
        alternative="two-sided",
        method="exact",
    )
    summary["paired_test"] = {
        "n": len(common),
        "W": float(stat.statistic),
        "P": float(stat.pvalue),
    }
    summary["alternative_foreground_mean_ratios"] = {
        g: {
            i: float(
                np.mean(
                    [
                        r["Omega_ratio"]
                        for r in perms
                        if r["Subgenome"] == g and r["Configuration"] == i
                    ]
                )
            )
            for i in range(16)
        }
        for g in ["A", "B"]
    }
    (a.out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
