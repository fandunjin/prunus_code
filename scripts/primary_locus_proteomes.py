#!/usr/bin/env python3
"""Select one valid representative per biological locus from audited frozen inputs.

Manifest TSV: track, frozen_pep, id_map, restored_pep, cds (optional).
restored_pep must have exactly the same ordered valid protein sequences as
frozen_pep; this restores provenance without collapsing identical paralogs.
CDS lengths come from cds_len= headers or an explicit CDS FASTA, never from
an amino-acid-length approximation. Longest CDS wins, then protein length,
then original identifier. Original stable IDs are retained for annotation joins.
"""
import argparse, csv, gzip, hashlib, json, re
from collections import defaultdict
from pathlib import Path


def fasta(path):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt", encoding="utf-8") as f:
        h = None
        s = []
        for l in f:
            if l.startswith(">"):
                if h is not None:
                    yield h, "".join(s).upper().rstrip("*")
                h = l[1:].strip()
                s = []
            else:
                s.append(l.strip())
        if h is not None:
            yield h, "".join(s).upper().rstrip("*")


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def rows(p):
    with open(p, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write(p, rs, fields):
    with open(p, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rs)


def locus(header):
    for pat in [
        r"(?:^|\s)gene=([^\s]+)",
        r"(?:^|\s)locus=([^\s]+)",
        r"^(gene-[^|\s]+)\|",
    ]:
        m = re.search(pat, header)
        if m:
            return m.group(1)
    token = header.split()[0]
    if re.fullmatch(r"Gtr\S+\.t\d+", token):
        return re.sub(r"\.t\d+$", "", token)
    raise ValueError("No verified parent-locus identity: " + header)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--min-cds", type=int, default=300)
    a = ap.parse_args()
    if a.out.exists():
        raise SystemExit("Output exists; use a new immutable run directory")
    (a.out / "pep").mkdir(parents=True)
    (a.out / "id_maps").mkdir()
    report = []
    for spec in rows(a.manifest):
        tr = spec["track"]
        fr = list(fasta(spec["frozen_pep"]))
        mp = rows(spec["id_map"])
        assert len(mp) == len(fr) and all(
            x["stable_id"] == y[0].split()[0] for x, y in zip(mp, fr)
        ), tr
        restored = spec.get("restored_pep", "")
        if restored:
            src = [(h, s) for h, s in fasta(restored) if "*" not in s]
            assert len(src) == len(fr) and all(
                x[1] == y[1] for x, y in zip(fr, src)
            ), f"{tr}: restored sequences differ"
        else:
            src = [(m["original_header"], s) for m, (_, s) in zip(mp, fr)]
        cds = {}
        if spec.get("cds"):
            for h, s in fasta(spec["cds"]):
                key = h.split()[0]
                assert key not in cds, (tr, key)
                cds[key] = len(s)
        by = defaultdict(list)
        allrows = []
        for m, (h, s) in zip(mp, src):
            assert "*" not in s and len(s) == int(m["protein_length"]), tr
            orig = h.split()[0]
            g = locus(h)
            mt = re.search(r"(?:^|\s)cds_len=(\d+)", h)
            n = (
                int(mt.group(1))
                if mt
                else cds.get(orig, cds.get(re.sub(r"\.p$", "", orig)))
            )
            if n is None:
                raise ValueError(f"{tr}: missing CDS length for {orig}")
            r = dict(
                stable_id=m["stable_id"],
                original_id=orig,
                original_header=h,
                protein_length=len(s),
                gene_locus=g,
                cds_length=n,
                selected=False,
                reason="short_cds" if n < a.min_cds else "alternative_isoform",
            )
            allrows.append(r)
            if n >= a.min_cds:
                by[g].append((r, s))
        chosen = []
        for g, items in by.items():
            r, s = min(
                items,
                key=lambda x: (
                    -x[0]["cds_length"],
                    -x[0]["protein_length"],
                    x[0]["original_id"],
                ),
            )
            r["selected"] = True
            r["reason"] = "representative"
            chosen.append((r, s))
        chosen.sort(key=lambda x: x[0]["stable_id"])
        assert len(chosen) == len({r["gene_locus"] for r, s in chosen})
        pep = a.out / f"pep/{tr}.fa"
        with pep.open("w", encoding="ascii") as f:
            for r, s in chosen:
                f.write(">" + r["stable_id"] + "\n")
                for i in range(0, len(s), 60):
                    f.write(s[i : i + 60] + "\n")
        fields = list(allrows[0])
        write(a.out / f"id_maps/{tr}.id_map.tsv", [r for r, s in chosen], fields)
        write(a.out / f"{tr}.selection_audit.tsv", allrows, fields)
        rr = dict(
            track=tr,
            input_records=len(fr),
            input_loci=len({r["gene_locus"] for r in allrows}),
            representatives=len(chosen),
            short_records=sum(r["reason"] == "short_cds" for r in allrows),
            alternative_records=sum(
                r["reason"] == "alternative_isoform" for r in allrows
            ),
            frozen_sha256=sha(spec["frozen_pep"]),
            restored_sha256=sha(restored) if restored else "",
            primary_sha256=sha(pep),
        )
        report.append(rr)
        print(json.dumps(rr), flush=True)
    write(a.out / "inventory.tsv", report, list(report[0]))
    (a.out / "contract.json").write_text(
        json.dumps(
            {
                "manifest_sha256": sha(a.manifest),
                "min_cds": a.min_cds,
                "selection": "longest valid CDS; protein length; original ID",
                "unit": "parent gene locus",
                "identical_paralogs": "kept as distinct loci",
            },
            indent=2,
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
