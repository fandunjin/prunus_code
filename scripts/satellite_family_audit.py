"""Audit every 163-170 bp monomer against published Prunus_170 references.
References are doubled to allow a tandem monomer's arbitrary phase. Exact
sequence copies are queried once, but all locus occurrences remain in counts.
"""

from pathlib import Path
import csv, json, subprocess, hashlib
from collections import Counter
import argparse

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--source", type=Path, required=True)
ap.add_argument("--out", type=Path, required=True)
args = ap.parse_args()
S = args.source
R = args.out
R.mkdir(exist_ok=True, parents=True)
targets = [
    "pziA1",
    "pziA2",
    "pziB1",
    "pziB2",
    "pspA1",
    "pspA2",
    "pspB1",
    "pspB2",
    "ppaA",
    "ppaB",
    "pcam",
    "ppe",
    "gtr",
    "rch",
]


def file_sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


contract = {
    "source": str(S.resolve()),
    "reference_sha256": {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in [
            S / "analysis/published_monomers.fa",
            S / "analysis/published_monomers_wang2022.fa",
        ]
    },
    "monomer_inputs_sha256": {
        t: file_sha(S / f"trash/{t}/{t}.fa_repeats_with_seq.csv") for t in targets
    },
    "width": [163, 170],
    "identity": 70,
    "query_span": 100,
    "evalue": 1e-5,
}
cp = R / "input_contract.json"
if cp.exists() and json.loads(cp.read_text()) != contract:
    raise SystemExit("Resume contract differs; choose a new output directory")
cp.write_text(json.dumps(contract, indent=2))


def fasta(p):
    h = None
    seq = []
    for l in p.open():
        if l.startswith(">"):
            if h:
                yield h, "".join(seq).upper()
            h = l[1:].split()[0]
            seq = []
        else:
            seq.append(l.strip())
    if h:
        yield h, "".join(seq).upper()


refs = list(fasta(S / "analysis/published_monomers.fa")) + list(
    fasta(S / "analysis/published_monomers_wang2022.fa")
)
with (R / "published_circular_references.fa").open("w") as f:
    seen = set()
    for h, s in refs:
        if s in seen:
            continue
        seen.add(s)
        f.write(">" + h + "\n" + s + s + "\n")
subprocess.run(
    [
        "makeblastdb",
        "-in",
        str(R / "published_circular_references.fa"),
        "-dbtype",
        "nucl",
        "-out",
        str(R / "references"),
    ],
    check=True,
    stdout=subprocess.DEVNULL,
)
summaries = (
    list(csv.DictReader((R / "census.tsv").open(), delimiter="\t"))
    if (R / "census.tsv").exists()
    else []
)
completed = {r["target"] for r in summaries}
for target in targets:
    if target in completed:
        continue
    dest = R / target
    dest.mkdir(exist_ok=True)
    src = S / f"trash/{target}/{target}.fa_repeats_with_seq.csv"
    seqs = Counter()
    allbp = 0
    alln = 0
    for r in csv.DictReader(src.open()):
        width = int(float(r["width"]))
        s = r["sequence"].upper()
        n = abs(int(float(r["end"])) - int(float(r["start"]))) + 1
        if width < 50:
            continue  # short boundary fragments are outside the declared census
        lo = min(int(float(r["start"])), int(float(r["end"])))
        hi = max(int(float(r["start"])), int(float(r["end"])))
        n = hi - max(1, lo) + 1
        assert len(s) == n and n <= width, (target, r)
        width = n
        if width < 50:
            continue
        if width >= 50:
            alln += 1
            allbp += n
        if 163 <= width <= 170:
            seqs[s] += 1
    with (dest / "queries.fa").open("w") as f:
        for i, s in enumerate(seqs):
            f.write(f">q{i}\n{s}\n")
    if seqs:
        cmd = [
            "blastn",
            "-task",
            "blastn",
            "-query",
            str(dest / "queries.fa"),
            "-db",
            str(R / "references"),
            "-out",
            str(dest / "hits.tsv"),
            "-outfmt",
            "6 qseqid sseqid pident length qstart qend evalue bitscore",
            "-word_size",
            "7",
            "-evalue",
            "1e-5",
            "-perc_identity",
            "70",
            "-dust",
            "yes",
            "-num_threads",
            "4",
            "-max_target_seqs",
            "20",
        ]
        with (dest / "blast.log").open("w") as f:
            subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.STDOUT)
        accepted = {
            int(r[0][1:])
            for r in csv.reader((dest / "hits.tsv").open(), delimiter="\t")
            if float(r[2]) >= 70
            and int(r[3]) >= 100
            and abs(int(r[5]) - int(r[4])) + 1 >= 100
        }
    else:
        accepted = set()
    famseq = {s for i, s in enumerate(seqs) if i in accepted}
    with (dest / "sequence_classification.tsv").open("w") as f:
        wr = csv.writer(f, delimiter="\t")
        wr.writerow(["sequence", "occurrences", "width", "published_homology"])
        for s, n in seqs.items():
            wr.writerow([s, n, len(s), s in famseq])
    famn = sum(seqs[s] for s in famseq)
    fambp = sum(len(s) * seqs[s] for s in famseq)
    candbp = sum(len(s) * n for s, n in seqs.items())
    with (dest / "monomer_intervals.tsv").open("w") as f:
        wr = csv.writer(f, delimiter="\t")
        wr.writerow(
            [
                "chromosome",
                "start0",
                "end",
                "array_id",
                "width",
                "class_token",
                "published_homology",
            ]
        )
        for r in csv.DictReader(src.open()):
            lo = max(1, min(int(float(r["start"])), int(float(r["end"]))))
            hi = max(int(float(r["start"])), int(float(r["end"])))
            n = hi - lo + 1
            if n < 50:
                continue
            wr.writerow(
                [
                    r["seqID"],
                    lo - 1,
                    hi,
                    r["arrayID"],
                    n,
                    r["class"],
                    r["sequence"].upper() in famseq,
                ]
            )
    rs = dict(
        target=target,
        all_monomer_records_ge50=alln,
        all_monomer_bp_ge50=allbp,
        width_candidate_unique_sequences=len(seqs),
        width_candidate_bp=candbp,
        homology_monomer_records=famn,
        homology_monomer_bp=fambp,
        width_candidates_without_supported_homology_bp=candbp - fambp,
    )
    summaries.append(rs)
    with (R / "census.tsv").open("w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rs), delimiter="\t")
        wr.writeheader()
        wr.writerows(summaries)
    print(json.dumps(rs), flush=True)
(R / "DONE.json").write_text(
    json.dumps(
        {
            "targets": len(summaries),
            "threshold": "BLASTN identity >=70%, alignment >=100 columns AND >=100 query bases, E<=1e-5, published references only, doubled for cyclic phase",
            "unit": "monomer occurrence; coordinate intervals are 0-based half-open; gaps/unmatched monomers not counted as family",
        },
        indent=2,
    )
)
