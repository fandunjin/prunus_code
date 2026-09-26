#!/usr/bin/env python3
import argparse
import hashlib
from pathlib import Path


def read_fasta(path):
    name = None
    header = None
    chunks = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, header, "".join(chunks)
                header = line[1:]
                name = header.split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, header, "".join(chunks)


def wrap(seq, width=80):
    for i in range(0, len(seq), width):
        yield seq[i:i + width]


def gap_runs(seq):
    start = None
    for i, base in enumerate(seq, start=1):
        if base in "Nn":
            if start is None:
                start = i
        elif start is not None:
            yield start, i - 1
            start = None
    if start is not None:
        yield start, len(seq)


def agp_rows(name, seq):
    rows = []
    pos = 1
    part = 1
    component_index = 1
    for gap_start, gap_end in list(gap_runs(seq)) + [(len(seq) + 1, len(seq))]:
        if pos <= gap_start - 1:
            comp_start = pos
            comp_end = gap_start - 1
            comp_len = comp_end - comp_start + 1
            comp_id = f"{name}_part{component_index:06d}"
            rows.append(
                [name, comp_start, comp_end, part, "W", comp_id, 1, comp_len, "+"]
            )
            part += 1
            component_index += 1
        if gap_start <= gap_end:
            gap_len = gap_end - gap_start + 1
            rows.append(
                [name, gap_start, gap_end, part, "N", gap_len, "scaffold", "yes", "paired-ends"]
            )
            part += 1
        pos = gap_end + 1
    return rows


def parse_args():
    parser = argparse.ArgumentParser(
        description="Finalize a Juicebox-adjusted 2n chromosome FASTA."
    )
    parser.add_argument("--input-fasta", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--chr-count", type=int, required=True)
    parser.add_argument("--pair-mode", default="adjacent")
    parser.add_argument("--representative-copy", type=int, choices=[1, 2], default=1)
    return parser.parse_args()


def main():
    args = parse_args()
    input_fasta = Path(args.input_fasta)
    out_dir = Path(args.out_dir)
    prefix = args.prefix
    chr_count = args.chr_count
    pair_mode = args.pair_mode
    representative_copy = args.representative_copy

    if pair_mode != "adjacent":
        raise SystemExit("ERROR: only PAIR_MODE=adjacent is implemented.")
    if chr_count % 2 != 0:
        raise SystemExit("ERROR: CHR_COUNT must be even for adjacent homolog pairing.")

    final_fasta = out_dir / f"{prefix}.final.fa"
    chr_fasta = out_dir / f"{prefix}.chromosomes.fa"
    unplaced_fasta = out_dir / f"{prefix}.unplaced.fa"
    rep_fasta = out_dir / f"{prefix}.n16_representative_copy{representative_copy}.fa"
    name_map = out_dir / f"{prefix}.rename_map.tsv"
    homolog_pairs = out_dir / f"{prefix}.homolog_pairs.tsv"
    chrom_sizes = out_dir / f"{prefix}.chrom.sizes"
    chr_sizes = out_dir / f"{prefix}.chromosomes.sizes"
    seq_stats = out_dir / f"{prefix}.sequence_stats.tsv"
    gap_stats = out_dir / f"{prefix}.gap_stats.tsv"
    agp_file = out_dir / f"{prefix}.agp"
    md5_file = out_dir / f"{prefix}.md5"

    records = list(read_fasta(input_fasta))
    if len(records) < chr_count:
        raise SystemExit(f"ERROR: found only {len(records)} records, expected at least {chr_count}.")

    seen_old = set()
    renamed = []
    for idx, (old_name, old_header, seq) in enumerate(records, start=1):
        if old_name in seen_old:
            raise SystemExit(f"ERROR: duplicate input FASTA name: {old_name}")
        seen_old.add(old_name)
        if not seq:
            raise SystemExit(f"ERROR: empty sequence in input FASTA: {old_name}")

        if idx <= chr_count:
            new_name = f"chr{idx:02d}"
            category = "chromosome"
            pair_id = (idx + 1) // 2
            copy = 1 if idx % 2 == 1 else 2
        else:
            new_name = f"unplaced{idx - chr_count:03d}"
            category = "unplaced"
            pair_id = ""
            copy = ""
        renamed.append(
            {
                "old_name": old_name,
                "old_header": old_header,
                "new_name": new_name,
                "category": category,
                "pair_id": pair_id,
                "copy": copy,
                "seq": seq,
            }
        )

    with final_fasta.open("w") as final_out, chr_fasta.open("w") as chr_out, unplaced_fasta.open("w") as unplaced_out:
        for rec in renamed:
            target_handles = [final_out]
            if rec["category"] == "chromosome":
                target_handles.append(chr_out)
            else:
                target_handles.append(unplaced_out)
            for handle in target_handles:
                handle.write(f">{rec['new_name']}\n")
                for line in wrap(rec["seq"]):
                    handle.write(line + "\n")

    with rep_fasta.open("w") as out:
        rep_index = 1
        for rec in renamed:
            if rec["category"] == "chromosome" and rec["copy"] == representative_copy:
                out.write(f">chr{rep_index:02d}\n")
                for line in wrap(rec["seq"]):
                    out.write(line + "\n")
                rep_index += 1
        if rep_index != chr_count // 2 + 1:
            raise SystemExit("ERROR: representative n=16 FASTA did not receive the expected record count.")

    with name_map.open("w") as out:
        out.write("old_name\tnew_name\tcategory\tlength\tpair_id\tcopy\told_header\n")
        for rec in renamed:
            out.write(
                f"{rec['old_name']}\t{rec['new_name']}\t{rec['category']}\t"
                f"{len(rec['seq'])}\t{rec['pair_id']}\t{rec['copy']}\t{rec['old_header']}\n"
            )

    with homolog_pairs.open("w") as out:
        out.write("pair_id\tcopy1_chr\tcopy2_chr\tcopy1_original\tcopy2_original\trepresentative_chr\n")
        chroms = [rec for rec in renamed if rec["category"] == "chromosome"]
        for i in range(0, len(chroms), 2):
            copy1 = chroms[i]
            copy2 = chroms[i + 1]
            pair_id = i // 2 + 1
            representative = copy1 if representative_copy == 1 else copy2
            out.write(
                f"pair{pair_id:02d}\t{copy1['new_name']}\t{copy2['new_name']}\t"
                f"{copy1['old_name']}\t{copy2['old_name']}\t{representative['new_name']}\n"
            )

    with chrom_sizes.open("w") as all_out, chr_sizes.open("w") as chr_out:
        for rec in renamed:
            line = f"{rec['new_name']}\t{len(rec['seq'])}\n"
            all_out.write(line)
            if rec["category"] == "chromosome":
                chr_out.write(line)

    with seq_stats.open("w") as out:
        out.write("seq_id\tcategory\tlength\tn_bases\tn_percent\tgap_count\tlongest_gap\tnon_n_bases\n")
        for rec in renamed:
            gaps = [(start, end) for start, end in gap_runs(rec["seq"])]
            n_bases = sum(end - start + 1 for start, end in gaps)
            longest_gap = max((end - start + 1 for start, end in gaps), default=0)
            length = len(rec["seq"])
            out.write(
                f"{rec['new_name']}\t{rec['category']}\t{length}\t{n_bases}\t"
                f"{(n_bases / length * 100):.6f}\t{len(gaps)}\t{longest_gap}\t{length - n_bases}\n"
            )

    with gap_stats.open("w") as out:
        out.write("seq_id\tgap_start\tgap_end\tgap_length\n")
        for rec in renamed:
            for start, end in gap_runs(rec["seq"]):
                out.write(f"{rec['new_name']}\t{start}\t{end}\t{end - start + 1}\n")

    with agp_file.open("w") as out:
        out.write("##agp-version 2.1\n")
        for rec in renamed:
            for row in agp_rows(rec["new_name"], rec["seq"]):
                out.write("\t".join(map(str, row)) + "\n")

    with md5_file.open("w") as out:
        for path in [final_fasta, chr_fasta, unplaced_fasta, rep_fasta]:
            digest = hashlib.md5(path.read_bytes()).hexdigest()
            out.write(f"{digest}  {path.name}\n")

    print(f"Wrote final FASTA: {final_fasta}")
    print(f"Wrote chromosome FASTA: {chr_fasta}")
    print(f"Wrote representative n=16 FASTA: {rep_fasta}")
    print(f"Wrote rename map: {name_map}")
    print(f"Wrote AGP: {agp_file}")


if __name__ == "__main__":
    main()
