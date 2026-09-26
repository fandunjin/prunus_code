#!/usr/bin/env python3
import argparse
import csv
import re


def read_fasta(path):
    name = None
    header = None
    chunks = []
    with open(path) as handle:
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


def write_selected(records, selected, out_path):
    with open(out_path, "w") as out:
        for out_name, source_name in selected:
            seq = records[source_name]
            out.write(f">{out_name}\n")
            for line in wrap(seq):
                out.write(line + "\n")


def build_fastas(args):
    records = {name: seq for name, header, seq in read_fasta(args.final_fasta)}
    pairs = []
    with open(args.homolog_pairs) as handle:
        header = next(handle).rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            pairs.append((fields[idx["pair_id"]], fields[idx["copy1_chr"]], fields[idx["copy2_chr"]]))

    if len(pairs) != args.chr_count:
        raise SystemExit(f"ERROR: expected {args.chr_count} homolog pairs, found {len(pairs)}")

    copy1_selected = []
    copy2_selected = []
    with open(args.copy_map, "w") as out:
        out.write("pair_id\tnormalized_chr\tcopy1_source\tcopy2_source\n")
        for i, (pair_id, copy1, copy2) in enumerate(pairs, start=1):
            out_name = f"chr{i:02d}"
            if copy1 not in records or copy2 not in records:
                raise SystemExit(f"ERROR: pair source not found in final FASTA: {pair_id} {copy1} {copy2}")
            copy1_selected.append((out_name, copy1))
            copy2_selected.append((out_name, copy2))
            out.write(f"{pair_id}\t{out_name}\t{copy1}\t{copy2}\n")

    write_selected(records, copy1_selected, args.copy1_out)
    write_selected(records, copy2_selected, args.copy2_out)

    ref_records = list(read_fasta(args.ref_fasta))
    selected_ref = []
    seen = set()
    for name, header, seq in ref_records:
        chr_num = None
        if args.ref_chr_mode == "ppa_super":
            match = re.search(r"contig: SUPER_(\d+)(?:[, \t]|$)", header)
            if match:
                chr_num = int(match.group(1))
        elif args.ref_chr_mode == "first_number":
            match = re.search(r"(\d+)", header)
            if match:
                chr_num = int(match.group(1))
        elif args.ref_chr_mode == "first_16_records":
            chr_num = len(selected_ref) + 1
        else:
            raise SystemExit(f"ERROR: unsupported REF_CHR_MODE: {args.ref_chr_mode}")

        if chr_num is None or chr_num < 1 or chr_num > args.chr_count or chr_num in seen:
            continue
        seen.add(chr_num)
        selected_ref.append((chr_num, seq))

    if len(selected_ref) != args.chr_count:
        raise SystemExit(
            f"ERROR: found {len(selected_ref)}/{args.chr_count} reference chromosomes with mode {args.ref_chr_mode}"
        )

    with open(args.ref_out, "w") as out:
        for chr_num, seq in sorted(selected_ref):
            out.write(f">chr{chr_num:02d}\n")
            for line in wrap(seq):
                out.write(line + "\n")


def summarize_copy_alignment(args):
    summary = {f"chr{i:02d}": {"+": 0, "-": 0, "other_chr": 0} for i in range(1, args.chr_count + 1)}

    with open(args.paf) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            qname, strand, tname = fields[0], fields[4], fields[5]
            aln_len = int(fields[10])
            if qname in summary and tname == qname:
                summary[qname][strand] += aln_len
            elif qname in summary:
                summary[qname]["other_chr"] += aln_len

    with open(args.out_tsv, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(
            ["chr", "same_chr_plus_aligned_bp", "same_chr_minus_aligned_bp", "other_chr_aligned_bp", "dominant_orientation"]
        )
        for i in range(1, args.chr_count + 1):
            chrom = f"chr{i:02d}"
            plus = summary[chrom]["+"]
            minus = summary[chrom]["-"]
            dominant = "+" if plus >= minus else "-"
            writer.writerow([chrom, plus, minus, summary[chrom]["other_chr"], dominant])


def orient_query(args):
    ref_names = [f"chr{i:02d}" for i in range(1, args.chr_count + 1)]
    query_names = [name for name, header, seq in read_fasta(args.query_in)]
    score = {}

    with open(args.paf) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            qname, strand, tname, aln_len = fields[0], fields[4], fields[5], int(fields[10])
            if qname not in query_names or tname not in ref_names or strand not in {"+", "-"}:
                continue
            key = (qname, tname)
            if key not in score:
                score[key] = {"+": 0, "-": 0}
            score[key][strand] += aln_len

    ranked = []
    for (query, ref), strands in score.items():
        total = strands["+"] + strands["-"]
        ranked.append((total, query, ref, strands["+"] , strands["-"]))
    ranked.sort(reverse=True)

    used_queries = set()
    used_refs = set()
    assignments = {}
    for total, query, ref, plus, minus in ranked:
        if query in used_queries or ref in used_refs:
            continue
        assignments[query] = {
            "ref": ref,
            "orientation": "+" if plus >= minus else "-",
            "plus": plus,
            "minus": minus,
            "total": total,
        }
        used_queries.add(query)
        used_refs.add(ref)

    missing_queries = [name for name in query_names if name not in assignments]
    missing_refs = [name for name in ref_names if name not in used_refs]
    if missing_queries or missing_refs:
        raise SystemExit(
            "ERROR: could not build one-to-one query/reference chromosome assignment. "
            f"Missing queries: {','.join(missing_queries)}; missing refs: {','.join(missing_refs)}"
        )

    with open(args.orientation_tsv, "w") as out:
        out.write("query_chr\tassigned_ref_chr\tnew_query_chr\torientation\tplus_aligned_bp\tminus_aligned_bp\ttotal_aligned_bp\n")
        for query in query_names:
            item = assignments[query]
            out.write(
                f"{query}\t{item['ref']}\t{item['ref']}\t{item['orientation']}\t"
                f"{item['plus']}\t{item['minus']}\t{item['total']}\n"
            )

    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    with open(args.query_out, "w") as out:
        for name, header, seq in read_fasta(args.query_in):
            item = assignments[name]
            if item["orientation"] == "-":
                seq = seq.translate(comp)[::-1]
            out.write(f">{item['ref']}\n")
            for line in wrap(seq):
                out.write(line + "\n")


def parse_args():
    parser = argparse.ArgumentParser(description="Helper commands for final homolog/SyRI workflow.")
    subparsers = parser.add_subparsers(dest="command")

    build = subparsers.add_parser("build-fastas")
    build.add_argument("--final-fasta", required=True)
    build.add_argument("--homolog-pairs", required=True)
    build.add_argument("--ref-fasta", required=True)
    build.add_argument("--ref-out", required=True)
    build.add_argument("--copy1-out", required=True)
    build.add_argument("--copy2-out", required=True)
    build.add_argument("--copy-map", required=True)
    build.add_argument("--ref-chr-mode", required=True)
    build.add_argument("--chr-count", type=int, required=True)
    build.set_defaults(func=build_fastas)

    summarize = subparsers.add_parser("summarize-copy-alignment")
    summarize.add_argument("--paf", required=True)
    summarize.add_argument("--out-tsv", required=True)
    summarize.add_argument("--chr-count", type=int, required=True)
    summarize.set_defaults(func=summarize_copy_alignment)

    orient = subparsers.add_parser("orient-query")
    orient.add_argument("--paf", required=True)
    orient.add_argument("--orientation-tsv", required=True)
    orient.add_argument("--query-in", required=True)
    orient.add_argument("--query-out", required=True)
    orient.add_argument("--chr-count", type=int, required=True)
    orient.set_defaults(func=orient_query)

    return parser.parse_args()


def main():
    args = parse_args()
    if not hasattr(args, "func"):
        raise SystemExit("ERROR: a subcommand is required.")
    args.func(args)


if __name__ == "__main__":
    main()
