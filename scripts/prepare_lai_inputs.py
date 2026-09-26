#!/usr/bin/env python3
import argparse
import gzip
import shutil
import subprocess
from pathlib import Path


def find_one(base_dir, patterns, label):
    candidates = []
    for pattern in patterns:
        candidates.extend(sorted(base_dir.glob(pattern)))
        candidates.extend(sorted(base_dir.rglob(pattern)))
    unique = []
    seen = set()
    for path in candidates:
        resolved = path.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        if path.is_file() and path.stat().st_size > 0:
            unique.append(path)
    if not unique:
        raise SystemExit("ERROR: could not find {} in {}".format(label, base_dir))
    return unique[0]


def nonempty(path):
    return path.is_file() and path.stat().st_size > 0


def find_genome_for_edta_gff(base_dir, gff_path):
    """Find the genome path matching an EDTA *.mod.EDTA.TEanno.gff3 file.

    The EDTA wiki runs LAI with the original genome filename used before EDTA
    adds the .mod.EDTA.TEanno.gff3 suffix. If that original FASTA is absent,
    fall back to the .mod FASTA because it has matching coordinates in many
    EDTA/LTR_retriever work directories.
    """
    name = gff_path.name
    suffixes = [".mod.EDTA.TEanno.gff3", ".EDTA.TEanno.gff3"]
    for suffix in suffixes:
        if name.endswith(suffix):
            stem = name[:-len(suffix)]
            original = base_dir / stem
            mod = base_dir / "{}.mod".format(stem)
            if nonempty(mod):
                return mod
            if nonempty(original):
                return original

    candidates = []
    for pattern in ["*.fa", "*.fasta", "*.fna"]:
        candidates.extend(sorted(base_dir.glob(pattern)))
    filtered = []
    for path in candidates:
        name = path.name
        if "TElib" in name or "masked" in name:
            continue
        if nonempty(path):
            filtered.append(path)
    if not filtered:
        raise SystemExit("ERROR: could not find a genome FASTA in {}".format(base_dir))
    return filtered[0]


def find_native_lai_inputs(base_dir):
    pass_patterns = ["*.real.pass.list"]
    all_patterns = ["*all_ltr.out", "*.all_ltr.out"]
    pass_candidates = []
    all_candidates = []
    for pattern in pass_patterns:
        pass_candidates.extend(sorted(base_dir.glob(pattern)))
        pass_candidates.extend(sorted(base_dir.rglob(pattern)))
    for pattern in all_patterns:
        all_candidates.extend(sorted(base_dir.glob(pattern)))
        all_candidates.extend(sorted(base_dir.rglob(pattern)))

    pass_candidates = [p for p in pass_candidates if nonempty(p)]
    all_candidates = [p for p in all_candidates if nonempty(p) and not p.name.endswith(".LAI")]
    if not pass_candidates or not all_candidates:
        return None
    pass_by_stem = {}
    for path in pass_candidates:
        name = path.name
        if name.endswith(".real.pass.list"):
            pass_by_stem[name[:-len(".real.pass.list")]] = path
    all_by_stem = {}
    for path in all_candidates:
        name = path.name
        if name.endswith(".all_ltr.out"):
            all_by_stem[name[:-len(".all_ltr.out")]] = path
        elif name.endswith("all_ltr.out"):
            all_by_stem[name[:-len("all_ltr.out")].rstrip("._-")] = path
    shared_stems = sorted(set(pass_by_stem).intersection(all_by_stem))
    if not shared_stems:
        return None

    genome_candidates = []
    for pattern in ["*.fa.mod", "*.fasta.mod", "*.fna.mod", "*.fa", "*.fasta", "*.fna"]:
        genome_candidates.extend(sorted(base_dir.glob(pattern)))
    genome_candidates = [
        p for p in genome_candidates
        if nonempty(p) and "TElib" not in p.name and "masked" not in p.name
    ]
    if not genome_candidates:
        return None

    stem = shared_stems[0]

    def rank_genome(path):
        name = path.name.lower()
        if name.endswith(".mod"):
            return 0
        return 1

    return {
        "genome": sorted(genome_candidates, key=lambda p: (rank_genome(p), p.name))[0],
        "intact": pass_by_stem[stem],
        "all": all_by_stem[stem],
    }


def run_gff2bed(gff_path, gff2bed_path, bed_path):
    with bed_path.open("w", encoding="utf-8") as out:
        subprocess.run(
            ["perl", str(gff2bed_path), str(gff_path), "structural"],
            check=True,
            stdout=out,
        )
    if not bed_path.is_file() or bed_path.stat().st_size == 0:
        raise SystemExit("ERROR: gff2bed produced an empty BED file: {}".format(bed_path))


def link_or_copy(src, dst):
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    try:
        dst.symlink_to(src.resolve())
    except OSError:
        shutil.copy2(src, dst)
    if not nonempty(dst):
        raise SystemExit("ERROR: failed to prepare LAI input {}".format(dst))


def normalize_genome_fasta(src, dst):
    """Write a LAI-safe FASTA with one whitespace-free ID per record.

    LAI/LTR_retriever versions differ in whether they match a FASTA record
    against the complete header or only its first token.  Haplotype FASTAs
    carry provenance text after the ID (for example ``C1A source_seqid=...``),
    so staging a one-token-header copy avoids silent empty LTR extraction.
    """
    opener = gzip.open if str(src).lower().endswith(".gz") else open
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_name(dst.name + ".tmp")
    seen = set()
    records = 0
    bases = 0
    current_id = None
    sequence = []

    def flush(out):
        nonlocal bases
        if current_id is None:
            return
        seq = "".join(sequence).replace(" ", "").replace("\t", "")
        if not seq:
            raise SystemExit("ERROR: empty sequence for {} in {}".format(current_id, src))
        bases += len(seq)
        for i in range(0, len(seq), 80):
            out.write(seq[i:i + 80] + "\n")

    try:
        with opener(src, "rt", encoding="utf-8", errors="replace") as inp, \
                tmp.open("w", encoding="utf-8") as out:
            for raw in inp:
                raw = raw.strip()
                if not raw:
                    continue
                if raw.startswith(">"):
                    flush(out)
                    token = raw[1:].split()
                    if not token:
                        raise SystemExit("ERROR: empty FASTA header in {}".format(src))
                    seq_id = token[0]
                    if seq_id in seen:
                        raise SystemExit(
                            "ERROR: duplicate FASTA sequence ID after normalization: {}".format(seq_id)
                        )
                    seen.add(seq_id)
                    records += 1
                    out.write(">{}\n".format(seq_id))
                    current_id = seq_id
                    sequence = []
                else:
                    if current_id is None:
                        raise SystemExit("ERROR: sequence data before FASTA header in {}".format(src))
                    sequence.append(raw)
            flush(out)
        if records == 0 or bases == 0:
            raise SystemExit("ERROR: genome FASTA is empty: {}".format(src))
        tmp.replace(dst)
    finally:
        if tmp.exists():
            tmp.unlink()


def bed_to_pass_and_rmout(bed_path, pass_path, rmout_path):
    pass_count = 0
    rmout_count = 0
    with bed_path.open(encoding="utf-8", errors="replace") as handle, \
            pass_path.open("w", encoding="utf-8") as pass_out, \
            rmout_path.open("w", encoding="utf-8") as rmout_out:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split()
            if len(fields) < 13:
                continue
            chrom = fields[0]
            start = fields[1]
            end = fields[2]
            anno = fields[3]
            intact_id = fields[6]
            direction = fields[8]
            superfamily = fields[12]
            joined = " ".join(fields).lower()
            if "ltr" in joined and "struc" in joined:
                pass_out.write("{}:{}..{}\t{}\n".format(chrom, start, end, intact_id))
                pass_count += 1
            rmout_out.write(
                "10000 0.001 0.001 0.001 {} {} {} NA {} {} {}\n".format(
                    chrom, start, end, direction, anno, superfamily
                )
            )
            rmout_count += 1
    if pass_count == 0:
        raise SystemExit("ERROR: no intact structural LTR records were written to {}".format(pass_path))
    if rmout_count == 0:
        raise SystemExit("ERROR: no LTR records were written to {}".format(rmout_path))


def parse_args():
    parser = argparse.ArgumentParser(description="Prepare LAI inputs from EDTA outputs.")
    parser.add_argument("--edta-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--gff2bed", required=True)
    parser.add_argument(
        "--method",
        choices=["auto", "native", "edta-gff3"],
        default="auto",
        help="Use native pass.list/all_ltr.out when available, or rebuild LAI inputs from EDTA GFF3.",
    )
    return parser.parse_args()


def resolve_executable_or_file(value, label):
    path = Path(value)
    if path.is_file():
        return path
    found = shutil.which(value)
    if found:
        return Path(found)
    raise SystemExit("ERROR: {} not found: {}".format(label, value))


def main():
    args = parse_args()
    edta_dir = Path(args.edta_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    gff = find_one(edta_dir, ["*.mod.EDTA.TEanno.gff3", "*.EDTA.TEanno.gff3"], "EDTA TE annotation GFF3")
    method = "edta-gff3"
    native = None
    if args.method in {"auto", "native"}:
        native = find_native_lai_inputs(edta_dir)
        if native:
            method = "native"
        elif args.method == "native":
            raise SystemExit("ERROR: native LAI inputs were not found in {}".format(edta_dir))

    bed = out_dir / "{}.EDTA.TEanno.struc.bed".format(args.prefix)
    if method == "native":
        genome_source = native["genome"]
        genome = out_dir / "{}.genome.fa".format(args.prefix)
        pass_list = out_dir / "{}.real.pass.list".format(args.prefix)
        rmout = out_dir / "{}.all_ltr.out".format(args.prefix)
        normalize_genome_fasta(genome_source, genome)
        link_or_copy(native["intact"], pass_list)
        link_or_copy(native["all"], rmout)
    else:
        genome_source = find_genome_for_edta_gff(edta_dir, gff)
        genome = out_dir / "{}.genome.fa".format(args.prefix)
        normalize_genome_fasta(genome_source, genome)
        pass_list = edta_dir / "{}.real.pass.list".format(args.prefix)
        rmout = edta_dir / "{}.all_ltr.out".format(args.prefix)
        gff2bed = resolve_executable_or_file(args.gff2bed, "gff2bed.pl")
        run_gff2bed(gff, gff2bed, bed)
        bed_to_pass_and_rmout(bed, pass_list, rmout)
        linked_pass_list = out_dir / pass_list.name
        linked_rmout = out_dir / rmout.name
        link_or_copy(pass_list, linked_pass_list)
        link_or_copy(rmout, linked_rmout)
        pass_list = linked_pass_list
        rmout = linked_rmout

    manifest = out_dir / "{}.lai_inputs.tsv".format(args.prefix)
    with manifest.open("w", encoding="utf-8") as out:
        out.write("key\tpath\n")
        out.write("genome\t{}\n".format(genome.resolve()))
        out.write("intact\t{}\n".format(pass_list.resolve()))
        out.write("all\t{}\n".format(rmout.resolve()))
        out.write("te_gff3\t{}\n".format(gff.resolve()))
        out.write("structural_bed\t{}\n".format("" if method == "native" else bed.resolve()))
        out.write("method\t{}\n".format(method))

    env_file = out_dir / "{}.lai_inputs.env".format(args.prefix)
    lai_result = out_dir / "{}.LAI".format(rmout.name)
    with env_file.open("w", encoding="utf-8") as out:
        out.write("GENOME='{}'\n".format(genome.resolve()))
        out.write("INTACT='{}'\n".format(pass_list.resolve()))
        out.write("ALL='{}'\n".format(rmout.resolve()))
        out.write("LAI_RESULT='{}'\n".format(lai_result.resolve()))
        out.write("TE_GFF3='{}'\n".format(gff.resolve()))
        out.write("STRUCTURAL_BED='{}'\n".format("" if method == "native" else bed.resolve()))
        out.write("LAI_INPUT_METHOD='{}'\n".format(method))


if __name__ == "__main__":
    main()
