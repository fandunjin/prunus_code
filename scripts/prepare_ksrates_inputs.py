#!/usr/bin/env python3
"""Prepare primary-haplotype CDS and GFF inputs for three ksrates analyses."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


FOCAL_COMPONENTS = {
    "pzi": ["pzi_hapA1", "pzi_hapB1"],
    "psp": ["psp_hapA1", "psp_hapB1"],
    "ppa": ["ppa_hapA", "ppa_hapB"],
}

LATIN_NAMES = {
    "pzi": "Prunus zippeliana",
    "psp": "Prunus spinulosa",
    "ppa": "Prunus padus",
    "pcam": "Prunus campanulata",
    "ppe": "Prunus persica",
    "gtr": "Gillenia trifoliata",
}

REFERENCES = ["pcam", "ppe", "gtr"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jcvi-input-root", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--bootstrap-iterations", type=int, default=1000)
    parser.add_argument("--pipeline-version", default="v2.1.1")
    parser.add_argument("--threads", type=int, default=32)
    return parser.parse_args()


def require(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit("Missing input: {}".format(path))


def write_tsv(path: Path, rows: Iterable[dict], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def safe_attribute(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.:-]", "_", value)


def fasta_ids(path: Path) -> List[str]:
    ids = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


def append_prefixed_fasta(source: Path, destination, prefix: str, seen: set) -> int:
    count = 0
    with source.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                original = line[1:].strip().split()[0]
                output_id = safe_attribute("{}__{}".format(prefix, original))
                if output_id in seen:
                    raise SystemExit("Duplicate output FASTA ID: {}".format(output_id))
                seen.add(output_id)
                destination.write(">{}\n".format(output_id))
                count += 1
            else:
                destination.write(line)
    return count


def append_gff_from_bed(
    bed: Path,
    destination,
    component: str,
    cds_ids: set,
) -> Tuple[int, int]:
    written = 0
    missing = 0
    with bed.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise SystemExit("Malformed BED line {}:{}".format(bed, line_number))
            seqid, start0, end, original_id, _score, strand = fields[:6]
            if original_id not in cds_ids:
                missing += 1
                continue
            output_id = "{}__{}".format(component, original_id)
            output_seqid = "{}__{}".format(component, seqid)
            destination.write(
                "{}\tJCVI\tmRNA\t{}\t{}\t.\t{}\t.\tID={};Name={}\n".format(
                    safe_attribute(output_seqid),
                    int(start0) + 1,
                    int(end),
                    strand,
                    safe_attribute(output_id),
                    safe_attribute(output_id),
                )
            )
            written += 1
    return written, missing


def multiline_mapping(items: Sequence[Tuple[str, str]], indent: int = 20) -> str:
    padding = " " * indent
    return (",\n" + padding).join("{} : {}".format(key, value) for key, value in items)


def write_ksrates_config(
    path: Path,
    focal: str,
    sequence_paths: Dict[str, Path],
    gff_path: Path,
    bootstrap_iterations: int,
) -> None:
    taxa = [focal, "pcam", "ppe", "gtr"]
    fasta_mapping = multiline_mapping([(taxon, str(sequence_paths[taxon].resolve())) for taxon in taxa])
    latin_mapping = multiline_mapping([(taxon, LATIN_NAMES[taxon]) for taxon in taxa])
    text = """[SPECIES]
focal_species = {focal}
newick_tree = (gtr,({focal},(pcam,ppe)));
latin_names =       {latin_mapping}
fasta_filenames =   {fasta_mapping}
gff_filename = {gff}
peak_database_path = {db_dir}/{focal}_ortholog_peak_db.tsv
ks_list_database_path = {db_dir}/{focal}_ortholog_ks_list_db.tsv

[ANALYSIS SETTING]
paranome = yes
collinearity = yes
reciprocal_retention = no
gff_feature = mRNA
gff_attribute = id
max_number_outgroups = 3
consensus_mode_for_multiple_outgroups = mean among outgroups

[PARAMETERS]
x_axis_max_limit_paralogs_plot = 3
bin_width_paralogs = 0.05
y_axis_max_limit_paralogs_plot = None
num_bootstrap_iterations = {bootstraps}
divergence_colors = Goldenrod, MediumBlue, Gray
x_axis_max_limit_orthologs_plots = 5
bin_width_orthologs = 0.05
max_ks_paralogs = 3
max_ks_orthologs = 5
""".format(
        focal=focal,
        latin_mapping=latin_mapping,
        fasta_mapping=fasta_mapping,
        gff=str(gff_path.resolve()),
        db_dir=str((path.parent.parent / "databases").resolve()),
        bootstraps=bootstrap_iterations,
    )
    path.write_text(text, encoding="utf-8")


def write_expert_config(path: Path) -> None:
    path.write_text(
        """[EXPERT PARAMETERS]
logging_level = info
preserve_ks_tmp_files = no
kde_bandwidth_modifier = 0.4
plot_adjustment_arrows = yes
max_mixture_model_iterations = 600
num_mixture_model_initializations = 30
extra_paralogs_analyses_methods = yes
max_mixture_model_components = 3
max_mixture_model_ks = 3
max_gene_family_size = 200
num_reciprocally_retained_gfs = 2000
use_bottom_gfs_instead_of_top = no
use_original_orthomcl_version = no
min_ks_anchor_pairs = 0.05
""",
        encoding="utf-8",
    )


def write_nextflow_config(path: Path, threads: int, pipeline_version: str) -> None:
    text = """profiles {{
    apptainer {{
        apptainer.enabled = true
        apptainer.autoMounts = true
    }}
    singularity {{
        singularity.enabled = true
        singularity.autoMounts = true
    }}
}}

executor {{
    name = 'local'
    cpus = {threads}
    queueSize = 2
}}

process {{
    container = 'vibpsb/ksrates:{version}'
    errorStrategy = 'terminate'

    withName: 'wgdParalogs' {{
        cpus = {heavy_threads}
        memory = '80 GB'
    }}
    withName: 'wgdOrthologs' {{
        cpus = {heavy_threads}
        memory = '80 GB'
    }}
    withName: 'estimatePeaks' {{
        memory = '16 GB'
    }}
    withName: 'paralogsAnalyses' {{
        memory = '32 GB'
    }}
    withName: 'plotOrthologDistrib' {{
        memory = '16 GB'
    }}
}}
""".format(threads=threads, heavy_threads=max(1, threads // 2), version=pipeline_version)
    path.write_text(text, encoding="utf-8")


def main() -> int:
    args = parse_args()
    input_root = Path(args.jcvi_input_root)
    outdir = Path(args.outdir)
    sequences_dir = outdir / "sequences"
    annotations_dir = outdir / "annotations"
    configs_dir = outdir / "configs"
    databases_dir = outdir / "databases"
    for directory in (sequences_dir, annotations_dir, configs_dir, databases_dir):
        directory.mkdir(parents=True, exist_ok=True)

    all_labels = sorted({x for values in FOCAL_COMPONENTS.values() for x in values} | set(REFERENCES))
    sources = {}
    for label in all_labels:
        cds = input_root / label / "{}.cds".format(label)
        bed = input_root / label / "{}.bed".format(label)
        require(cds)
        require(bed)
        sources[label] = (cds, bed)

    sequence_paths: Dict[str, Path] = {}
    qc_rows = []
    for focal, components in FOCAL_COMPONENTS.items():
        fasta_path = sequences_dir / "{}.fasta".format(focal)
        gff_path = annotations_dir / "{}.gff3".format(focal)
        seen = set()
        sequence_count = 0
        gff_count = 0
        missing_count = 0
        with fasta_path.open("w", encoding="utf-8") as fasta_out, gff_path.open("w", encoding="utf-8") as gff_out:
            gff_out.write("##gff-version 3\n")
            for component in components:
                cds, bed = sources[component]
                ids = set(fasta_ids(cds))
                sequence_count += append_prefixed_fasta(cds, fasta_out, component, seen)
                written, missing = append_gff_from_bed(bed, gff_out, component, ids)
                gff_count += written
                missing_count += missing
        if sequence_count != gff_count or missing_count:
            raise SystemExit(
                "FASTA/GFF mismatch for {}: fasta={} gff={} missing={}".format(
                    focal, sequence_count, gff_count, missing_count
                )
            )
        sequence_paths[focal] = fasta_path
        qc_rows.append(
            {
                "species": focal,
                "display_name": LATIN_NAMES[focal],
                "components": ",".join(components),
                "cds_records": sequence_count,
                "gff_records": gff_count,
                "missing_bed_ids": missing_count,
                "status": "PASS",
            }
        )

    for reference in REFERENCES:
        cds, _bed = sources[reference]
        fasta_path = sequences_dir / "{}.fasta".format(reference)
        seen = set()
        with fasta_path.open("w", encoding="utf-8") as fasta_out:
            count = append_prefixed_fasta(cds, fasta_out, reference, seen)
        sequence_paths[reference] = fasta_path
        qc_rows.append(
            {
                "species": reference,
                "display_name": LATIN_NAMES[reference],
                "components": reference,
                "cds_records": count,
                "gff_records": "NA",
                "missing_bed_ids": "NA",
                "status": "PASS",
            }
        )

    for focal in FOCAL_COMPONENTS:
        write_ksrates_config(
            configs_dir / "config_{}.txt".format(focal),
            focal,
            sequence_paths,
            annotations_dir / "{}.gff3".format(focal),
            args.bootstrap_iterations,
        )
    write_expert_config(configs_dir / "config_expert.txt")
    write_nextflow_config(configs_dir / "nextflow.config", args.threads, args.pipeline_version)
    write_tsv(
        outdir / "input_qc.tsv",
        qc_rows,
        ["species", "display_name", "components", "cds_records", "gff_records", "missing_bed_ids", "status"],
    )
    (outdir / "species_tree.newick").write_text(
        "# Each focal analysis uses (gtr,(focal,(pcam,ppe)));\n",
        encoding="utf-8",
    )
    (outdir / "PREPARED").write_text("ksrates_inputs_prepared\n", encoding="utf-8")
    print("[OK] ksrates inputs prepared: {}".format(outdir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
