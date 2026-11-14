import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/repeat_coverage.py"


class RepeatCoverageTests(unittest.TestCase):
    def test_tandem_priority_score_span_and_complete_denominator(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / "genome.fa"
            fasta.write_text(">chr1.full\n" + "A" * 160 + "N" * 40 + "\n")
            Path(str(fasta) + ".fai").write_text("chr1.full\t200\t11\t200\t201\n")
            gff = root / "repeats.gff3"
            gff.write_text(
                "chr\tEDTA\tLTR/Gypsy\t1\t150\t10\t+\t.\tName=TE1;classification=LTR/Gypsy\n"
                "chr\tEDTA\tLINE/L1\t51\t180\t20\t+\t.\tName=TE2;classification=LINE/L1\n"
                "chr\tEDTA\tLTR/Copia\t141\t190\t20\t+\t.\tName=TE3;classification=LTR/Copia\n"
            )
            lineage = root / "lineage.tsv"
            lineage.write_text(
                "record_id\tclass\torder\tfamily\tsublineage\tclade\n"
                "TE1\tClass_I\tLTR\tTy3/gypsy\tchromovirus\tTekay\n"
                "TE3\tClass_I\tLTR\tTy1/copia\tAle\t\n"
            )
            census = root / "census"
            (census / "pziA1").mkdir(parents=True)
            (census / "pziA1/monomer_intervals.tsv").write_text(
                "chromosome\tstart0\tend\tpublished_homology\nchr1.full\t20\t70\tTrue\nchr1.full\t70\t120\tFalse\n"
            )
            manifest = root / "targets.tsv"
            manifest.write_text(
                f"target\tfasta\tedta_gff\tlineage_table\nPZI-A1\t{fasta}\t{gff}\t{lineage}\n"
            )
            out = root / "out"
            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--manifest",
                    str(manifest),
                    "--census-dir",
                    str(census),
                    "--out",
                    str(out),
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            with (out / "exclusive_summary.tsv").open() as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(int(row["genome_bp"]), 200)
            self.assertEqual(int(row["exclusive_repeat_bp"]), 190)
            self.assertEqual(int(row["family_bp"]), 50)
            self.assertEqual(int(row["other_tandem_bp"]), 50)
            with (out / "exclusive_composition.tsv").open() as handle:
                counts = {
                    r["subfamily"]: int(r["bp"])
                    for r in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(counts["Tekay"], 20)
            self.assertEqual(counts["L1"], 60)
            self.assertEqual(counts["Ale"], 10)
            self.assertEqual(
                json.loads((out / "COVERAGE_DONE.json").read_text())["targets"], 1
            )


if __name__ == "__main__":
    unittest.main()
