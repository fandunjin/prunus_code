import unittest, tempfile, csv, subprocess, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from primary_locus_proteomes import locus
from locus_family_analysis import bh, contrast


class LocusTests(unittest.TestCase):
    def test_parent_identity_is_explicit(self):
        self.assertEqual(locus("LOC_1-mRNA-2 gene=LOC_1 cds_len=600"), "LOC_1")
        self.assertEqual(locus("gene-LOC2|rna-XM_2 description"), "gene-LOC2")
        self.assertEqual(locus("Gtr01g123.t2"), "Gtr01g123")
        with self.assertRaises(ValueError):
            locus("synthetic__g0000001")

    def test_isoforms_collapse_identical_paralogs_survive(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            pep = d / "frozen.fa"
            mp = d / "map.tsv"
            manifest = d / "manifest.tsv"
            out = d / "out"
            seqs = ["M" + "A" * 99, "M" + "A" * 198, "M" + "A" * 198, "M" + "A" * 97]
            specs = [
                ("s1", "LOC_1-mRNA-1", "LOC_1", 303),
                ("s2", "LOC_1-mRNA-2", "LOC_1", 600),
                ("s3", "LOC_2-mRNA-1", "LOC_2", 600),
                ("s4", "LOC_3-mRNA-1", "LOC_3", 297),
            ]
            pep.write_text(
                "".join(">" + r[0] + "\n" + s + "\n" for r, s in zip(specs, seqs))
            )
            with mp.open("w", newline="") as f:
                wr = csv.writer(f, delimiter="\t")
                wr.writerow(
                    ["stable_id", "original_id", "original_header", "protein_length"]
                )
                for (sid, oid, g, n), s in zip(specs, seqs):
                    wr.writerow([sid, oid, f"{oid} gene={g} cds_len={n}", len(s)])
            with manifest.open("w", newline="") as f:
                wr = csv.writer(f, delimiter="\t")
                wr.writerow(["track", "frozen_pep", "id_map", "restored_pep", "cds"])
                wr.writerow(["test", pep, mp, "", ""])
            command = [
                sys.executable,
                str(
                    Path(__file__).resolve().parents[1]
                    / "scripts/primary_locus_proteomes.py"
                ),
                "--manifest",
                str(manifest),
                "--out",
                str(out),
            ]
            p = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stderr)
            text = (out / "pep/test.fa").read_text()
            self.assertNotIn(">s1\n", text)
            self.assertIn(">s2\n", text)
            self.assertIn(">s3\n", text)
            self.assertNotIn(">s4\n", text)
            # A completed freeze cannot be silently overwritten by another run.
            self.assertNotEqual(
                subprocess.run(command, capture_output=True).returncode, 0
            )

    def test_full_precision_fdr_and_boundary_null(self):
        expected = [0.03, 0.04, 0.04]
        for got, want in zip(bh([0.01, 0.04, 0.03]), expected):
            self.assertAlmostEqual(got, want)
        exp, en, p, pl = contrast(0, 100, 0, 100)
        self.assertEqual((exp, p, pl), (0, 1, 1))
        _, en, p, _ = contrast(20, 100, 10, 100)
        self.assertEqual(en, 2)
        self.assertLess(p, 0.01)


if __name__ == "__main__":
    unittest.main()
