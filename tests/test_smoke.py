from __future__ import absolute_import

import os
import shutil
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from prelnc import __version__
from prelnc.cli import hexamer_main, predict_main
from prelnc.config import FEATURE_SETS


def _data(name):
    return os.path.join(ROOT, "tests", "data", name)


class HexamerTableTest(unittest.TestCase):
    def test_matches_reference_implementation(self):
        """The golden file was produced by the original make_hexamer_tab.py."""
        expected = {}
        with open(_data("tiny_hexamer_expected.tsv")) as handle:
            next(handle)
            for line in handle:
                fields = line.split()
                expected[fields[0]] = (float(fields[1]), float(fields[2]))

        workdir = tempfile.mkdtemp()
        try:
            output = os.path.join(workdir, "hexamer.tsv")
            hexamer_main([
                "-c", _data("tiny_coding.fa"),
                "-n", _data("tiny_noncoding.fa"),
                "-o", output,
            ])
            with open(output) as handle:
                header = next(handle).split()
                rows = [line.split() for line in handle]
        finally:
            shutil.rmtree(workdir)

        self.assertEqual(header, ["hexamer", "coding", "noncoding"])
        self.assertEqual(len(rows), len(expected))
        self.assertEqual(set(row[0] for row in rows), set(expected))

        worst = 0.0
        for kmer, coding, noncoding in rows:
            worst = max(
                worst,
                abs(float(coding) - expected[kmer][0]),
                abs(float(noncoding) - expected[kmer][1]),
            )
        self.assertLess(worst, 1e-9)


class PredictionSmokeTest(unittest.TestCase):
    def test_predict_example(self):
        workdir = tempfile.mkdtemp()
        try:
            output = os.path.join(workdir, "example.tsv")
            predict_main([
                "-i", os.path.join(ROOT, "example", "example.fa"),
                "--species", "human",
                "-o", output,
                "-q",
            ])
            with open(output) as handle:
                lines = handle.read().strip().split("\n")
            with open(os.path.join(workdir, "example.features.tsv")) as handle:
                header = next(handle).rstrip("\n").split("\t")
                values = next(handle).rstrip("\n").split("\t")
            self.assertEqual(sorted(os.listdir(workdir)), ["example.features.tsv", "example.tsv"])
        finally:
            shutil.rmtree(workdir)

        self.assertEqual(lines[0].split("\t"), ["transcript_id", "prediction_score", "prediction"])
        self.assertEqual(lines[1].split("\t")[0], "ENST00000658521.1")
        self.assertGreaterEqual(float(lines[1].split("\t")[1]), 0.0)
        self.assertLessEqual(float(lines[1].split("\t")[1]), 1.0)
        self.assertIn(lines[1].split("\t")[2], ("lncRNA", "coding"))

        self.assertEqual(header, ["transcript_id"] + list(FEATURE_SETS["animal"]))
        self.assertEqual(len(values), len(FEATURE_SETS["animal"]) + 1)

    def test_existing_output_is_not_overwritten(self):
        workdir = tempfile.mkdtemp()
        try:
            output = os.path.join(workdir, "example.tsv")
            with open(output, "w") as handle:
                handle.write("keep me\n")
            args = ["-i", os.path.join(ROOT, "example", "example.fa"), "--species", "human", "-q"]
            predict_main(args + ["-o", output])
            predict_main(args + ["-o", output])

            with open(output) as handle:
                self.assertEqual(handle.read(), "keep me\n")
            self.assertTrue(os.path.exists(os.path.join(workdir, "example_1.tsv")))
            self.assertTrue(os.path.exists(os.path.join(workdir, "example_2.tsv")))
            self.assertTrue(os.path.exists(os.path.join(workdir, "example_2.features.tsv")))
        finally:
            shutil.rmtree(workdir)

    def test_default_output_goes_to_current_directory(self):
        workdir = tempfile.mkdtemp()
        previous = os.getcwd()
        try:
            os.chdir(workdir)
            predict_main([
                "-i", os.path.join(ROOT, "example", "example.fa"),
                "--species", "human",
                "-q",
            ])
            produced = sorted(os.listdir(workdir))
        finally:
            os.chdir(previous)
            shutil.rmtree(workdir)

        self.assertEqual(len(produced), 2)
        main_files = [name for name in produced if not name.endswith(".features.tsv")]
        feature_files = [name for name in produced if name.endswith(".features.tsv")]
        self.assertEqual(len(main_files), 1)
        self.assertEqual(len(feature_files), 1)
        self.assertTrue(main_files[0].startswith("prelnc_human_prediction_"))
        self.assertTrue(main_files[0].endswith(".tsv"))
        self.assertEqual(feature_files[0], main_files[0][: -len(".tsv")] + ".features.tsv")


class VersionTest(unittest.TestCase):
    def test_version_matches_package_metadata(self):
        self.assertEqual(__version__.count("."), 2)


if __name__ == "__main__":
    unittest.main()
