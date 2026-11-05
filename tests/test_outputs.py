from __future__ import absolute_import

import os
import shutil
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from prelnc.outputs import (
    OutputError,
    delimiter_for,
    resolve_outputs,
    write_table,
)


class DelimiterTest(unittest.TestCase):
    def test_tsv_family_uses_tabs(self):
        for extension in (".tsv", ".tab", ".txt", ".TSV"):
            self.assertEqual(delimiter_for("out" + extension), "\t")

    def test_csv_uses_comma(self):
        self.assertEqual(delimiter_for("out.csv"), ",")

    def test_unknown_extension_is_rejected(self):
        with self.assertRaises(OutputError):
            delimiter_for("out.xlsx")
        with self.assertRaises(OutputError):
            delimiter_for("out")


class ResolveOutputsTest(unittest.TestCase):
    def setUp(self):
        self.workdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.workdir)

    def _touch(self, name):
        path = os.path.join(self.workdir, name)
        with open(path, "w") as handle:
            handle.write("existing\n")
        return path

    def test_free_name_is_used_as_is(self):
        target = os.path.join(self.workdir, "out.tsv")
        plan = resolve_outputs(target, companions=(".features.tsv",))
        self.assertEqual(plan["main"], os.path.realpath(target) if os.path.islink(target) else target)
        self.assertFalse(plan["renamed"])
        self.assertEqual(plan["companions"], [os.path.join(self.workdir, "out.features.tsv")])

    def test_existing_file_gets_numbered_copy(self):
        self._touch("out.tsv")
        target = os.path.join(self.workdir, "out.tsv")
        plan = resolve_outputs(target, companions=(".features.tsv",))
        self.assertTrue(plan["renamed"])
        self.assertEqual(plan["main"], os.path.join(self.workdir, "out_1.tsv"))
        self.assertEqual(plan["companions"], [os.path.join(self.workdir, "out_1.features.tsv")])

    def test_companion_conflict_also_triggers_rename(self):
        # Only the companion exists: the pair must stay consistent.
        self._touch("out.features.tsv")
        plan = resolve_outputs(os.path.join(self.workdir, "out.tsv"), companions=(".features.tsv",))
        self.assertEqual(plan["main"], os.path.join(self.workdir, "out_1.tsv"))
        self.assertEqual(plan["companions"], [os.path.join(self.workdir, "out_1.features.tsv")])

    def test_force_keeps_requested_path(self):
        self._touch("out.tsv")
        plan = resolve_outputs(os.path.join(self.workdir, "out.tsv"), force=True)
        self.assertFalse(plan["renamed"])
        self.assertEqual(plan["main"], os.path.join(self.workdir, "out.tsv"))

    def test_default_name_is_timestamped_in_cwd(self):
        plan = resolve_outputs(default_stem="prelnc_test_prediction", cwd=self.workdir)
        name = os.path.basename(plan["main"])
        self.assertTrue(name.startswith("prelnc_test_prediction_"))
        self.assertTrue(name.endswith(".tsv"))
        self.assertEqual(os.path.dirname(plan["main"]), os.path.abspath(self.workdir))

    def test_missing_extension_is_rejected(self):
        with self.assertRaises(OutputError):
            resolve_outputs(os.path.join(self.workdir, "out"))

    def test_non_table_output_skips_extension_check(self):
        plan = resolve_outputs(
            os.path.join(self.workdir, "model.pkl"),
            table=False,
            default_extension=".pkl",
            companions=(".feature.tsv",),
        )
        self.assertEqual(plan["main"], os.path.join(self.workdir, "model.pkl"))
        self.assertEqual(plan["companions"], [os.path.join(self.workdir, "model.feature.tsv")])


class WriteTableTest(unittest.TestCase):
    def setUp(self):
        self.workdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.workdir)

    def test_tsv_content(self):
        path = os.path.join(self.workdir, "t.tsv")
        write_table(path, ["a", "b"], [[1, 2], [3, 4]])
        with open(path) as handle:
            self.assertEqual(handle.read(), "a\tb\n1\t2\n3\t4\n")

    def test_csv_content(self):
        path = os.path.join(self.workdir, "t.csv")
        write_table(path, ["a", "b"], [[1, 2]])
        with open(path) as handle:
            self.assertEqual(handle.read(), "a,b\n1,2\n")

    def test_value_containing_delimiter_is_quoted(self):
        path = os.path.join(self.workdir, "t.csv")
        write_table(path, ["a"], [["x,y"]])
        with open(path) as handle:
            self.assertEqual(handle.read(), 'a\n"x,y"\n')

    def test_no_temporary_file_is_left_behind(self):
        path = os.path.join(self.workdir, "t.tsv")
        write_table(path, ["a"], [[1]])
        self.assertEqual(sorted(os.listdir(self.workdir)), ["t.tsv"])


if __name__ == "__main__":
    unittest.main()
