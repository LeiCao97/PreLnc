from __future__ import absolute_import

import io
import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from prelnc.progress import Progress, _display_width, pad, print_summary, truncate


class ProgressTest(unittest.TestCase):
    def test_disabled_progress_writes_nothing(self):
        stream = io.StringIO()
        progress = Progress(100, "test", stream=stream, enabled=False)
        for _ in range(100):
            progress.advance()
        progress.finish()
        self.assertEqual(stream.getvalue(), "")

    def test_empty_total_is_silent(self):
        stream = io.StringIO()
        progress = Progress(0, "test", stream=stream)
        progress.advance()
        progress.finish()
        self.assertEqual(stream.getvalue(), "")

    def test_redirected_stream_prints_at_most_twenty_lines(self):
        stream = io.StringIO()
        progress = Progress(100, "CDS", stream=stream)
        for _ in range(100):
            progress.advance()
        progress.finish()
        lines = [line for line in stream.getvalue().split("\n") if line.strip()]
        self.assertLessEqual(len(lines), 21)
        self.assertTrue(lines[-1].endswith("100/100 items · elapsed 0s · eta 0s"), lines[-1])
        self.assertEqual(stream.getvalue().count("\r"), 0)

    def test_small_batch_still_reports_every_item(self):
        stream = io.StringIO()
        progress = Progress(3, "hexamer", stream=stream)
        for _ in range(3):
            progress.advance()
        progress.finish()
        lines = [line for line in stream.getvalue().split("\n") if line.strip()]
        self.assertEqual(len(lines), 3)

    def test_summary_block_is_aligned(self):
        stream = io.StringIO()
        print_summary(
            [("Input", "a.fa (2 sequences)"), ("Predictions", "out.tsv")],
            title="Prediction finished",
            stream=stream,
        )
        lines = stream.getvalue().split("\n")
        self.assertEqual(lines[0], "[PreLnc] Prediction finished")
        self.assertIn("a.fa (2 sequences)", lines[1])
        self.assertIn("out.tsv", lines[2])
        # Values must start at the same display column.
        self.assertEqual(
            _display_width(lines[1][: lines[1].index(":")]),
            _display_width(lines[2][: lines[2].index(":")]),
        )

    def test_summary_has_no_non_ascii_text(self):
        stream = io.StringIO()
        print_summary([("Input", "a.fa"), ("Elapsed", "1s")], title="Prediction finished", stream=stream)
        for line in stream.getvalue().split("\n"):
            self.assertTrue(all(ord(char) < 128 for char in line), line)


class WidthTest(unittest.TestCase):
    def test_pad_counts_wide_characters_as_two_columns(self):
        self.assertEqual(pad("结果表", 8), "结果表  ")
        self.assertEqual(pad("abc", 8), "abc     ")

    def test_truncate_respects_display_width(self):
        self.assertEqual(truncate("结果表", 6), "结果表")
        self.assertEqual(truncate("结果表", 5), "结果")
        self.assertEqual(truncate("abcdef", 4), "abcd")


if __name__ == "__main__":
    unittest.main()
