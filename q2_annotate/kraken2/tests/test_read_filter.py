# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import tempfile
import unittest
import os
from pathlib import Path

import pandas as pd

from qiime2.plugin.testing import TestPluginBase

from q2_annotate.kraken2.read_filter import (
    _assert_distinct_input_output_paths,
    _collect_matching_taxon_ids,
    _extract_matching_read_ids,
    _extract_matching_read_ids_from_output,
    _filter_paired_end_fastq,
    _filter_single_end_fastq,
    _manifest_to_sample_fastqs,
    _normalize_read_id,
    _resolve_fastq_path,
    _validate_read_sample_ids,
)


class TestReadFilterHelpers(TestPluginBase):
    package = "q2_annotate.kraken2.tests"

    @staticmethod
    def _load_report_df(fp: str) -> pd.DataFrame:
        df = pd.read_csv(fp, sep="\t", header=None)
        df.columns = [
            "perc_frags_covered",
            "n_frags_covered",
            "n_frags_assigned",
            "rank",
            "taxon_id",
            "name",
        ]
        return df

    @staticmethod
    def _load_output_df(fp: str) -> pd.DataFrame:
        df = pd.read_csv(fp, sep="\t", header=None)
        df.columns = [
            "classified",
            "sequence_id",
            "taxon_id",
            "read_length",
            "lca_or_minimizers",
        ]
        return df

    def test_normalize_read_id(self):
        self.assertEqual(_normalize_read_id("@abc/1"), "abc")
        self.assertEqual(_normalize_read_id("@abc/2"), "abc")
        self.assertEqual(_normalize_read_id("@abc 1:N:0:1"), "abc")
        self.assertEqual(_normalize_read_id("abc"), "abc")

    def test_collect_matching_taxon_ids_exact(self):
        report = self._load_report_df(
            self.get_data_path(
                "filter/output-report-alignment/reports/sample-1.report.txt"
            )
        )

        observed = _collect_matching_taxon_ids(
            report_df=report,
            taxonomy="Bacteria",
            include_descendants=True,
            contains=False,
        )

        self.assertIn("2", observed)
        self.assertIn("1767", observed)
        self.assertIn("1138383", observed)

    def test_collect_matching_taxon_ids_no_descendants(self):
        report = self._load_report_df(
            self.get_data_path(
                "filter/output-report-alignment/reports/sample-1.report.txt"
            )
        )

        observed = _collect_matching_taxon_ids(
            report_df=report,
            taxonomy="Bacteria",
            include_descendants=False,
            contains=False,
        )
        self.assertEqual(observed, {"2"})

    def test_extract_matching_read_ids(self):
        output = self._load_output_df(
            self.get_data_path(
                "filter/output-report-alignment/outputs/sample-1.output.txt"
            )
        )

        observed = _extract_matching_read_ids(output_df=output, taxon_ids={"1767"})

        self.assertIn("record-11", observed)
        self.assertIn("record-14", observed)
        self.assertNotIn("record-1", observed)

    def test_extract_matching_read_ids_from_output(self):
        output_fp = Path(
            self.get_data_path(
                "filter/output-report-alignment/outputs/sample-1.output.txt"
            )
        )
        observed = _extract_matching_read_ids_from_output(
            output_fp=output_fp, taxon_ids={"1767"}
        )

        self.assertIn("record-11", observed)
        self.assertIn("record-14", observed)
        self.assertNotIn("record-1", observed)

    def test_filter_single_end_fastq(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            input_fp = Path(tmpdir) / "input.fastq"
            output_fp = Path(tmpdir) / "output.fastq"

            input_fp.write_text(
                "@read-1\nAAAA\n+\n!!!!\n"
                "@read-2\nCCCC\n+\n####\n"
                "@read-3\nGGGG\n+\n$$$$\n"
            )

            _filter_single_end_fastq(
                input_fp=input_fp,
                output_fp=output_fp,
                matched_read_ids={"read-2"},
                exclude=False,
            )

            observed = output_fp.read_text()
            expected = "@read-2\nCCCC\n+\n####\n"
            self.assertEqual(observed, expected)

    def test_filter_paired_end_fastq(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            forward_input_fp = Path(tmpdir) / "forward.fastq"
            reverse_input_fp = Path(tmpdir) / "reverse.fastq"
            forward_output_fp = Path(tmpdir) / "forward.filtered.fastq"
            reverse_output_fp = Path(tmpdir) / "reverse.filtered.fastq"

            forward_input_fp.write_text(
                "@read-1/1\nAAAA\n+\n!!!!\n"
                "@read-2/1\nCCCC\n+\n####\n"
                "@read-3/1\nGGGG\n+\n$$$$\n"
            )
            reverse_input_fp.write_text(
                "@read-1/2\nTTTT\n+\n!!!!\n"
                "@read-2/2\nGGGG\n+\n####\n"
                "@read-3/2\nCCCC\n+\n$$$$\n"
            )

            _filter_paired_end_fastq(
                forward_input_fp=forward_input_fp,
                reverse_input_fp=reverse_input_fp,
                forward_output_fp=forward_output_fp,
                reverse_output_fp=reverse_output_fp,
                matched_read_ids={"read-2"},
                exclude=False,
            )

            self.assertEqual(
                forward_output_fp.read_text(), "@read-2/1\nCCCC\n+\n####\n"
            )
            self.assertEqual(
                reverse_output_fp.read_text(), "@read-2/2\nGGGG\n+\n####\n"
            )

    def test_filter_paired_end_fastq_mismatch(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            forward_input_fp = Path(tmpdir) / "forward.fastq"
            reverse_input_fp = Path(tmpdir) / "reverse.fastq"
            forward_output_fp = Path(tmpdir) / "forward.filtered.fastq"
            reverse_output_fp = Path(tmpdir) / "reverse.filtered.fastq"

            forward_input_fp.write_text("@read-1/1\nAAAA\n+\n!!!!\n")
            reverse_input_fp.write_text("@read-X/2\nTTTT\n+\n!!!!\n")

            with self.assertRaisesRegex(ValueError, "out of sync"):
                _filter_paired_end_fastq(
                    forward_input_fp=forward_input_fp,
                    reverse_input_fp=reverse_input_fp,
                    forward_output_fp=forward_output_fp,
                    reverse_output_fp=reverse_output_fp,
                    matched_read_ids={"read-1"},
                    exclude=False,
                )

    def test_assert_distinct_input_output_paths(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "reads.fastq"
            filepath.write_text("")

            with self.assertRaisesRegex(ValueError, "identical"):
                _assert_distinct_input_output_paths(filepath, filepath)

    def test_assert_distinct_input_output_paths_hardlink(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            input_fp = Path(tmpdir) / "reads.fastq"
            output_fp = Path(tmpdir) / "reads_copy.fastq"
            input_fp.write_text("@r1\nAAAA\n+\n!!!!\n")
            os.link(input_fp, output_fp)

            with self.assertRaisesRegex(ValueError, "same inode"):
                _assert_distinct_input_output_paths(input_fp, output_fp)

    def test_resolve_fastq_path_relative(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            reads_dir = Path(tmpdir)
            observed = _resolve_fastq_path("reads.fastq.gz", reads_dir)
            self.assertEqual(observed, reads_dir / "reads.fastq.gz")

    def test_manifest_to_sample_fastqs_long_format(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            reads_dir = Path(tmpdir)
            manifest = pd.DataFrame(
                [
                    {
                        "sample-id": "sample1",
                        "filename": "sample1_R1.fastq.gz",
                        "direction": "forward",
                    },
                    {
                        "sample-id": "sample1",
                        "filename": "sample1_R2.fastq.gz",
                        "direction": "reverse",
                    },
                    {
                        "sample-id": "sample2",
                        "filename": "sample2_R1.fastq.gz",
                        "direction": "forward",
                    },
                    {
                        "sample-id": "sample2",
                        "filename": "sample2_R2.fastq.gz",
                        "direction": "reverse",
                    },
                ]
            )

            observed = _manifest_to_sample_fastqs(manifest, reads_dir)
            self.assertEqual(
                observed["sample1"]["forward"], reads_dir / "sample1_R1.fastq.gz"
            )
            self.assertEqual(
                observed["sample1"]["reverse"], reads_dir / "sample1_R2.fastq.gz"
            )
            self.assertEqual(
                observed["sample2"]["forward"], reads_dir / "sample2_R1.fastq.gz"
            )
            self.assertEqual(
                observed["sample2"]["reverse"], reads_dir / "sample2_R2.fastq.gz"
            )

    def test_manifest_to_sample_fastqs_mixed_single_and_paired_raises(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            reads_dir = Path(tmpdir)
            manifest = pd.DataFrame(
                [
                    {
                        "sample-id": "sample1",
                        "filename": "sample1_R1.fastq.gz",
                        "direction": "forward",
                    },
                    {
                        "sample-id": "sample1",
                        "filename": "sample1_R2.fastq.gz",
                        "direction": "reverse",
                    },
                    {
                        "sample-id": "sample2",
                        "filename": "sample2_R1.fastq.gz",
                        "direction": "forward",
                    },
                ]
            )

            with self.assertRaisesRegex(ValueError, "mix of single-end and paired-end"):
                _manifest_to_sample_fastqs(manifest, reads_dir)

    def test_validate_read_sample_ids_allows_extra_read_samples(self):
        _validate_read_sample_ids(
            read_sample_ids={"sample1", "sample2"},
            report_sample_ids={"sample1"},
            output_sample_ids={"sample1"},
        )

    def test_validate_read_sample_ids_missing_classification_sample(self):
        with self.assertRaisesRegex(ValueError, "missing from reads"):
            _validate_read_sample_ids(
                read_sample_ids={"sample1"},
                report_sample_ids={"sample1", "sample2"},
                output_sample_ids={"sample1", "sample2"},
            )


if __name__ == "__main__":
    unittest.main()
