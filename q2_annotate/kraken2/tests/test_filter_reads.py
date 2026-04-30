# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import gzip
import os
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from qiime2.plugin.testing import TestPluginBase
from q2_types.kraken2 import (
    Kraken2OutputDirectoryFormat,
    Kraken2ReportDirectoryFormat,
    Kraken2ReportFormat,
)
from q2_types.per_sample_sequences import CasavaOneEightSingleLanePerSampleDirFmt

from q2_annotate.kraken2.filter_reads import (
    _assert_distinct_input_output_paths,
    _collect_matching_taxon_ids,
    _extract_matching_read_ids_from_output,
    _filter_paired_end_fastq,
    _filter_single_end_fastq,
    _normalize_read_id,
    _validate_read_sample_ids,
    filter_reads_kraken2,
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
        report = Path(
            self.get_data_path(
                "filter/output-report-alignment/reports/sample-1.report.txt"
            )
        )

        observed = _collect_matching_taxon_ids(
            report=report,
            taxonomy="Bacteria",
            include_descendants=True,
            contains=False,
        )

        self.assertSetEqual(
            {
                "1138383",
                "120793",
                "1760",
                "1762",
                "1763",
                "1767",
                "1783272",
                "2",
                "201174",
                "85007",
            },
            observed,
        )

    def test_collect_matching_taxon_ids_no_descendants(self):
        report = Path(
            self.get_data_path(
                "filter/output-report-alignment/reports/sample-1.report.txt"
            )
        )

        observed = _collect_matching_taxon_ids(
            report=report,
            taxonomy="Bacteria",
            include_descendants=False,
            contains=False,
        )
        self.assertSetEqual({"2"}, observed)

    def test_extract_matching_read_ids_from_output(self):
        output_fp = Path(
            self.get_data_path(
                "filter/output-report-alignment/outputs/sample-1.output.txt"
            )
        )
        observed = _extract_matching_read_ids_from_output(
            output_fp=output_fp, taxon_ids={"1767"}
        )

        self.assertSetEqual({"record-11", "record-14"}, observed)

    def test_filter_single_end_fastq(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            input_fp = Path(self.get_data_path("read-filter/fastq/single-end.fastq"))
            output_fp = Path(tmpdir) / "output.fastq"

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
            forward_input_fp = Path(
                self.get_data_path("read-filter/fastq/forward.fastq")
            )
            reverse_input_fp = Path(
                self.get_data_path("read-filter/fastq/reverse.fastq")
            )
            forward_output_fp = Path(tmpdir) / "forward.filtered.fastq"
            reverse_output_fp = Path(tmpdir) / "reverse.filtered.fastq"

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
            forward_input_fp = Path(
                self.get_data_path("read-filter/fastq/forward.fastq")
            )
            reverse_input_fp = Path(
                self.get_data_path("read-filter/fastq/reverse-mismatched.fastq")
            )
            forward_output_fp = Path(tmpdir) / "forward.filtered.fastq"
            reverse_output_fp = Path(tmpdir) / "reverse.filtered.fastq"

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
            filepath.write_text("abc")

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

    def test_validate_read_sample_ids_rejects_extra_read_samples(self):
        with self.assertRaisesRegex(ValueError, "missing from Kraken2"):
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


class TestFilterKraken2Reads(TestPluginBase):
    package = "q2_annotate.kraken2.tests"

    @staticmethod
    def _read_gzip_text(filepath: Path) -> str:
        with gzip.open(filepath, mode="rt") as fh:
            return fh.read()

    def setUp(self):
        super().setUp()

        self.reports = Kraken2ReportDirectoryFormat(
            self.get_data_path("read-filter/reports"), mode="r"
        )
        self.outputs = Kraken2OutputDirectoryFormat(
            self.get_data_path("read-filter/outputs"), mode="r"
        )

    def test_filter_paired_end_reads_by_taxonomy(self):
        reads = CasavaOneEightSingleLanePerSampleDirFmt(
            self.get_data_path("read-filter/paired-end-reads"), mode="r"
        )

        observed = filter_reads_kraken2(
            reads=reads,
            reports=self.reports,
            outputs=self.outputs,
            taxonomy="Bacteria",
        )

        self.assertIsInstance(observed, CasavaOneEightSingleLanePerSampleDirFmt)
        self.assertEqual(
            self._read_gzip_text(
                Path(observed.path) / "sample1_0_L001_R1_001.fastq.gz"
            ),
            "@record-1/1\nAAAA\n+\n!!!!\n@record-3/1\nGGGG\n+\n$$$$\n",
        )
        self.assertEqual(
            self._read_gzip_text(
                Path(observed.path) / "sample1_0_L001_R2_001.fastq.gz"
            ),
            "@record-1/2\nTTTT\n+\n!!!!\n@record-3/2\nCCCC\n+\n$$$$\n",
        )

    def test_filter_single_end_reads_with_inverse_filter(self):
        reads = CasavaOneEightSingleLanePerSampleDirFmt(
            self.get_data_path("read-filter/single-end-reads"), mode="r"
        )

        observed = filter_reads_kraken2(
            reads=reads,
            reports=self.reports,
            outputs=self.outputs,
            taxonomy="Bacteria",
            exclude=True,
        )

        self.assertEqual(
            self._read_gzip_text(
                Path(observed.path) / "sample1_0_L001_R1_001.fastq.gz"
            ),
            "@record-2\nCCCC\n+\n####\n",
        )

    def test_filter_reads_raises_for_missing_taxonomy(self):
        reads = CasavaOneEightSingleLanePerSampleDirFmt(
            self.get_data_path("read-filter/single-end-reads"), mode="r"
        )

        with self.assertRaisesRegex(ValueError, "was not found"):
            filter_reads_kraken2(
                reads=reads,
                reports=self.reports,
                outputs=self.outputs,
                taxonomy="Archaea",
            )


if __name__ == "__main__":
    unittest.main()
