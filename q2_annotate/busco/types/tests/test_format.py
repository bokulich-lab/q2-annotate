# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from qiime2.core.exceptions import ValidationError
from qiime2.plugin.testing import TestPluginBase
from q2_annotate.busco.types import BUSCOResultsFormat, BuscoDatabaseDirFmt


class TestBUSCOFormats(TestPluginBase):
    package = "q2_annotate.busco.types.tests"

    def test_busco_results_format_ok(self):
        results = BUSCOResultsFormat(self.get_data_path("busco_results.tsv"), mode="r")
        results.validate(level="min")
        results.validate(level="max")

    def test_busco_results_format_ok_header_2(self):
        results = BUSCOResultsFormat(
            self.get_data_path("busco_results_header_2.tsv"), mode="r"
        )
        results.validate(level="min")
        results.validate(level="max")

    def test_busco_results_format_with_optional_columns(self):
        results = BUSCOResultsFormat(self.get_data_path("unbinned_res.tsv"), mode="r")
        results.validate()

    def test_busco_results_format_error_header(self):
        results = BUSCOResultsFormat(
            self.get_data_path("busco_results_broken_header.tsv"), mode="r"
        )
        with self.assertRaisesRegex(ValidationError, "Invalid header"):
            results.validate()

    def test_busco_results_format_error_values(self):
        results = BUSCOResultsFormat(
            self.get_data_path("busco_results_broken_values.tsv"), mode="r"
        )
        with self.assertRaisesRegex(ValidationError, "Line 4 has 14 columns"):
            results.validate()

    def test_BuscoDatabaseDirFmt_valid(self):
        dirpath = self.get_data_path("busco_db_valid")
        format = BuscoDatabaseDirFmt(dirpath, mode="r")
        format.validate()

    def test_BuscoDatabaseDirFmt_invalid(self):
        dirpath = self.get_data_path("busco_db_invalid")
        format = BuscoDatabaseDirFmt(dirpath, mode="r")
        with self.assertRaisesRegex(ValidationError, "Missing one or more files"):
            format.validate()
