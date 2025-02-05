# ----------------------------------------------------------------------------
# Copyright (c) 2023-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os.path

import pandas as pd
import qiime2
from q2_types.kraken2 import Kraken2ReportDirectoryFormat
from qiime2.plugin.testing import TestPluginBase

from q2_moshpit.kraken2.filter import _validate_parameters, \
    _find_empty_reports, _create_filtered_results, \
    filter_kraken_reports


class TestFilterKrakenReports(TestPluginBase):
    package = "q2_moshpit.kraken2.tests"

    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        instance = cls()

        cls.report_mags = Kraken2ReportDirectoryFormat(
            instance.get_data_path("reports-mags"), "r"
        )
        cls.report_mags_unclassified_missing_frac = (
            Kraken2ReportDirectoryFormat(
                instance.get_data_path(
                    "reports-mags-unclassified-missing-frac"),
                "r"))

        cls.file_dict_report_mags = cls.report_mags.file_dict()
        cls.file_dict_report_unclassified = (
            cls.report_mags_unclassified_missing_frac.file_dict())

        cls.metadata_df = pd.read_csv(
            instance.get_data_path("metadata/metadata.tsv"),
            sep="\t",
            index_col="id"
        )

        cls.metadata1 = qiime2.Metadata(cls.metadata_df)

        cls.metadata_df.drop(
            "8894435a-c836-4c18-b475-8b38a9ab6c6b", inplace=True)
        cls.metadata2 = qiime2.Metadata(cls.metadata_df)

    def test_find_empty_reports(self):
        empty_reports = _find_empty_reports(
            file_dict={"": self.file_dict_report_mags}
        )
        self.assertEqual(
            empty_reports, {"8894435a-c836-4c18-b475-8b38a9ab6c6b"}
        )

    def test_find_empty_reports_missing_frac(self):
        empty_reports = _find_empty_reports(
            file_dict={"": self.file_dict_report_unclassified}
        )
        self.assertEqual(
            empty_reports, {"8894435a-c836-4c18-b475-8b38a9ab6c6b"}
        )

    def test_create_filter_results_reports(self):
        results = _create_filtered_results(
            file_dict={"sample1": self.file_dict_report_unclassified},
            ids_to_keep={"8894435a-c836-4c18-b475-8b38a9ab6c6b"}
        )
        self.assertTrue(
            os.path.exists(
                os.path.join(
                    str(results),
                    "sample1",
                    "8894435a-c836-4c18-b475-8b38a9ab6c6b.report.txt"
                )
            )
        )

    def test_filter_kraken_reports_error(self):
        with self.assertRaisesRegex(
                ValueError, "No IDs remain after filtering."
        ):
            filter_kraken_reports(
                reports=self.report_mags_unclassified_missing_frac,
                metadata=self.metadata1,
                exclude_ids=True
            )

    def test_filter_kraken_reports_metadata(self):
        results = filter_kraken_reports(
            reports=self.report_mags_unclassified_missing_frac,
            metadata=self.metadata2,
        )
        self.assertTrue(
            os.path.exists(
                os.path.join(
                    str(results),
                    "3b72d1a7-ddb0-4dc7-ac36-080ceda04aaa.report.txt"
                )
            )
        )

    def test_missing_metadata_and_remove_empty(self):
        with self.assertRaisesRegex(
                ValueError, r'--m-metadata-file.*--p-remove-empty'
        ):
            _validate_parameters(metadata=None, remove_empty=False,
                                 where=None, exclude_ids=False)

    def test_where_without_metadata(self):
        with self.assertRaisesRegex(
                ValueError, r'--p-where.*--m-metadata-file'
        ):
            _validate_parameters(metadata=None, remove_empty=True,
                                 where=True, exclude_ids=False)

    def test_exclude_ids_without_metadata(self):
        with self.assertRaisesRegex(
                ValueError, r'--p-exclude-ids.*--m-metadata-file'
        ):
            _validate_parameters(metadata=None, remove_empty=True,
                                 where=None, exclude_ids=True)
