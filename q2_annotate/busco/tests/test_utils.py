# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import json
import os
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock

import pandas as pd
from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt
from qiime2.plugin.testing import TestPluginBase

from q2_annotate.busco.types import BuscoDatabaseDirFmt
from q2_annotate.busco.utils import (
    _parse_busco_params,
    _parse_df_columns,
    _partition_dataframe,
    _get_feature_table,
    _calculate_summary_stats,
    _validate_lineage_dataset_input,
    _extract_json_data,
    _calculate_contamination_completeness,
    _process_busco_results,
    _calculate_unbinned_percentage,
    _count_contigs,
    _add_unbinned_metrics,
)


class TestBUSCOUtils(TestPluginBase):
    package = "q2_annotate.busco.tests"

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.mags = MultiMAGSequencesDirFmt(
            path=self.get_data_path("mags"),
            mode="r",
        )
        self.feature_data_mags = MAGSequencesDirFmt(
            path=self.get_data_path("mags/sample1"),
            mode="r",
        )
        self.df1 = pd.DataFrame(
            {
                "sample_id": ["sample1"] * 6 + ["sample2"] * 4 + ["sample3"] * 5,
                "mag_id": [f"mag{i}" for i in range(1, 16)],
                "value": range(15),
            }
        )
        self.df2 = pd.DataFrame(
            {
                "sample_id": ["sample1"] * 6 + ["sample2"] * 6 + ["sample3"] * 3,
                "mag_id": [f"mag{i}" for i in range(1, 16)],
                "value": range(15),
            }
        )
        self.df3 = pd.DataFrame(
            {
                "mag_id": ["mag1", "mag2", "mag3"],
                "sample_id": ["sample1", "sample2", "sample3"],
                "dataset": ["dataset1", "dataset2", "dataset3"],
                "single": [1, 2, 3],
                "duplicated": [4, 5, 6],
                "fragmented": [7, 8, 9],
                "missing": [10, 11, 12],
                "complete": [13, 14, 15],
                "n_markers": [16, 17, 18],
                "contigs_n50": [19, 20, 21],
                "percent_gaps": [22, 23, 24],
                "scaffolds": [25, 26, 27],
                "length": [28, 29, 30],
                "completeness": [31, 32, 33],
                "contamination": [34, 35, 36],
                "unbinned_contigs": [10, 11, 12],
            }
        )
        self.df4 = pd.DataFrame(
            {
                "id": ["mag1", "mag2", "mag3"],
                "percent_gaps": ["10%", "20%", "30%"],
                "single": ["1.0", "2.0", "3.0"],
                "duplicated": ["4.0", "5.0", "6.0"],
                "fragmented": ["7.0", "8.0", "9.0"],
                "missing": ["10.0", "11.0", "12.0"],
                "complete": ["13.0", "14.0", "15.0"],
                "n_markers": ["16", "17", "18"],
                "completeness": ["31.0", "32.0", "33.0"],
                "contamination": ["34.0", "35.0", "36.0"],
            }
        )
        self.df5 = pd.DataFrame(
            {
                "index": [0, 1, 2],
                "mag_id": ["mag1", "mag2", "mag3"],
                "percent_gaps": [10.0, 20.0, 30.0],
                "single": [1.0, 2.0, 3.0],
                "duplicated": [4.0, 5.0, 6.0],
                "fragmented": [7.0, 8.0, 9.0],
                "missing": [10.0, 11.0, 12.0],
                "complete": [13.0, 14.0, 15.0],
                "n_markers": [16, 17, 18],
                "completeness": [31.0, 32.0, 33.0],
                "contamination": [34.0, 35.0, 36.0],
            }
        )
        self.busco_results = {
            "dataset": "bacteria_odb10",
            "complete": 87.9,
            "complete_value": 80,
            "single": 86.3,
            "duplicated": 1.6,
            "duplicated_value": 5,
            "fragmented": 4.8,
            "missing": 7.3,
            "missing_value": 10,
            "n_markers": 100,
            "scaffold_n50": "975",
            "contigs_n50": "975",
            "percent_gaps": "0.000%",
            "scaffolds": "2",
            "length": "1935",
        }

    def test_parse_busco_params_1(self):
        observed = _parse_busco_params("evalue", 0.66)
        expected = ["--evalue", str(0.66)]
        self.assertSetEqual(set(observed), set(expected))

    def test_parse_busco_params_2(self):
        observed = _parse_busco_params("augustus", True)
        expected = ["--augustus"]
        self.assertSetEqual(set(observed), set(expected))

    def test_parse_busco_params_3(self):
        observed = _parse_busco_params("lineage_dataset", "bacteria-XYZ")
        expected = ["--lineage_dataset", "bacteria-XYZ"]
        self.assertSetEqual(set(observed), set(expected))

    def test_parse_df_columns(self):
        obs = _parse_df_columns(self.df4)
        exp = self.df5
        pd.testing.assert_frame_equal(obs, exp)

    def test_parse_df_columns_no_additional_metrics(self):
        self.df4.drop(columns=["completeness", "contamination"], inplace=True)
        obs = _parse_df_columns(self.df4)
        exp = self.df5
        exp = exp.drop(columns=["completeness", "contamination"])
        pd.testing.assert_frame_equal(obs, exp)

    def test_partition_dataframe_sample_data_max_rows_5(self):
        partitions = _partition_dataframe(self.df1, 5, True)
        self.assertEqual(len(partitions), 3)
        obs_shapes = [p.shape for p in partitions]
        exp_shapes = [(6, 3), (4, 3), (5, 3)]
        self.assertListEqual(obs_shapes, exp_shapes)

        partitions = _partition_dataframe(self.df2, 5, True)
        self.assertEqual(len(partitions), 3)
        obs_shapes = [p.shape for p in partitions]
        exp_shapes = [(6, 3), (6, 3), (3, 3)]
        self.assertListEqual(obs_shapes, exp_shapes)

    def test_partition_dataframe_sample_data_max_rows_10(self):
        partitions = _partition_dataframe(self.df1, 10, True)
        self.assertEqual(len(partitions), 2)
        obs_shapes = [p.shape for p in partitions]
        exp_shapes = [(10, 3), (5, 3)]
        self.assertListEqual(obs_shapes, exp_shapes)

        partitions = _partition_dataframe(self.df2, 10, True)
        self.assertEqual(len(partitions), 2)
        obs_shapes = [p.shape for p in partitions]
        exp_shapes = [(6, 3), (9, 3)]
        self.assertListEqual(obs_shapes, exp_shapes)

    def test_partition_dataframe_sample_data_max_rows_15(self):
        partitions = _partition_dataframe(self.df1, 15, True)
        self.assertEqual(len(partitions), 1)
        obs_shapes = [p.shape for p in partitions]
        exp_shapes = [
            (15, 3),
        ]
        self.assertListEqual(obs_shapes, exp_shapes)

        partitions = _partition_dataframe(self.df2, 15, True)
        self.assertEqual(len(partitions), 1)
        obs_shapes = [p.shape for p in partitions]
        exp_shapes = [
            (15, 3),
        ]
        self.assertListEqual(obs_shapes, exp_shapes)

    def test_partition_dataframe_feature_data_max_rows_5(self):
        n = 5
        df1 = self.df1.copy()
        df1 = df1.loc[df1["sample_id"] == "sample1"]
        partitions = _partition_dataframe(df1, n, False)
        self.assertEqual(len(partitions), 2)
        obs_shapes = [p.shape for p in partitions]
        exp_shapes = [(5, 3), (1, 3)]
        self.assertListEqual(obs_shapes, exp_shapes)

        df2 = self.df2.copy()
        df2 = df2.loc[df2["sample_id"] == "sample3"]
        partitions = _partition_dataframe(df2, n, False)
        self.assertEqual(len(partitions), 1)
        obs_shapes = [p.shape for p in partitions]
        exp_shapes = [(3, 3)]
        self.assertListEqual(obs_shapes, exp_shapes)

    def test_partition_dataframe_feature_data_max_rows_10(self):
        n = 10
        df1 = self.df1.copy()
        df1 = df1.loc[df1["sample_id"] == "sample1"]
        partitions = _partition_dataframe(df1, n, False)
        self.assertEqual(len(partitions), 1)
        obs_shapes = [p.shape for p in partitions]
        exp_shapes = [(6, 3)]
        self.assertListEqual(obs_shapes, exp_shapes)

        df2 = self.df2.copy()
        df2 = df2.loc[df2["sample_id"] == "sample2"]
        partitions = _partition_dataframe(df2, n, False)
        self.assertEqual(len(partitions), 1)
        obs_shapes = [p.shape for p in partitions]
        exp_shapes = [(6, 3)]
        self.assertListEqual(obs_shapes, exp_shapes)

    def test_partition_dataframe_feature_data_max_rows_15(self):
        n = 10
        df1 = self.df1.copy()
        df1 = df1.loc[df1["sample_id"] == "sample1"]
        partitions = _partition_dataframe(df1, n, False)
        self.assertEqual(len(partitions), 1)
        obs_shapes = [p.shape for p in partitions]
        exp_shapes = [(6, 3)]
        self.assertListEqual(obs_shapes, exp_shapes)

        df2 = self.df2.copy()
        df2 = df2.loc[df2["sample_id"] == "sample2"]
        partitions = _partition_dataframe(df2, n, False)
        self.assertEqual(len(partitions), 1)
        obs_shapes = [p.shape for p in partitions]
        exp_shapes = [(6, 3)]
        self.assertListEqual(obs_shapes, exp_shapes)

    def test_get_feature_table_sample_data(self):
        obs = json.loads(_get_feature_table(self.df3))
        with open(self.get_data_path("feature_table_sample_data.json"), "r") as f:
            exp = json.load(f)
        self.assertDictEqual(obs, exp)

    def test_get_feature_table_feature_data(self):
        df3 = self.df3.copy()
        df3 = df3.loc[df3["sample_id"] == "sample1"]
        obs = json.loads(_get_feature_table(df3))
        with open(self.get_data_path("feature_table_feature_data.json"), "r") as f:
            exp = json.load(f)
        self.assertDictEqual(obs, exp)

    def test_calculate_summary_stats(self):
        obs = _calculate_summary_stats(self.df3)
        exp = pd.DataFrame(
            {
                "min": pd.Series(
                    {
                        "single": 1,
                        "duplicated": 4,
                        "fragmented": 7,
                        "missing": 10,
                        "complete": 13,
                        "completeness": 31.0,
                        "contamination": 34.0,
                        "unbinned_contigs": 10.0,
                    }
                ),
                "median": pd.Series(
                    {
                        "single": 2.0,
                        "duplicated": 5.0,
                        "fragmented": 8.0,
                        "missing": 11.0,
                        "complete": 14.0,
                        "completeness": 32.0,
                        "contamination": 35.0,
                        "unbinned_contigs": 11.0,
                    }
                ),
                "mean": pd.Series(
                    {
                        "single": 2.0,
                        "duplicated": 5.0,
                        "fragmented": 8.0,
                        "missing": 11.0,
                        "complete": 14.0,
                        "completeness": 32.0,
                        "contamination": 35.0,
                        "unbinned_contigs": 11.0,
                    }
                ),
                "max": pd.Series(
                    {
                        "single": 3,
                        "duplicated": 6,
                        "fragmented": 9,
                        "missing": 12,
                        "complete": 15,
                        "completeness": 33.0,
                        "contamination": 36.0,
                        "unbinned_contigs": 12.0,
                    }
                ),
                "count": pd.Series(
                    {
                        "single": 3,
                        "duplicated": 3,
                        "fragmented": 3,
                        "missing": 3,
                        "complete": 3,
                        "completeness": 3.0,
                        "contamination": 3.0,
                        "unbinned_contigs": 3.0,
                    }
                ),
            }
        ).T.to_json(orient="table")

        self.assertEqual(obs, exp)

    def test_calculate_summary_stats_no_additional_metrics(self):
        self.df3.drop(
            columns=["completeness", "contamination", "unbinned_contigs"], inplace=True
        )
        obs = _calculate_summary_stats(self.df3)
        exp = pd.DataFrame(
            {
                "min": pd.Series(
                    {
                        "single": 1,
                        "duplicated": 4,
                        "fragmented": 7,
                        "missing": 10,
                        "complete": 13,
                    }
                ),
                "median": pd.Series(
                    {
                        "single": 2.0,
                        "duplicated": 5.0,
                        "fragmented": 8.0,
                        "missing": 11.0,
                        "complete": 14.0,
                    }
                ),
                "mean": pd.Series(
                    {
                        "single": 2.0,
                        "duplicated": 5.0,
                        "fragmented": 8.0,
                        "missing": 11.0,
                        "complete": 14.0,
                    }
                ),
                "max": pd.Series(
                    {
                        "single": 3,
                        "duplicated": 6,
                        "fragmented": 9,
                        "missing": 12,
                        "complete": 15,
                    }
                ),
                "count": pd.Series(
                    {
                        "single": 3,
                        "duplicated": 3,
                        "fragmented": 3,
                        "missing": 3,
                        "complete": 3,
                    }
                ),
            }
        ).T.to_json(orient="table")

        self.assertEqual(obs, exp)

    def test_validate_lineage_dataset_input_valid(self):
        # Give path to valid database
        p = self.get_data_path("busco_db")
        busco_db = BuscoDatabaseDirFmt(path=p, mode="r")
        _validate_lineage_dataset_input(
            lineage_dataset="lineage_1",
            busco_db=busco_db,
        )

    def test_validate_lineage_dataset_input_invalid(self):
        # Give path to valid database
        p = self.get_data_path("busco_db")
        busco_db = BuscoDatabaseDirFmt(path=p, mode="r")

        with self.assertRaisesRegex(ValueError, "is not present in input database."):
            # Run busco
            _validate_lineage_dataset_input(
                lineage_dataset="lineage2",
                busco_db=busco_db,
            )

    def test_extract_json_data(self):
        obs = _extract_json_data(
            self.get_data_path(
                "busco_output/short_summary.specific.iridoviridae_odb10."
                "24dee6fe-9b84-45bb-8145-de7b092533a1.fasta.json"
            )
        )

        self.assertEqual(obs, self.busco_results)

    def test_calculate_contamination_completeness_normal(self):
        completeness, contamination = _calculate_contamination_completeness(
            missing=10, total=100, duplicated=5, complete=50
        )
        self.assertEqual(completeness, 90.0)
        self.assertEqual(contamination, 10.0)

    def test_calculate_contamination_completeness_divide_by_zero(self):
        completeness, contamination = _calculate_contamination_completeness(
            missing=5, total=100, duplicated=5, complete=0
        )
        self.assertEqual(completeness, 95.0)
        self.assertEqual(contamination, None)

    @patch(
        "q2_annotate.busco.utils._calculate_contamination_completeness",
        return_value=(95.0, 10.0),
    )
    def test_with_completeness_contamination(self, mock_calc):
        results = self.busco_results.copy()
        output = _process_busco_results(
            results=results,
            sample_id="sample1",
            mag_id="mag1",
            file_name="mag1.fasta",
            additional_metrics=True,
        )

        expected = {
            "mag_id": "mag1",
            "sample_id": "sample1",
            "input_file": "mag1.fasta",
            "dataset": "bacteria_odb10",
            "complete": 87.9,
            "single": 86.3,
            "duplicated": 1.6,
            "fragmented": 4.8,
            "missing": 7.3,
            "n_markers": 100,
            "scaffold_n50": "975",
            "contigs_n50": "975",
            "percent_gaps": "0.000%",
            "scaffolds": "2",
            "length": "1935",
            "completeness": 95.0,
            "contamination": 10.0,
        }

        self.assertEqual(output, expected)

    def test_without_completeness_contamination(self):
        results = self.busco_results.copy()
        output = _process_busco_results(
            results=results,
            sample_id="sample1",
            mag_id="mag1",
            file_name="mag1.fasta",
            additional_metrics=False,
        )

        expected = {
            "mag_id": "mag1",
            "sample_id": "sample1",
            "input_file": "mag1.fasta",
            "dataset": "bacteria_odb10",
            "complete": 87.9,
            "single": 86.3,
            "duplicated": 1.6,
            "fragmented": 4.8,
            "missing": 7.3,
            "n_markers": 100,
            "scaffold_n50": "975",
            "contigs_n50": "975",
            "percent_gaps": "0.000%",
            "scaffolds": "2",
            "length": "1935",
        }

        self.assertEqual(output, expected)

    def test_count_binned_contigs(self):
        sample_path = Path(self.get_data_path("mags")) / "sample1"
        fasta_files = glob.glob(os.path.join(sample_path, "*.fasta"))
        count = _count_contigs([Path(x) for x in fasta_files])
        self.assertEqual(count, 7)

    def test_count_unbinned_contigs(self):
        sample_path = Path(self.get_data_path("unbinned")) / "sample1_contigs.fa"
        count = _count_contigs([sample_path])
        self.assertEqual(count, 3)

    def test_calculate_unbinned_percentage(self):
        mag_sample_path = Path(self.get_data_path("mags")) / "sample1"
        mag_sample_files = glob.glob(os.path.join(mag_sample_path, "*.fasta"))
        unbinned_sample_path = (
            Path(self.get_data_path("unbinned")) / "sample1_contigs.fa"
        )
        percentage, count = _calculate_unbinned_percentage(
            [Path(x) for x in mag_sample_files], [unbinned_sample_path]
        )

        # Type checks
        self.assertIsInstance(percentage, float)
        self.assertIsInstance(count, int)

        expected_count = 3
        expected_percentage = (3 / (3 + 7)) * 100

        self.assertEqual(count, expected_count)
        self.assertEqual(percentage, expected_percentage)

    def test_no_unbinned(self):
        mag_sample_path = Path(self.get_data_path("mags")) / "sample1"
        mag_sample_files = glob.glob(os.path.join(mag_sample_path, "*.fasta"))
        unbinned_sample_path = (
            Path(self.get_data_path("unbinned_empty")) / "sample1_contigs.fa"
        )

        percentage, count = _calculate_unbinned_percentage(
            [Path(x) for x in mag_sample_files], [unbinned_sample_path]
        )
        # Type and range checks
        self.assertIsInstance(percentage, float)
        self.assertIsInstance(count, int)
        self.assertEqual(count, 0)
        self.assertEqual(percentage, 0.0)

    def test_only_unbinned(self):
        mag_sample_path = Path(self.get_data_path("mags_empty")) / "sample1"
        mag_sample_files = glob.glob(os.path.join(mag_sample_path, "*.fasta"))
        unbinned_sample_path = (
            Path(self.get_data_path("unbinned")) / "sample1_contigs.fa"
        )

        percentage, count = _calculate_unbinned_percentage(
            [Path(x) for x in mag_sample_files], [unbinned_sample_path]
        )
        expected_count = 3
        self.assertEqual(count, expected_count)
        self.assertEqual(percentage, 100)

    @patch("q2_annotate.busco.utils._count_contigs")
    def test_add_unbinned_metrics(self, mock_count_contigs):
        # Set up test data with scaffolds column
        df = pd.DataFrame(
            {
                "sample_id": ["sample1", "sample1"],
                "mag_id": ["mag1", "mag2"],
                "scaffolds": [10, 20],  # 30 total binned contigs for sample1
                "busco_score": [95.0, 90.0],
            }
        )

        # Mock unbinned_contigs
        unbinned_mock = MagicMock()
        unbinned_mock.sample_dict.return_value = {"sample1": "fake_unbinned.fasta"}

        # Mock _count_contigs to return 5 unbinned contigs
        mock_count_contigs.return_value = 5

        # Call the function
        result = _add_unbinned_metrics(df, unbinned_mock)

        # Verify unbinned_contigs.sample_dict was called
        unbinned_mock.sample_dict.assert_called()

        # Check that the columns were added
        self.assertIn("unbinned_contigs", result.columns)
        self.assertIn("unbinned_contigs_count", result.columns)

        # Check the values - should be same for all MAGs in sample1
        # 5 unbinned / (30 binned + 5 unbinned) * 100 = ~14.29%
        expected_percentage = (5 / (30 + 5)) * 100
        for _, row in result.iterrows():
            self.assertAlmostEqual(
                row["unbinned_contigs"], expected_percentage, places=2
            )
            self.assertEqual(row["unbinned_contigs_count"], 5)

    @patch("q2_annotate.busco.utils._count_contigs")
    def test_add_unbinned_metrics_missing_sample_warning(self, mock_count_contigs):
        # BUSCO results with sample1 and sample2
        df = pd.DataFrame(
            {
                "sample_id": ["sample1", "sample1", "sample2"],
                "mag_id": ["mag1", "mag2", "mag3"],
                "scaffolds": [10, 20, 15],
            }
        )

        # Unbinned contigs only for sample2 (sample1 is missing)
        unbinned_mock = MagicMock()
        unbinned_mock.sample_dict.return_value = {"sample2": "fake_unbinned.fasta"}
        mock_count_contigs.return_value = 5

        # Should warn about missing sample1
        with self.assertWarns(Warning):
            result = _add_unbinned_metrics(df, unbinned_mock)

        # sample1 should have NaN for unbinned metrics
        sample1_rows = result[result["sample_id"] == "sample1"]
        for _, row in sample1_rows.iterrows():
            self.assertTrue(pd.isna(row["unbinned_contigs"]))
            self.assertTrue(pd.isna(row["unbinned_contigs_count"]))

        # sample2 should have proper values
        sample2_rows = result[result["sample_id"] == "sample2"]
        expected_percentage = (5 / (15 + 5)) * 100
        for _, row in sample2_rows.iterrows():
            self.assertAlmostEqual(
                row["unbinned_contigs"], expected_percentage, places=2
            )
            self.assertEqual(row["unbinned_contigs_count"], 5)

    @patch("q2_annotate.busco.utils._count_contigs")
    def test_add_unbinned_metrics_multiple_samples(self, mock_count_contigs):
        # Multiple samples
        df = pd.DataFrame(
            {
                "sample_id": ["sample1", "sample1", "sample2", "sample2"],
                "mag_id": ["mag1", "mag2", "mag3", "mag4"],
                "scaffolds": [10, 20, 15, 25],  # sample1: 30, sample2: 40
            }
        )

        unbinned_mock = MagicMock()
        unbinned_mock.sample_dict.return_value = {
            "sample1": "fake_unbinned_1.fasta",
            "sample2": "fake_unbinned_2.fasta",
        }

        # Return different counts for different samples
        mock_count_contigs.side_effect = [5, 10]  # 5 for sample1, 10 for sample2

        result = _add_unbinned_metrics(df, unbinned_mock)

        # Check sample1 - 5 unbinned / (30 + 5) * 100
        sample1_expected = (5 / 35) * 100
        sample1_rows = result[result["sample_id"] == "sample1"]
        for _, row in sample1_rows.iterrows():
            self.assertAlmostEqual(row["unbinned_contigs"], sample1_expected, places=2)
            self.assertEqual(row["unbinned_contigs_count"], 5)

        # Check sample2 - 10 unbinned / (40 + 10) * 100
        sample2_expected = (10 / 50) * 100
        sample2_rows = result[result["sample_id"] == "sample2"]
        for _, row in sample2_rows.iterrows():
            self.assertAlmostEqual(row["unbinned_contigs"], sample2_expected, places=2)
            self.assertEqual(row["unbinned_contigs_count"], 10)

    @patch("q2_annotate.busco.utils._count_contigs")
    def test_add_unbinned_metrics_zero_unbinned(self, mock_count_contigs):
        df = pd.DataFrame(
            {
                "sample_id": ["sample1"],
                "mag_id": ["mag1"],
                "scaffolds": [30],
            }
        )

        unbinned_mock = MagicMock()
        unbinned_mock.sample_dict.return_value = {"sample1": "fake_unbinned.fasta"}
        mock_count_contigs.return_value = 0  # No unbinned contigs

        result = _add_unbinned_metrics(df, unbinned_mock)

        # Should have 0% unbinned
        self.assertEqual(result["unbinned_contigs"].iloc[0], 0.0)
        self.assertEqual(result["unbinned_contigs_count"].iloc[0], 0)
