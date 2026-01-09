# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json
import pandas as pd
from qiime2.plugin.testing import TestPluginBase

from q2_annotate.busco.busco import (
    _prepare_histogram_data,
    _prepare_box_plot_data,
    _prepare_scatter_data,
    _prepare_detailed_data,
    _prepare_assembly_data,
)


class TestVegaDataPreparation(TestPluginBase):
    package = "q2_annotate.busco.tests"

    def setUp(self):
        super().setUp()
        # Basic test data with all metrics
        self.full_df = pd.DataFrame(
            {
                "sample_id": ["sample1", "sample1", "sample2"],
                "mag_id": ["mag1", "mag2", "mag3"],
                "dataset": ["bacteria_odb10", "bacteria_odb10", "bacteria_odb10"],
                "n_markers": [100, 100, 100],
                "single": [85.0, 90.0, 80.0],
                "duplicated": [5.0, 3.0, 7.0],
                "fragmented": [3.0, 2.0, 4.0],
                "missing": [7.0, 5.0, 9.0],
                "completeness": [93.0, 95.0, 91.0],
                "contamination": [5.9, 3.3, 8.8],
                "contigs_n50": [50000, 60000, 45000],
                "length": [2000000, 2500000, 1800000],
                "scaffold_n50": [51000, 61000, 46000],
                "percent_gaps": [0.5, 0.3, 0.7],
                "scaffolds": [50, 40, 60],
            }
        )

        # Test data without completeness/contamination
        self.basic_df = self.full_df.drop(columns=["completeness", "contamination"])

    def test_prepare_histogram_data_with_all_metrics(self):
        result = _prepare_histogram_data(self.full_df)

        # Parse JSON result
        data = json.loads(result)

        # Should have melted all metrics for all MAGs
        # 3 MAGs * 8 metrics = 24 records
        self.assertEqual(len(data), 24)

        # Check structure
        for record in data:
            self.assertIn("sample_id", record)
            self.assertIn("mag_id", record)
            self.assertIn("dataset", record)
            self.assertIn("n_markers", record)
            self.assertIn("category", record)
            self.assertIn("metric", record)

        # Check that all expected metrics are present
        categories = {record["category"] for record in data}
        expected = {
            "single",
            "duplicated",
            "fragmented",
            "missing",
            "completeness",
            "contamination",
            "contigs_n50",
            "length",
        }
        self.assertEqual(categories, expected)

    def test_prepare_histogram_data_without_completeness(self):
        result = _prepare_histogram_data(self.basic_df)

        data = json.loads(result)

        # 3 MAGs * 6 metrics (no completeness/contamination) = 18 records
        self.assertEqual(len(data), 18)

        categories = {record["category"] for record in data}
        expected = {
            "single",
            "duplicated",
            "fragmented",
            "missing",
            "contigs_n50",
            "length",
        }
        self.assertEqual(categories, expected)

        # Ensure completeness and contamination are not present
        self.assertNotIn("completeness", categories)
        self.assertNotIn("contamination", categories)

    def test_prepare_box_plot_data_with_all_metrics(self):
        result = _prepare_box_plot_data(self.full_df)

        # Should be a dict with metric names as keys
        self.assertIsInstance(result, dict)

        # Check all expected metrics are present
        expected_metrics = [
            "single",
            "duplicated",
            "fragmented",
            "missing",
            "completeness",
            "contamination",
            "contigs_n50",
            "length",
        ]
        for metric in expected_metrics:
            self.assertIn(metric, result)

        # Check structure of each metric's data
        for metric, data in result.items():
            self.assertIsInstance(data, list)
            self.assertEqual(len(data), 3)  # 3 MAGs
            for record in data:
                self.assertIn("sample_id", record)
                self.assertIn("mag_id", record)
                self.assertIn("value", record)

    def test_prepare_box_plot_data_without_completeness(self):
        result = _prepare_box_plot_data(self.basic_df)

        # Should not have completeness/contamination
        self.assertNotIn("completeness", result)
        self.assertNotIn("contamination", result)

        # Should have other metrics
        self.assertIn("single", result)
        self.assertIn("contigs_n50", result)

    def test_prepare_scatter_data_with_metrics(self):
        data_str, has_data, upper_x, upper_y = _prepare_scatter_data(self.full_df)

        # Should have data
        self.assertTrue(has_data)
        self.assertIsNotNone(data_str)

        # Parse JSON
        data = json.loads(data_str)
        self.assertEqual(len(data), 3)  # 3 MAGs

        # Check axis bounds are calculated
        self.assertGreater(upper_x, 0)
        self.assertGreater(upper_y, 0)
        # Upper bounds should be slightly higher than max values (91-95 range)
        self.assertAlmostEqual(upper_x, 95.0 * 1.1, places=1)
        self.assertLessEqual(upper_x, 110.0)
        self.assertAlmostEqual(upper_y, 8.8 * 1.1, places=1)
        self.assertLessEqual(upper_y, 110.0)

    def test_prepare_scatter_data_without_metrics(self):
        data_str, has_data, upper_x, upper_y = _prepare_scatter_data(self.basic_df)

        # Should not have data when completeness/contamination missing
        self.assertFalse(has_data)
        self.assertIsNone(data_str)

        # Default bounds
        self.assertEqual(upper_x, 110)
        self.assertEqual(upper_y, 110)

    def test_prepare_scatter_data_with_extreme_values(self):
        # Test with very low values
        low_df = self.full_df.copy()
        low_df["completeness"] = [2.0, 3.0, 1.0]
        low_df["contamination"] = [1.0, 2.0, 0.5]

        _, _, upper_x, upper_y = _prepare_scatter_data(low_df)

        # Should use minimum bound of 5.0
        self.assertEqual(upper_x, 5.0)
        self.assertEqual(upper_y, 5.0)

    def test_prepare_detailed_data(self):
        result = _prepare_detailed_data(self.full_df)

        data = json.loads(result)

        # 3 MAGs * 4 BUSCO categories = 12 records
        self.assertEqual(len(data), 12)

        # Check structure
        for record in data:
            self.assertIn("sample_id", record)
            self.assertIn("mag_id", record)
            self.assertIn("dataset", record)
            self.assertIn("n_markers", record)
            self.assertIn("category", record)
            self.assertIn("BUSCO_percentage", record)
            self.assertIn("frac_markers", record)

            # Check frac_markers format: "~XX/100"
            self.assertTrue(record["frac_markers"].startswith("~"))
            self.assertIn("/", record["frac_markers"])
            self.assertTrue(record["frac_markers"].endswith("/100"))

    def test_prepare_assembly_data_with_all_columns(self):
        result = _prepare_assembly_data(self.full_df)

        data = json.loads(result)
        self.assertEqual(len(data), 3)  # 3 MAGs

        # Check all expected columns are present
        for record in data:
            self.assertIn("sample_id", record)
            self.assertIn("mag_id", record)
            self.assertIn("scaffold_n50", record)
            self.assertIn("contigs_n50", record)
            self.assertIn("percent_gaps", record)
            self.assertIn("scaffolds", record)

    def test_prepare_assembly_data_with_missing_columns(self):
        # Create df without scaffold_n50
        df = self.full_df.drop(columns=["scaffold_n50"])
        result = _prepare_assembly_data(df)

        data = json.loads(result)

        # Should still work but without scaffold_n50
        for record in data:
            self.assertNotIn("scaffold_n50", record)
            self.assertIn("contigs_n50", record)
            self.assertIn("scaffolds", record)

    def test_nan_handling_in_histogram(self):
        # Add NaN values
        df_with_nan = self.full_df.copy()
        df_with_nan.loc[0, "completeness"] = float("nan")

        result = _prepare_histogram_data(df_with_nan)

        # Should replace NaN with null in JSON
        self.assertIn("null", result)
        self.assertNotIn("NaN", result)

    def test_nan_handling_in_scatter(self):
        # Test scatter with NaN values
        df_with_nan = self.full_df.copy()
        df_with_nan.loc[0, "completeness"] = float("nan")

        data_str, has_data, upper_x, upper_y = _prepare_scatter_data(df_with_nan)

        # Should still work
        self.assertTrue(has_data)
        self.assertIn("null", data_str)
        self.assertNotIn("NaN", data_str)
