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

    def _prep_histogram_expected(self, metrics):
        """Helper to generate expected histogram data for given metrics."""
        base = {
            "sample_id": ["sample1", "sample1", "sample2"],
            "mag_id": ["mag1", "mag2", "mag3"],
            "dataset": "bacteria_odb10",
            "n_markers": 100,
        }
        expected = []
        for metric in metrics:
            for i, mag in enumerate(base["mag_id"]):
                expected.append({
                    "sample_id": base["sample_id"][i],
                    "mag_id": mag,
                    "dataset": base["dataset"],
                    "n_markers": base["n_markers"],
                    "category": metric,
                    "metric": self.full_df.loc[i, metric] if metric in self.full_df.columns else self.basic_df.loc[i, metric],
                })
        return expected

    def _prep_box_plot_expected(self, metrics):
        """Helper to generate expected box plot data for given metrics."""
        base = {
            "sample_id": ["sample1", "sample1", "sample2"],
            "mag_id": ["mag1", "mag2", "mag3"],
        }
        expected = {}
        for metric in metrics:
            expected[metric] = []
            for i, mag in enumerate(base["mag_id"]):
                expected[metric].append({
                    "sample_id": base["sample_id"][i],
                    "mag_id": mag,
                    "value": self.full_df.loc[i, metric] if metric in self.full_df.columns else self.basic_df.loc[i, metric],
                })
        return expected

    def test_prepare_histogram_data_with_all_metrics(self):
        result = _prepare_histogram_data(self.full_df)
        data = json.loads(result)

        metrics = ["single", "duplicated", "fragmented", "missing", "completeness",
                   "contamination", "contigs_n50", "length"]
        expected = self._prep_histogram_expected(metrics)

        self.assertEqual(data, expected)

    def test_prepare_histogram_data_without_completeness(self):
        result = _prepare_histogram_data(self.basic_df)
        data = json.loads(result)

        metrics = ["single", "duplicated", "fragmented", "missing", "contigs_n50", "length"]
        expected = self._prep_histogram_expected(metrics)

        self.assertEqual(data, expected)

    def test_prepare_box_plot_data_with_all_metrics(self):
        result = _prepare_box_plot_data(self.full_df)

        metrics = ["single", "duplicated", "fragmented", "missing", "completeness",
                   "contamination", "contigs_n50", "length"]
        expected = self._prep_box_plot_expected(metrics)

        self.assertEqual(result, expected)

    def test_prepare_box_plot_data_without_completeness(self):
        result = _prepare_box_plot_data(self.basic_df)

        metrics = ["single", "duplicated", "fragmented", "missing", "contigs_n50", "length"]
        expected = self._prep_box_plot_expected(metrics)

        self.assertEqual(result, expected)

    def test_prepare_scatter_data_with_metrics(self):
        data_str, has_data, upper_x, upper_y = _prepare_scatter_data(self.full_df)

        self.assertTrue(has_data)
        data = json.loads(data_str)

        # Function returns all columns, so check structure and key fields
        self.assertEqual(len(data), 3)

        # Check first record has all expected fields
        self.assertIn("sample_id", data[0])
        self.assertIn("mag_id", data[0])
        self.assertIn("completeness", data[0])
        self.assertIn("contamination", data[0])

        # Check specific values for key fields
        self.assertEqual(data[0]["sample_id"], "sample1")
        self.assertEqual(data[0]["mag_id"], "mag1")
        self.assertEqual(data[0]["completeness"], 93.0)
        self.assertEqual(data[0]["contamination"], 5.9)

        self.assertEqual(data[1]["completeness"], 95.0)
        self.assertEqual(data[1]["contamination"], 3.3)

        self.assertEqual(data[2]["completeness"], 91.0)
        self.assertEqual(data[2]["contamination"], 8.8)

        self.assertAlmostEqual(upper_x, 95.0 * 1.1, places=1)
        self.assertAlmostEqual(upper_y, 8.8 * 1.1, places=1)

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

        base = {"sample_id": ["sample1", "sample1", "sample2"],
                "mag_id": ["mag1", "mag2", "mag3"]}
        categories = ["single", "duplicated", "fragmented", "missing"]

        expected = []
        for cat in categories:
            for i, mag in enumerate(base["mag_id"]):
                val = self.full_df.loc[i, cat]
                expected.append({
                    "sample_id": base["sample_id"][i],
                    "mag_id": mag,
                    "dataset": "bacteria_odb10",
                    "n_markers": 100,
                    "category": cat,
                    "BUSCO_percentage": val,
                    "frac_markers": f"~{int(val)}/100",
                })

        self.assertEqual(data, expected)

    def test_prepare_assembly_data_with_all_columns(self):
        result = _prepare_assembly_data(self.full_df)
        data = json.loads(result)

        expected = [
            {"sample_id": "sample1", "mag_id": "mag1", "scaffold_n50": 51000, "contigs_n50": 50000, "percent_gaps": 0.5, "scaffolds": 50},
            {"sample_id": "sample1", "mag_id": "mag2", "scaffold_n50": 61000, "contigs_n50": 60000, "percent_gaps": 0.3, "scaffolds": 40},
            {"sample_id": "sample2", "mag_id": "mag3", "scaffold_n50": 46000, "contigs_n50": 45000, "percent_gaps": 0.7, "scaffolds": 60},
        ]
        self.assertEqual(data, expected)

    def test_prepare_assembly_data_with_missing_columns(self):
        df = self.full_df.drop(columns=["scaffold_n50"])
        result = _prepare_assembly_data(df)
        data = json.loads(result)

        expected = [
            {"sample_id": "sample1", "mag_id": "mag1", "contigs_n50": 50000, "percent_gaps": 0.5, "scaffolds": 50},
            {"sample_id": "sample1", "mag_id": "mag2", "contigs_n50": 60000, "percent_gaps": 0.3, "scaffolds": 40},
            {"sample_id": "sample2", "mag_id": "mag3", "contigs_n50": 45000, "percent_gaps": 0.7, "scaffolds": 60},
        ]
        self.assertEqual(data, expected)

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
