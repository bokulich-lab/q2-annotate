# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json
import tempfile
import pandas as pd
from qiime2.plugin.testing import TestPluginBase

from q2_annotate.busco.plots_detailed import _draw_detailed_plots
from q2_annotate.busco.plots_summary import (
    _draw_marker_summary_histograms,
    _draw_selectable_summary_histograms,
    _draw_completeness_vs_contamination,
)


class TestBUSCOPlots(TestPluginBase):
    package = "q2_annotate.busco.tests"

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.df_sample_data = pd.read_csv(
            self.get_data_path("summaries/all_renamed_with_lengths.csv")
        )
        self.df_feature_data = pd.read_csv(
            self.get_data_path("summaries/all_renamed_with_lengths_feature_data.csv")
        )

    def test_draw_detailed_plots_sample_data(self):
        obs = _draw_detailed_plots(
            df=self.df_sample_data,
            is_sample_data=True,
            width=100,
            height=250,
            label_font_size=10,
            title_font_size=15,
            spacing=5,
        )

        self.assertIsInstance(obs, dict)
        self.assertIn("config", obs)
        self.assertIn("hconcat", obs)
        self.assertIn("data", obs["hconcat"][0])
        for i in (0, 1):
            self.assertEqual(obs["hconcat"][i]["spacing"], 5)
            self.assertEqual(obs["hconcat"][i]["spec"]["height"]["step"], 250)
            self.assertEqual(obs["hconcat"][i]["spec"]["width"], 100)

        config = obs["config"]
        for key in ["axis", "header", "legend"]:
            self.assertEqual(config[key]["labelFontSize"], 10)
            self.assertEqual(config[key]["titleFontSize"], 15)

        facet_row = obs["hconcat"][0]["facet"]["row"]
        self.assertEqual(facet_row["title"], "Sample ID and MAG ID")
        self.assertEqual(facet_row["header"]["labelFontSize"], 15)

    def test_draw_detailed_plots_feature_data(self):
        obs = _draw_detailed_plots(
            df=self.df_feature_data,
            is_sample_data=False,
            width=100,
            height=250,
            label_font_size=10,
            title_font_size=15,
            spacing=5,
        )

        self.assertIsInstance(obs, dict)
        self.assertIn("config", obs)
        self.assertIn("hconcat", obs)
        self.assertIn("data", obs["hconcat"][0])
        for i in (0, 1):
            self.assertEqual(obs["hconcat"][i]["spacing"], 5)
            self.assertEqual(obs["hconcat"][i]["spec"]["height"]["step"], 250)
            self.assertEqual(obs["hconcat"][i]["spec"]["width"], 100)

        config = obs["config"]
        for key in ["axis", "header", "legend"]:
            self.assertEqual(config[key]["labelFontSize"], 10)
            self.assertEqual(config[key]["titleFontSize"], 15)

        facet_row = obs["hconcat"][0]["facet"]["row"]
        self.assertEqual(facet_row["title"], "MAG ID")
        self.assertEqual(facet_row["header"]["labelFontSize"], 0)

    def test_draw_marker_summary_histograms(self):
        obs = _draw_marker_summary_histograms(data=self.df_sample_data)

        self.assertIsInstance(obs, dict)
        self.assertIn("config", obs)
        self.assertIn("vconcat", obs)
        self.assertEqual(len(obs["vconcat"][0]["hconcat"]), 5)
        self.assertEqual(len(obs["vconcat"][1]["hconcat"]), 5)

        exp_titles = ["Single", "Duplicated", "Fragmented", "Missing", "Completeness"]
        obs_titles = [x["encoding"]["x"]["title"] for x in obs["vconcat"][0]["hconcat"]]
        self.assertListEqual(exp_titles, obs_titles)

        exp_titles = [
            "Contamination",
            "Scaffolds",
            "Contigs n50",
            "Scaffold n50",
            "Length",
        ]
        obs_titles = [x["encoding"]["x"]["title"] for x in obs["vconcat"][1]["hconcat"]]
        self.assertListEqual(exp_titles, obs_titles)

        config = obs["config"]
        for key in ["axis", "header", "legend"]:
            self.assertEqual(config[key]["labelFontSize"], 12)
            self.assertEqual(config[key]["titleFontSize"], 15)

    def test_draw_marker_summary_histograms_no_additional_metrics(self):
        self.df_sample_data.drop(
            columns=["completeness", "contamination"], inplace=True
        )
        obs = _draw_marker_summary_histograms(data=self.df_sample_data)

        self.assertIsInstance(obs, dict)
        self.assertIn("config", obs)
        self.assertIn("vconcat", obs)
        self.assertEqual(len(obs["vconcat"][0]["hconcat"]), 4)
        self.assertEqual(len(obs["vconcat"][1]["hconcat"]), 4)

        exp_titles = ["Single", "Duplicated", "Fragmented", "Missing"]
        obs_titles = [x["encoding"]["x"]["title"] for x in obs["vconcat"][0]["hconcat"]]
        self.assertListEqual(exp_titles, obs_titles)

        exp_titles = ["Scaffolds", "Contigs n50", "Scaffold n50", "Length"]
        obs_titles = [x["encoding"]["x"]["title"] for x in obs["vconcat"][1]["hconcat"]]
        self.assertListEqual(exp_titles, obs_titles)

        config = obs["config"]
        for key in ["axis", "header", "legend"]:
            self.assertEqual(config[key]["labelFontSize"], 12)
            self.assertEqual(config[key]["titleFontSize"], 15)

    def test_draw_selectable_summary_histograms(self):
        obs = _draw_selectable_summary_histograms(data=self.df_sample_data)

        self.assertIsInstance(obs, dict)
        self.assertIn("config", obs)
        self.assertIn("mark", obs)
        self.assertIn("transform", obs)
        self.assertIn("filter", obs["transform"][0])
        self.assertEqual(obs["mark"]["type"], "bar")

        config = obs["config"]
        for key in ["axis", "header", "legend"]:
            self.assertEqual(config[key]["labelFontSize"], 12)
            self.assertEqual(config[key]["titleFontSize"], 15)

        cat = {entry["category"] for entry in obs["datasets"][obs["data"]["name"]]}
        self.assertTrue("completeness" in cat)
        self.assertTrue("contamination" in cat)

    def test_draw_selectable_summary_histograms_no_additional_metrics(self):
        self.df_sample_data.drop(
            columns=["completeness", "contamination"], inplace=True
        )
        obs = _draw_selectable_summary_histograms(data=self.df_sample_data)

        self.assertIsInstance(obs, dict)
        self.assertIn("config", obs)
        self.assertIn("mark", obs)
        self.assertIn("transform", obs)
        self.assertIn("filter", obs["transform"][0])
        self.assertEqual(obs["mark"]["type"], "bar")

        config = obs["config"]
        for key in ["axis", "header", "legend"]:
            self.assertEqual(config[key]["labelFontSize"], 12)
            self.assertEqual(config[key]["titleFontSize"], 15)

        cat = {entry["category"] for entry in obs["datasets"][obs["data"]["name"]]}
        self.assertFalse("completeness" in cat)
        self.assertFalse("contamination" in cat)

    def test_scatter_sample_data(self):
        obs = _draw_completeness_vs_contamination(self.df_sample_data)
        # altair .interactive() adds parameter name that increments with each call
        obs["params"][1]["name"] = "param_1"
        with open(self.get_data_path("scatterplot/sample_data.json"), "r") as f:
            exp = json.load(f)
        self.assertEqual(obs, exp)

    def test_scatter_feature_data(self):
        obs = _draw_completeness_vs_contamination(self.df_feature_data)
        # altair .interactive() adds parameter name that increments with each call
        obs["params"][1]["name"] = "param_1"
        with open(self.get_data_path("scatterplot/feature_data.json"), "r") as f:
            exp = json.load(f)
        self.assertEqual(obs, exp)
