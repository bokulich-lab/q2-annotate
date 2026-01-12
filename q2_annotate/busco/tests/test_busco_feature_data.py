# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json
from unittest.mock import patch, ANY, MagicMock, call

import pandas as pd
import pytest
import qiime2
from q2_types.feature_data_mag import MAGSequencesDirFmt
from qiime2.plugin.testing import TestPluginBase

from q2_annotate.busco.busco import _visualize_busco, evaluate_busco, _busco_helper


class TestBUSCOFeatureData(TestPluginBase):
    package = "q2_annotate.busco.tests"

    def setUp(self):
        super().setUp()
        self.mags = MAGSequencesDirFmt(
            path=self.get_data_path("mags/sample1"),
            mode="r",
        )

    @patch("q2_annotate.busco.busco._extract_json_data")
    @patch("q2_annotate.busco.busco._process_busco_results")
    @patch("q2_annotate.busco.busco._run_busco")
    @patch("q2_annotate.busco.busco.glob.glob")
    def test_busco_helper(self, mock_glob, mock_run, mock_process, mock_extract):
        with open(
            self.get_data_path("busco_results_json/busco_results_feature_data.json"),
            "r",
        ) as f:
            busco_list = json.load(f)

        mock_process.side_effect = busco_list

        obs = _busco_helper(self.mags, ["--lineage_dataset", "bacteria_odb10"], True)

        exp = pd.read_csv(
            self.get_data_path(
                "busco_results/results_all/busco_results_feature_data.tsv"
            ),
            sep="\t",
            keep_default_na=False,
        )
        exp["sample_id"] = exp["sample_id"].astype(object)
        pd.testing.assert_frame_equal(obs, exp)

        mock_run.assert_called_once_with(
            input_dir=ANY,
            output_dir=ANY,
            sample_id="feature_data",
            params=["--lineage_dataset", "bacteria_odb10"],
        )
        mag_ids = [
            "24dee6fe-9b84-45bb-8145-de7b092533a1",
            "ca7012fc-ba65-40c3-84f5-05aa478a7585",
            "fb0bc871-04f6-486b-a10e-8e0cb66f8de3",
        ]
        expected_calls = [
            call(ANY, "feature_data", mag_id, f"{mag_id}.fasta", True)
            for mag_id in mag_ids
        ]
        mock_process.assert_has_calls(expected_calls)

    @patch("q2_annotate.busco.busco._get_feature_table", return_value="table1")
    @patch("q2_annotate.busco.busco._calculate_summary_stats", return_value="stats1")
    @patch("q2templates.render")
    @patch("q2_annotate.busco.busco._cleanup_bootstrap")
    def test_visualize_busco(
        self,
        mock_clean,
        mock_render,
        mock_stats,
        mock_table,
    ):
        _visualize_busco(
            output_dir=self.temp_dir.name,
            results=pd.read_csv(
                self.get_data_path(
                    "summaries/all_renamed_with_lengths_feature_data.csv"
                )
            ),
        )

        mock_render.assert_called_once()
        context = mock_render.call_args[1]["context"]

        # Check that tabs are correct
        self.assertEqual(len(context["tabs"]), 3)
        self.assertEqual(context["tabs"][0]["title"], "QC overview")
        self.assertEqual(context["tabs"][1]["title"], "Per-MAG metrics")
        self.assertEqual(context["tabs"][2]["title"], "Feature details")

        # Check that mocked values are passed through
        self.assertEqual(context["table"], "table1")
        self.assertEqual(context["summary_stats_json"], "stats1")

        # Check unbinned flag
        self.assertFalse(context["unbinned"])

        # Check that all expected keys are present
        expected_keys = {
            "tabs", "vega_histogram_spec", "vega_unbinned_spec",
            "vega_box_plot_spec", "vega_scatter_spec", "vega_busco_detailed_spec",
            "vega_assembly_detailed_spec", "histogram_data", "box_plot_data",
            "scatter_data", "detailed_data", "assembly_data", "mag_ids_sorted",
            "metrics_json", "assembly_metrics_json", "sample_ids_json",
            "is_sample_data", "comp_cont", "upper_x", "upper_y",
            "unbinned", "unbinned_data", "table", "summary_stats_json", "page_size"
        }
        self.assertEqual(set(context.keys()), expected_keys)

        mock_clean.assert_called_with(self.temp_dir.name)

    # TODO: maybe this could be turned into an actual test
    def test_evaluate_busco_action(self):
        mags = qiime2.Artifact.import_data(
            "FeatureData[MAG]", self.get_data_path("mags/sample2")
        )
        busco_db = qiime2.Artifact.import_data(
            "ReferenceDB[BUSCO]", self.get_data_path("busco_db")
        )

        fake_partition = MagicMock()
        fake_partition.values.return_value = ["partition1", "partition2"]

        # Create mock actions
        partition_action_mock = MagicMock(return_value=(fake_partition,))
        evaluate_busco_mock = MagicMock(return_value=(0,))
        collate_mock = MagicMock(return_value=("collated_result",))
        visualize_mock = MagicMock(return_value=("visualization",))

        def get_action_side_effect(plugin, action):
            if action == "partition_feature_data_mags":
                return partition_action_mock
            elif action == "_evaluate_busco":
                return evaluate_busco_mock
            elif action == "collate_busco_results":
                return collate_mock
            elif action == "_visualize_busco":
                return visualize_mock

        mock_ctx = MagicMock()
        mock_ctx.get_action.side_effect = get_action_side_effect

        obs = evaluate_busco(
            ctx=mock_ctx,
            mags=mags,
            unbinned_contigs=None,
            db=busco_db,
            num_partitions=2,
            lineage_dataset="bacteria_odb10",
        )
        exp = ("collated_result", "visualization")
        self.assertTupleEqual(obs, exp)

    def test_evaluate_busco_action_with_unbinned(self):
        mags = qiime2.Artifact.import_data(
            "FeatureData[MAG]", self.get_data_path("mags/sample2")
        )
        unbinned = qiime2.Artifact.import_data(
            "SampleData[Contigs]", self.get_data_path("unbinned")
        )
        busco_db = qiime2.Artifact.import_data(
            "ReferenceDB[BUSCO]", self.get_data_path("busco_db")
        )

        fake_partition = MagicMock()
        fake_partition.values.return_value = ["partition1", "partition2"]

        # Create mock actions
        partition_action_mock = MagicMock(return_value=(fake_partition,))
        evaluate_busco_mock = MagicMock(return_value=(0,))
        collate_mock = MagicMock(return_value=("collated_result",))
        visualize_mock = MagicMock(return_value=("visualization",))

        def get_action_side_effect(plugin, action):
            if action == "partition_feature_data_mags":
                return partition_action_mock
            elif action == "_evaluate_busco":
                return evaluate_busco_mock
            elif action == "collate_busco_results":
                return collate_mock
            elif action == "_visualize_busco":
                return visualize_mock

        mock_ctx = MagicMock()
        mock_ctx.get_action.side_effect = get_action_side_effect

        with pytest.warns(match="unbinned contigs will be ignored"):
            obs = evaluate_busco(
                ctx=mock_ctx,
                mags=mags,
                unbinned_contigs=unbinned,
                db=busco_db,
                num_partitions=2,
                lineage_dataset="bacteria_odb10",
            )
        exp = ("collated_result", "visualization")
        self.assertTupleEqual(obs, exp)
