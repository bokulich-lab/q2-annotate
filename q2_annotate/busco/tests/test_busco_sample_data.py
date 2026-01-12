# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json

import qiime2
import pandas as pd
from q2_annotate.busco.busco import (
    _run_busco,
    _busco_helper,
    _evaluate_busco,
    _visualize_busco,
    evaluate_busco,
)
from unittest.mock import patch, ANY, call, MagicMock
from qiime2.plugin.testing import TestPluginBase
from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt, ContigSequencesDirFmt
from q2_annotate.busco.types import BuscoDatabaseDirFmt


class TestBUSCOSampleData(TestPluginBase):
    package = "q2_annotate.busco.tests"

    def setUp(self):
        super().setUp()
        self.mags = MultiMAGSequencesDirFmt(
            path=self.get_data_path("mags"),
            mode="r",
        )
        self.unbinned = ContigSequencesDirFmt(
            path=self.get_data_path("unbinned"), mode="r"
        )

        self.busco_db = BuscoDatabaseDirFmt(
            path=self.get_data_path("busco_db"), mode="r"
        )

    @patch("q2_annotate.busco.busco.run_command")
    def test_run_busco(self, mock_run):
        _run_busco(
            input_dir="input_dir",
            output_dir="cwd/output_dir",
            sample_id="sample1",
            params=["--lineage_dataset", "bacteria_odb10", "--cpu", "7"],
        )

        mock_run.assert_called_once_with(
            [
                "busco",
                "--lineage_dataset",
                "bacteria_odb10",
                "--cpu",
                "7",
                "--in",
                "input_dir",
                "--out_path",
                "cwd/output_dir",
                "-o",
                "sample1",
            ],
            cwd="cwd",
        )

    @patch("q2_annotate.busco.busco._extract_json_data")
    @patch("q2_annotate.busco.busco._process_busco_results")
    @patch("q2_annotate.busco.busco._run_busco")
    @patch("q2_annotate.busco.busco.glob.glob")
    def test_busco_helper(self, mock_glob, mock_run, mock_process, mock_extract):
        with open(
            self.get_data_path("busco_results_json/busco_results.json"), "r"
        ) as f:
            busco_list = json.load(f)

        mock_process.side_effect = busco_list

        obs = _busco_helper(self.mags, ["--lineage_dataset", "bacteria_odb10"], True)
        exp = pd.read_csv(
            self.get_data_path("busco_results/results_all/busco_results.tsv"), sep="\t"
        )

        pd.testing.assert_frame_equal(obs, exp)

        mock_run.assert_has_calls(
            [
                call(
                    input_dir=ANY,
                    output_dir=ANY,
                    sample_id="sample1",
                    params=["--lineage_dataset", "bacteria_odb10"],
                ),
                call(
                    input_dir=ANY,
                    output_dir=ANY,
                    sample_id="sample2",
                    params=["--lineage_dataset", "bacteria_odb10"],
                ),
            ]
        )

        mag_ids = [
            "24dee6fe-9b84-45bb-8145-de7b092533a1",
            "ca7012fc-ba65-40c3-84f5-05aa478a7585",
            "fb0bc871-04f6-486b-a10e-8e0cb66f8de3",
            "d65a71fa-4279-4588-b937-0747ed5d604d",
            "db03f8b6-28e1-48c5-a47c-9c65f38f7357",
            "fa4d7420-d0a4-455a-b4d7-4fa66e54c9bf",
        ]
        samples = ["sample1", "sample1", "sample1", "sample2", "sample2", "sample2"]

        expected_calls = [
            call(ANY, sample, mag_id, f"{mag_id}.fasta", True)
            for sample, mag_id in zip(samples, mag_ids)
        ]
        mock_process.assert_has_calls(expected_calls)

    @patch("q2_annotate.busco.busco._busco_helper")
    def test_evaluate_busco_offline(self, mock_helper):
        _evaluate_busco(
            mags=self.mags,
            unbinned_contigs=self.unbinned,  # new
            db=self.busco_db,
            mode="some_mode",
            lineage_dataset="lineage_1",
        )
        mock_helper.assert_called_with(
            self.mags,
            [
                "--mode",
                "some_mode",
                "--lineage_dataset",
                "lineage_1",
                "--cpu",
                "1",
                "--contig_break",
                "10",
                "--evalue",
                "0.001",
                "--limit",
                "3",
                "--offline",
                "--download_path",
                str(self.busco_db),
            ],
            False,
        )

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
                self.get_data_path("summaries/all_renamed_with_lengths.csv")
            ),
        )

        mock_render.assert_called_once()
        context = mock_render.call_args[1]["context"]

        # Check that tabs are correct
        self.assertEqual(len(context["tabs"]), 3)
        self.assertEqual(context["tabs"][0]["title"], "QC overview")
        self.assertEqual(context["tabs"][1]["title"], "Sample details")
        self.assertEqual(context["tabs"][2]["title"], "Feature details")

        # Check that mocked values are passed through
        self.assertEqual(context["table"], "table1")
        self.assertEqual(context["summary_stats_json"], "stats1")

        # Check unbinned flag
        self.assertTrue(context["unbinned"])

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
            "SampleData[MAGs]", self.get_data_path("mags")
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
            if action == "partition_sample_data_mags":
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
            db=busco_db,
            unbinned_contigs=unbinned,
            num_partitions=2,
            lineage_dataset="bacteria_odb10",
        )
        exp = ("collated_result", "visualization")

        self.assertTupleEqual(obs, exp)

    def test_evaluate_busco_action_no_unbinned(self):
        mags = qiime2.Artifact.import_data(
            "SampleData[MAGs]", self.get_data_path("mags")
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
            if action == "partition_sample_data_mags":
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
            db=busco_db,
            unbinned_contigs=None,
            num_partitions=2,
            lineage_dataset="bacteria_odb10",
        )
        exp = ("collated_result", "visualization")

        self.assertTupleEqual(obs, exp)
