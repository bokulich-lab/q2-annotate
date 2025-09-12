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
        mock_process.assert_has_calls(
            [
                call(
                    ANY,
                    "feature_data",
                    "24dee6fe-9b84-45bb-8145-de7b092533a1",
                    "24dee6fe-9b84-45bb-8145-de7b092533a1.fasta",
                    True,
                ),
                call(
                    ANY,
                    "feature_data",
                    "ca7012fc-ba65-40c3-84f5-05aa478a7585",
                    "ca7012fc-ba65-40c3-84f5-05aa478a7585.fasta",
                    True,
                ),
                call(
                    ANY,
                    "feature_data",
                    "fb0bc871-04f6-486b-a10e-8e0cb66f8de3",
                    "fb0bc871-04f6-486b-a10e-8e0cb66f8de3.fasta",
                    True,
                ),
            ]
        )

    @patch(
        "q2_annotate.busco.busco._draw_detailed_plots",
        return_value={"fake1": {"plot": "NaN"}},
    )
    @patch(
        "q2_annotate.busco.busco._draw_marker_summary_histograms",
        return_value={"fake2": {"plot": "NaN"}},
    )
    @patch(
        "q2_annotate.busco.busco._draw_selectable_summary_histograms",
        return_value={"fake3": {"plot": "spec"}},
    )
    @patch(
        "q2_annotate.busco.busco._draw_completeness_vs_contamination",
        return_value={"fake4": {"plot": "NaN"}},
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
        mock_scatter,
        mock_selectable,
        mock_marker,
        mock_detailed,
    ):
        _visualize_busco(
            output_dir=self.temp_dir.name,
            results=pd.read_csv(
                self.get_data_path(
                    "summaries/all_renamed_with_lengths_feature_data.csv"
                )
            ),
        )

        mock_detailed.assert_called_once()
        mock_marker.assert_called_once()

        exp_context = {
            "tabs": [
                {"title": "QC overview", "url": "index.html"},
                {"title": "BUSCO plots", "url": "detailed_view.html"},
                {"title": "BUSCO table", "url": "table.html"},
            ],
            "vega_json": json.dumps(
                {
                    "partition_0": {
                        "subcontext": {"fake1": {"plot": "null"}},
                        "counters": {"from": 1, "to": 2},
                        "ids": [
                            "ab23d75d-547d-455a-8b51-16b46ddf7496",
                            "0e514d88-16c4-4273-a1df-1a360eb2c823",
                        ],
                    }
                }
            ),
            "vega_summary_json": json.dumps({"fake2": {"plot": "null"}}),
            "table": "table1",
            "summary_stats_json": "stats1",
            "scatter_json": json.dumps({"fake4": {"plot": "null"}}),
            "comp_cont": True,
            "unbinned": False,
            "vega_selectable_unbinned_json": None,
            "page_size": 100,
        }
        mock_render.assert_called_with(ANY, self.temp_dir.name, context=exp_context)
        mock_clean.assert_called_with(self.temp_dir.name)

    # TODO: maybe this could be turned into an actual test
    @patch("q2_annotate.busco.busco._validate_parameters")
    def test_evaluate_busco_action(self, mock_validate):
        mags = qiime2.Artifact.import_data(
            "FeatureData[MAG]", self.get_data_path("mags/sample2")
        )
        busco_db = qiime2.Artifact.import_data(
            "ReferenceDB[BUSCO]", self.get_data_path("busco_db")
        )

        fake_partition = MagicMock(
            values=MagicMock(return_value=["partition1", "partition2"])
        )

        def fake_filter_contigs(*args, **kwargs):
            return ("filtered_unbinned",)

        mock_action = MagicMock(
            side_effect=[
                lambda x, y, z, **kwargs: (0,),  # evaluate_busco
                lambda x: ("collated_result",),  # collate
                lambda x: ("visualization",),  # visualize
                fake_filter_contigs,  # filter unbinned
                lambda x, y: (fake_partition,),  # partition
            ]
        )
        mock_ctx = MagicMock(get_action=mock_action)
        obs = evaluate_busco(
            ctx=mock_ctx,
            mags=mags,
            unbinned_contigs=None,
            db=busco_db,
            num_partitions=2,
        )
        exp = ("collated_result", "visualization")
        self.assertTupleEqual(obs, exp)

    @patch("q2_annotate.busco.busco._validate_parameters")
    def test_evaluate_busco_action_with_unbinned(self, mock_validate):
        mags = qiime2.Artifact.import_data(
            "FeatureData[MAG]", self.get_data_path("mags/sample2")
        )
        unbinned = qiime2.Artifact.import_data(
            "SampleData[Contigs]", self.get_data_path("unbinned")
        )
        busco_db = qiime2.Artifact.import_data(
            "ReferenceDB[BUSCO]", self.get_data_path("busco_db")
        )

        fake_partition = MagicMock(
            values=MagicMock(return_value=["partition1", "partition2"])
        )

        def fake_filter_contigs(*args, **kwargs):
            return ("filtered_unbinned",)

        mock_action = MagicMock(
            side_effect=[
                lambda x, y, z, **kwargs: (0,),  # evaluate_busco
                lambda x: ("collated_result",),  # collate
                lambda x: ("visualization",),  # visualize
                fake_filter_contigs,  # filter unbinned
                lambda x, y: (fake_partition,),  # partition
            ]
        )
        mock_ctx = MagicMock(get_action=mock_action)

        with pytest.warns(match="unbinned contigs will be ignored"):
            obs = evaluate_busco(
                ctx=mock_ctx,
                mags=mags,
                unbinned_contigs=unbinned,
                db=busco_db,
                num_partitions=2,
            )
        exp = ("collated_result", "visualization")
        self.assertTupleEqual(obs, exp)
