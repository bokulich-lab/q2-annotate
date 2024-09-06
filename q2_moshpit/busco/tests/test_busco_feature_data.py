# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json
import os
import shutil
import qiime2
import pandas as pd
from q2_moshpit.busco.busco import (
    _run_busco, _visualize_busco, evaluate_busco
)
from unittest.mock import patch, ANY, MagicMock
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data_mag import MAGSequencesDirFmt


class TestBUSCOFeatureData(TestPluginBase):
    package = "q2_moshpit.busco.tests"

    def setUp(self):
        super().setUp()
        self.mags = MAGSequencesDirFmt(
            path=self.get_data_path('mags/sample1'),
            mode="r",
        )

    def _prepare_summaries(self):
        os.makedirs(os.path.join(self.temp_dir.name, "sample1"))
        shutil.copy(
            self.get_data_path('summaries/batch_summary_1.txt'),
            os.path.join(
                self.temp_dir.name, "sample1", 'batch_summary.txt'
            )
        )

    @patch('q2_moshpit.busco.busco.run_command')
    def test_run_busco(self, mock_run):
        self._prepare_summaries()

        obs = _run_busco(
            output_dir=self.temp_dir.name,
            mags=self.mags,
            params=['--lineage_dataset', 'bacteria_odb10', '--cpu', '7']
        )
        exp = {
            'sample1': f"{self.temp_dir.name}/sample1/batch_summary.txt",
        }

        self.assertDictEqual(obs, exp)
        mock_run.assert_called_once_with([
            'busco', '--lineage_dataset', 'bacteria_odb10',
            '--cpu', '7', '--in', self.get_data_path('mags/sample1'),
            '--out_path', self.temp_dir.name, '-o', 'sample1'
        ], cwd=os.path.dirname(self.temp_dir.name))

    @patch(
        "q2_moshpit.busco.busco._draw_detailed_plots",
        return_value={"fake1": {"plot": "spec"}}
    )
    @patch(
        "q2_moshpit.busco.busco._draw_marker_summary_histograms",
        return_value={"fake2": {"plot": "spec"}}
    )
    @patch(
        "q2_moshpit.busco.busco._draw_selectable_summary_histograms",
        return_value={"fake3": {"plot": "spec"}}
    )
    @patch(
        "q2_moshpit.busco.busco._get_feature_table", return_value="table1"
    )
    @patch(
        "q2_moshpit.busco.busco._calculate_summary_stats",
        return_value="stats1"
    )
    @patch("q2templates.render")
    @patch("q2_moshpit.busco.busco._cleanup_bootstrap")
    def test_visualize_busco(
            self, mock_clean, mock_render, mock_stats, mock_table,
            mock_selectable, mock_marker, mock_detailed
    ):
        _visualize_busco(
            output_dir=self.temp_dir.name,
            busco_results=pd.read_csv(
                self.get_data_path(
                    'summaries/all_renamed_with_lengths_feature_data.csv'
                )
            )
        )

        mock_detailed.assert_called_once()
        mock_marker.assert_called_once()

        exp_context = {
            "tabs": [
                {"title": "QC overview", "url": "index.html"},
                {"title": "BUSCO plots", "url": "detailed_view.html"},
                {"title": "BUSCO table", "url": "table.html"}
            ],
            "vega_json": json.dumps(
                {
                    "partition_0": {
                        "subcontext": {"fake1": {"plot": "spec"}},
                        "counters": {"from": 1, "to": 2},
                        "ids": [
                            "ab23d75d-547d-455a-8b51-16b46ddf7496",
                            "0e514d88-16c4-4273-a1df-1a360eb2c823"
                        ]
                    }
                }
            ),
            "vega_summary_json": json.dumps({"fake2": {"plot": "spec"}}),
            "table": "table1",
            "summary_stats_json": "stats1",
            "page_size": 100
        }
        mock_render.assert_called_with(
            ANY, self.temp_dir.name, context=exp_context
        )
        mock_clean.assert_called_with(self.temp_dir.name)

    # TODO: maybe this could be turned into an actual test
    def test_evaluate_busco_action(self):
        mock_action = MagicMock(side_effect=[
            lambda x, **kwargs: (0, ),
            lambda x: ("collated_result", ),
            lambda x: ("visualization", ),
            lambda x, y: ({"mag1": {}, "mag2": {}}, )
        ])
        mock_ctx = MagicMock(get_action=mock_action)
        mags = qiime2.Artifact.import_data(
            'FeatureData[MAG]',
            self.get_data_path('mags/sample2')
        )
        busco_db = qiime2.Artifact.import_data(
            'ReferenceDB[BuscoDB]',
            self.get_data_path('busco_db')
        )
        obs = evaluate_busco(
            ctx=mock_ctx,
            bins=mags,
            busco_db=busco_db,
            num_partitions=2
        )
        exp = ("collated_result", "visualization")
        self.assertTupleEqual(obs, exp)
