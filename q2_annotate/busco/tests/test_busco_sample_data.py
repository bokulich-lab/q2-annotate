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
from q2_annotate.busco.busco import (
    _run_busco, _busco_helper, _evaluate_busco,
    _visualize_busco, evaluate_busco
)
from unittest.mock import patch, ANY, call, MagicMock
from qiime2.plugin.testing import TestPluginBase
from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt
from q2_annotate.busco.types import BuscoDatabaseDirFmt


class TestBUSCOSampleData(TestPluginBase):
    package = "q2_annotate.busco.tests"

    def setUp(self):
        super().setUp()
        self.mags = MultiMAGSequencesDirFmt(
            path=self.get_data_path('mags'),
            mode="r",
        )
        self.busco_db = BuscoDatabaseDirFmt(
            path=self.get_data_path("busco_db"),
            mode="r"
        )

    def _prepare_summaries(self):
        for s in ['1', '2']:
            os.makedirs(os.path.join(self.temp_dir.name, f"sample{s}"))
            shutil.copy(
                self.get_data_path(f'summaries/batch_summary_{s}.txt'),
                os.path.join(
                    self.temp_dir.name, f"sample{s}", 'batch_summary.txt'
                )
            )

    @patch('q2_annotate.busco.busco.run_command')
    def test_run_busco(self, mock_run):
        self._prepare_summaries()

        obs = _run_busco(
            output_dir=self.temp_dir.name,
            mags=self.mags,
            params=['--lineage_dataset', 'bacteria_odb10', '--cpu', '7']
        )
        exp = {
            'sample1': f"{self.temp_dir.name}/sample1/batch_summary.txt",
            'sample2': f"{self.temp_dir.name}/sample2/batch_summary.txt",
        }

        self.assertDictEqual(obs, exp)
        mock_run.assert_has_calls([
            call(
                [
                    'busco', '--lineage_dataset', 'bacteria_odb10',
                    '--cpu', '7', '--in', self.get_data_path('mags/sample1'),
                    '--out_path', self.temp_dir.name, '-o', 'sample1'
                ],
                cwd=os.path.dirname(self.temp_dir.name)
            ),
            call(
                [
                    'busco', '--lineage_dataset', 'bacteria_odb10',
                    '--cpu', '7', '--in', self.get_data_path('mags/sample2'),
                    '--out_path', self.temp_dir.name, '-o', 'sample2'
                ],
                cwd=os.path.dirname(self.temp_dir.name)
            ),
        ])

    @patch('q2_annotate.busco.busco._run_busco')
    @patch('q2_annotate.busco.busco._get_mag_lengths')
    def test_busco_helper(self, mock_len, mock_run):
        self._prepare_summaries()
        mock_run.return_value = {
            'sample1': f"{self.temp_dir.name}/sample1/batch_summary.txt",
            'sample2': f"{self.temp_dir.name}/sample2/batch_summary.txt",
        }
        mock_len.return_value = pd.Series(
            {
                'ab23d75d-547d-455a-8b51-16b46ddf7496': 3177889,
                '0e514d88-16c4-4273-a1df-1a360eb2c823': 19625518,
                '8098e3a8-df4a-46af-83e2-6c2443d74cb9': 912062,
                'd4637408-80ab-49d5-ab32-c66509c3a544': 2678999,
                '503c2f56-3e4f-4ce7-9b61-b63bc7fe0592': 21035714
            }, name='length'
        )

        obs = _busco_helper(
            self.mags, ['--lineage_dataset', 'bacteria_odb10']
        )
        exp = pd.read_csv(
            self.get_data_path('summaries/all_renamed_with_lengths.csv'),
        )

        pd.testing.assert_frame_equal(obs, exp)
        mock_run.assert_called_with(
            output_dir=ANY, mags=self.mags,
            params=['--lineage_dataset', 'bacteria_odb10']
        )

    @patch("q2_annotate.busco.busco._busco_helper")
    def test_evaluate_busco_offline(self, mock_helper):
        _evaluate_busco(
            mags=self.mags,
            db=self.busco_db,
            mode="some_mode",
            lineage_dataset="lineage_1"
        )
        mock_helper.assert_called_with(
            self.mags,
            [
                '--mode', 'some_mode', '--lineage_dataset', 'lineage_1',
                '--cpu', '1', '--contig_break', '10', '--evalue', '0.001',
                '--limit', '3', '--offline', "--download_path",
                f"{str(self.busco_db)}/busco_downloads"
            ]
        )

    @patch(
        "q2_annotate.busco.busco._draw_detailed_plots",
        return_value={"fake1": {"plot": "spec"}}
    )
    @patch(
        "q2_annotate.busco.busco._draw_marker_summary_histograms",
        return_value={"fake2": {"plot": "spec"}}
    )
    @patch(
        "q2_annotate.busco.busco._draw_selectable_summary_histograms",
        return_value={"fake3": {"plot": "spec"}}
    )
    @patch(
        "q2_annotate.busco.busco._get_feature_table", return_value="table1"
    )
    @patch(
        "q2_annotate.busco.busco._calculate_summary_stats",
        return_value="stats1"
    )
    @patch("q2templates.render")
    @patch("q2_annotate.busco.busco._cleanup_bootstrap")
    def test_visualize_busco(
            self, mock_clean, mock_render, mock_stats, mock_table,
            mock_selectable, mock_marker, mock_detailed
    ):
        _visualize_busco(
            output_dir=self.temp_dir.name,
            results=pd.read_csv(
                self.get_data_path('summaries/all_renamed_with_lengths.csv')
            )
        )

        mock_detailed.assert_called_once()
        mock_marker.assert_called_once()
        mock_selectable.assert_called_once()

        exp_context = {
            "tabs": [
                {"title": "QC overview", "url": "index.html"},
                {"title": "Sample details", "url": "detailed_view.html"},
                {"title": "Feature details", "url": "table.html"}
            ],
            "vega_json": json.dumps(
                {"partition_0": {
                    "subcontext": {"fake1": {"plot": "spec"}},
                    "counters": {"from": 1, "to": 2},
                    "ids": ["sample1", "sample2"]}}
            ),
            "vega_summary_json": json.dumps({"fake2": {"plot": "spec"}}),
            "vega_summary_selectable_json": json.dumps(
                {"fake3": {"plot": "spec"}}
            ),
            "table": "table1",
            "summary_stats_json": "stats1",
            "page_size": 100
        }
        mock_render.assert_called_with(
            ANY, self.temp_dir.name, context=exp_context
        )
        mock_clean.assert_called_with(self.temp_dir.name)

    # TODO: maybe this could be turned into an actual test
    @patch('q2_annotate.busco.busco._validate_parameters')
    def test_evaluate_busco_action(self, mock_validate):
        mock_action = MagicMock(side_effect=[
            lambda x, **kwargs: (0, ),
            lambda x: ("collated_result", ),
            lambda x: ("visualization", ),
            lambda x, y: ({"mag1": {}, "mag2": {}}, )
        ])
        mock_ctx = MagicMock(get_action=mock_action)
        mags = qiime2.Artifact.import_data(
            'SampleData[MAGs]',
            self.get_data_path('mags')
        )
        busco_db = qiime2.Artifact.import_data(
            'ReferenceDB[BuscoDB]',
            self.get_data_path('busco_db')
        )
        obs = evaluate_busco(
            ctx=mock_ctx,
            mags=mags,
            db=busco_db,
            num_partitions=2
        )
        exp = ("collated_result", "visualization")
        self.assertTupleEqual(obs, exp)
