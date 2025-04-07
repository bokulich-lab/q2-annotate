# ----------------------------------------------------------------------------
# Copyright (c) 2022-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from pathlib import Path

import pandas as pd
from pandas.testing import assert_frame_equal

from qiime2.plugin.testing import TestPluginBase
from q2_types.kraken2 import (
    Kraken2ReportFormat,
    Kraken2OutputFormat,
    Kraken2ReportDirectoryFormat,
    Kraken2OutputDirectoryFormat,
)

from q2_annotate.kraken2.filter import _report_df_to_tree, _dump_tree_to_report
from q2_annotate.kraken2.merge import _merge_trees, _merge_kraken2_results


class TestTreeMerging(TestPluginBase):
    package = "q2_annotate.kraken2.tests"

    def _get_tree_df(self, fp):
        fp = self.get_data_path(Path('merge/tree-merging') / fp)
        return Kraken2ReportFormat(fp, mode='r').view(pd.DataFrame)

    def _get_tree(self, fp):
        return _report_df_to_tree(self._get_tree_df(fp))

    def test_merge_trees_no_unclassified_nodes(self):
        '''
        Tests that two reports are merged as expected. Tree 1 and tree 2 each
        have some unique nodes and some shared nodes, so tests that node
        grafting is performed properly and tests that overlapping nodes are
        combined properly.

        Note that the ordering of taxa within the expected and observed reports
        is different, hence the `sort_values` and `reindex`. Despite different
        taxon ordering, both reports represent the same information.
        '''
        tree1 = self._get_tree(
            'no-unclassified/tree-1.report.txt'
        )
        tree2 = self._get_tree(
            'no-unclassified/tree-2.report.txt'
        )
        obs_df = _dump_tree_to_report(*_merge_trees(tree1, tree2))

        exp_df = self._get_tree_df(
            'no-unclassified/merged-tree.report.txt'
        )

        assert_frame_equal(
            obs_df.sort_values(by='perc_frags_covered').reset_index(drop=True),
            exp_df.sort_values(by='perc_frags_covered').reset_index(drop=True),
            check_dtype=False
        )

    def test_merge_trees_one_unclassified_nodes(self):
        '''
        Same as `test_merge_trees_no_unclassified_nodes` except one report
        has an unclassified node.
        '''
        tree1 = self._get_tree(
            'one-unclassified/tree-1.report.txt'
        )
        tree2 = self._get_tree(
            'one-unclassified/tree-2.report.txt'
        )
        obs_df = _dump_tree_to_report(*_merge_trees(tree1, tree2))

        exp_df = self._get_tree_df(
            'one-unclassified/merged-tree.report.txt'
        )

        assert_frame_equal(
            obs_df.sort_values(by='perc_frags_covered').reset_index(drop=True),
            exp_df.sort_values(by='perc_frags_covered').reset_index(drop=True),
            check_dtype=False
        )

    def test_merge_trees_two_unclassified_nodes(self):
        '''
        Same as `test_merge_trees_no_unclassified_nodes` except both reports
        have an unclassified node.
        '''
        tree1 = self._get_tree(
            'two-unclassified/tree-1.report.txt'
        )
        tree2 = self._get_tree(
            'two-unclassified/tree-2.report.txt'
        )
        obs_df = _dump_tree_to_report(*_merge_trees(tree1, tree2))

        exp_df = self._get_tree_df(
            'two-unclassified/merged-tree.report.txt'
        )

        assert_frame_equal(
            obs_df.sort_values(by='perc_frags_covered').reset_index(drop=True),
            exp_df.sort_values(by='perc_frags_covered').reset_index(drop=True),
            check_dtype=False
        )


class TestResultMerging(TestPluginBase):
    package = "q2_annotate.kraken2.tests"

    def setUp(self):
        super().setUp()

        reports_fp = self.get_data_path(
            Path('merge') / 'result-merging' / 'reports'
        )
        outputs_fp = self.get_data_path(
            Path('merge') / 'result-merging' / 'outputs'
        )
        self.first_reports = Kraken2ReportDirectoryFormat(
            Path(reports_fp) / 'artifact-1', mode='r'
        )
        self.second_reports = Kraken2ReportDirectoryFormat(
            Path(reports_fp) / 'artifact-2', mode='r'
        )

        self.first_outputs = Kraken2OutputDirectoryFormat(
            Path(outputs_fp) / 'artifact-1', mode='r'
        )
        self.second_outputs = Kraken2OutputDirectoryFormat(
            Path(outputs_fp) / 'artifact-2', mode='r'
        )

        expected_fp = self.get_data_path(
            Path('merge') / 'result-merging' / 'expected'
        )
        self.expected_reports = Kraken2ReportDirectoryFormat(
            Path(expected_fp) / 'reports', mode='r'
        )
        self.expected_outputs = Kraken2OutputDirectoryFormat(
            Path(expected_fp) / 'outputs', mode='r'
        )

    def _assert_reports_equal(
        self, first: Kraken2ReportFormat, second: Kraken2ReportFormat
    ):
        first_df = first.view(
            pd.DataFrame
        ).sort_values(by='perc_frags_covered').reset_index(drop=True)
        second_df = second.view(
            pd.DataFrame
        ).sort_values(by='perc_frags_covered').reset_index(drop=True)

        assert_frame_equal(first_df, second_df)

    def _assert_outputs_equal(
        self, first: Kraken2OutputFormat, second: Kraken2OutputFormat
    ):
        first_df = first.view(
            pd.DataFrame
        ).sort_values(by='sequence_id').reset_index(drop=True)
        second_df = second.view(
            pd.DataFrame
        ).sort_values(by='sequence_id').reset_index(drop=True)

        assert_frame_equal(first_df, second_df)

    def _assert_formats_equal(self, first_fp: str, second_fp: str, type: str):
        if type == 'reports':
            first = Kraken2ReportFormat(first_fp, mode='r')
            second = Kraken2ReportFormat(second_fp, mode='r')
            self._assert_reports_equal(first, second)
        elif type == 'outputs':
            first = Kraken2OutputFormat(first_fp, mode='r')
            second = Kraken2OutputFormat(second_fp, mode='r')
            self._assert_outputs_equal(first, second)

    def _assert_directory_formats_equal(
        self, first, second, type: str, mags=False
    ):
        first_file_dict = first.file_dict()
        second_file_dict = second.file_dict()

        self.assertEqual(
            list(first_file_dict.keys()), list(second_file_dict.keys())
        )

        if mags:
            for sample_id in first_file_dict:
                self.assertEqual(
                    list(first_file_dict[sample_id].keys()),
                    list(second_file_dict[sample_id].keys())
                )

        if mags:
            for sample_id in first_file_dict:
                for filename in first_file_dict[sample_id]:
                    first_fp = first_file_dict[sample_id][filename]
                    second_fp = second_file_dict[sample_id][filename]
                    self._assert_formats_equal(first_fp, second_fp, type)
        else:
            for sample_id in first_file_dict:
                first_fp = first_file_dict[sample_id]
                second_fp = second_file_dict[sample_id]
                self._assert_formats_equal(first_fp, second_fp, type)

    def test_result_merging(self):
        obs_reports, obs_outputs = _merge_kraken2_results(
            [self.first_reports, self.second_reports],
            [self.first_outputs, self.second_outputs]
        )

        self._assert_directory_formats_equal(
            obs_reports, self.expected_reports, type='reports'
        )
        self._assert_directory_formats_equal(
            obs_outputs, self.expected_outputs, type='outputs'
        )
