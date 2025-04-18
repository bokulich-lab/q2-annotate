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

    def test_result_merging_reads(self):
        '''
        Tests that two kraken2 reports and two kraken2 outputs input artifacts
        all with a single and identical sample ID get properly merged into
        a single kraken2 reports output and a single kraken2 outputs output
        (phew).
        '''
        reports_dir = self.get_data_path(
            Path('merge') / 'result-merging' / 'reads' / 'reports'
        )
        outputs_dir = self.get_data_path(
            Path('merge') / 'result-merging' / 'reads' / 'outputs'
        )
        first_reports = Kraken2ReportDirectoryFormat(
            Path(reports_dir) / 'artifact-1', mode='r'
        )
        second_reports = Kraken2ReportDirectoryFormat(
            Path(reports_dir) / 'artifact-2', mode='r'
        )

        first_outputs = Kraken2OutputDirectoryFormat(
            Path(outputs_dir) / 'artifact-1', mode='r'
        )
        second_outputs = Kraken2OutputDirectoryFormat(
            Path(outputs_dir) / 'artifact-2', mode='r'
        )

        expected_dir = self.get_data_path(
            Path('merge') / 'result-merging' / 'reads' / 'expected'
        )
        exp_reports = Kraken2ReportDirectoryFormat(
            Path(expected_dir) / 'reports', mode='r'
        )
        exp_outputs = Kraken2OutputDirectoryFormat(
            Path(expected_dir) / 'outputs', mode='r'
        )

        obs_reports, obs_outputs = _merge_kraken2_results(
            [first_reports, second_reports], [first_outputs, second_outputs]
        )

        self._assert_directory_formats_equal(
            obs_reports, exp_reports, type='reports'
        )
        self._assert_directory_formats_equal(
            obs_outputs, exp_outputs, type='outputs'
        )

    def test_result_merging_mags(self):
        '''
        Tests that mag report/output directories for the same sample ID have
        their constituent report files merged into a single directory in the
        output data.
        '''
        reports_dir = self.get_data_path(
            Path('merge') / 'result-merging' / 'mags' / 'reports'
        )
        outputs_dir = self.get_data_path(
            Path('merge') / 'result-merging' / 'mags' / 'outputs'
        )
        first_reports = Kraken2ReportDirectoryFormat(
            Path(reports_dir) / 'artifact-1', mode='r'
        )
        second_reports = Kraken2ReportDirectoryFormat(
            Path(reports_dir) / 'artifact-2', mode='r'
        )
        first_outputs = Kraken2OutputDirectoryFormat(
            Path(outputs_dir) / 'artifact-1', mode='r'
        )
        second_outputs = Kraken2OutputDirectoryFormat(
            Path(outputs_dir) / 'artifact-2', mode='r'
        )

        expected_dir = self.get_data_path(
            Path('merge') / 'result-merging' / 'mags' / 'expected'
        )
        exp_reports = Kraken2ReportDirectoryFormat(
            Path(expected_dir) / 'reports', mode='r'
        )
        exp_outputs = Kraken2OutputDirectoryFormat(
            Path(expected_dir) / 'outputs', mode='r'
        )

        obs_reports, obs_outputs = _merge_kraken2_results(
            [first_reports, second_reports], [first_outputs, second_outputs]
        )

        self._assert_directory_formats_equal(
            obs_reports, exp_reports, type='reports', mags=True
        )
        self._assert_directory_formats_equal(
            obs_outputs, exp_outputs, type='outputs', mags=True
        )

    def test_duplicate_mag_ids_cause_error(self):
        '''
        Tests that if two directories with the same sample ID within a report
        directory format or an output directory format have a report with the
        same hash (representing the same MAG) then an error is raised.
        '''
        reports_dir = self.get_data_path(
            Path('merge') / 'result-merging' / 'mags' / 'reports'
        )
        outputs_dir = self.get_data_path(
            Path('merge') / 'result-merging' / 'mags' / 'outputs'
        )
        first_reports = Kraken2ReportDirectoryFormat(
            Path(reports_dir) / 'artifact-1', mode='r'
        )
        second_reports = Kraken2ReportDirectoryFormat(
            Path(reports_dir) / 'artifact-1', mode='r'
        )
        first_outputs = Kraken2OutputDirectoryFormat(
            Path(outputs_dir) / 'artifact-1', mode='r'
        )
        second_outputs = Kraken2OutputDirectoryFormat(
            Path(outputs_dir) / 'artifact-1', mode='r'
        )

        with self.assertRaisesRegex(ValueError, 'Two MAGs with the same uuid'):
            _merge_kraken2_results(
                [first_reports, second_reports], [first_outputs, second_outputs]
            )

    def test_merging_with_minimizers_and_sample_id_overlap_errors(self):
        '''
        Tests that when reports contain minimizer information and two or more
        reports with the same sample ID from different input artifacts are to
        be merged, an error is raised.
        '''
        reports_dir = self.get_data_path(
            Path('merge') / 'result-merging' / 'minimizers' / 'reports'
        )
        outputs_dir = self.get_data_path(
            Path('merge') / 'result-merging' / 'minimizers' / 'outputs'
        )
        first_reports = Kraken2ReportDirectoryFormat(
            Path(reports_dir) / 'artifact-1', mode='r'
        )
        second_reports = Kraken2ReportDirectoryFormat(
            Path(reports_dir) / 'artifact-2', mode='r'
        )
        first_outputs = Kraken2OutputDirectoryFormat(
            Path(outputs_dir) / 'artifact-1', mode='r'
        )
        second_outputs = Kraken2OutputDirectoryFormat(
            Path(outputs_dir) / 'artifact-2', mode='r'
        )

        with self.assertRaisesRegex(ValueError, 'option to capture minimzer'):
            _merge_kraken2_results(
                [first_reports, second_reports], [first_outputs, second_outputs]
            )

    def test_merging_with_minimzer_and_no_sample_id_overlap_succeeds(self):
        '''
        Tests that reports with minimizer information present across multiple
        input artifacts can still be successfully condensed into a single
        set of output artifacts, as long as there is no direct merging of
        reports on a per-sample ID basis.
        '''
        reports_dir = self.get_data_path(
            Path('merge') / 'result-merging' / 'minimizers' / 'reports'
        )
        outputs_dir = self.get_data_path(
            Path('merge') / 'result-merging' / 'minimizers' / 'outputs'
        )
        first_reports = Kraken2ReportDirectoryFormat(
            Path(reports_dir) / 'artifact-1', mode='r'
        )
        second_reports = Kraken2ReportDirectoryFormat(
            Path(reports_dir) / 'artifact-3', mode='r'
        )
        first_outputs = Kraken2OutputDirectoryFormat(
            Path(outputs_dir) / 'artifact-1', mode='r'
        )
        second_outputs = Kraken2OutputDirectoryFormat(
            Path(outputs_dir) / 'artifact-3', mode='r'
        )

        expected_dir = self.get_data_path(
            Path('merge') / 'result-merging' / 'minimizers' / 'expected'
        )
        exp_reports = Kraken2ReportDirectoryFormat(
            Path(expected_dir) / 'reports', mode='r'
        )
        exp_outputs = Kraken2OutputDirectoryFormat(
            Path(expected_dir) / 'outputs', mode='r'
        )

        obs_reports, obs_outputs = _merge_kraken2_results(
            [first_reports, second_reports], [first_outputs, second_outputs]
        )

        self._assert_directory_formats_equal(
            obs_reports, exp_reports, type='reports'
        )
        self._assert_directory_formats_equal(
            obs_outputs, exp_outputs, type='outputs'
        )
