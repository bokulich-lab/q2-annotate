# ----------------------------------------------------------------------------
# Copyright (c) 2022-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pandas as pd
from pandas.testing import assert_frame_equal

from qiime2.plugin.testing import TestPluginBase
from q2_types.kraken2 import (
    Kraken2ReportFormat,
    Kraken2OutputFormat,
    Kraken2ReportDirectoryFormat,
)

from q2_annotate.kraken2.merge import _merge_trees
from q2_annotate.kraken2.filter import _report_df_to_tree, _dump_tree_to_report

class TestTreeMerging(TestPluginBase):
    package = "q2_annotate.kraken2.tests"

    def _get_tree_df(self, fp):
        fp = self.get_data_path(fp)
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
            'report-merging/no-unclassified/tree-1.report.txt'
        )
        tree2 = self._get_tree(
            'report-merging/no-unclassified/tree-2.report.txt'
        )
        obs_df = _dump_tree_to_report(*_merge_trees(tree1, tree2))

        exp_df = self._get_tree_df(
            'report-merging/no-unclassified/merged-tree.report.txt'
        )

        assert_frame_equal(
            obs_df.sort_values(by='perc_frags_covered').reset_index(drop=True),
            exp_df.sort_values(by='perc_frags_covered').reset_index(drop=True),
            check_dtype=False,
        )

    def test_merge_trees_one_unclassified_nodes(self):
        '''
        Same as `test_merge_trees_no_unclassified_nodes` except one report
        has an unclassified node.
        '''
        tree1 = self._get_tree(
            'report-merging/one-unclassified/tree-1.report.txt'
        )
        tree2 = self._get_tree(
            'report-merging/one-unclassified/tree-2.report.txt'
        )
        obs_df = _dump_tree_to_report(*_merge_trees(tree1, tree2))

        exp_df = self._get_tree_df(
            'report-merging/one-unclassified/merged-tree.report.txt'
        )

        assert_frame_equal(
            obs_df.sort_values(by='perc_frags_covered').reset_index(drop=True),
            exp_df.sort_values(by='perc_frags_covered').reset_index(drop=True),
            check_dtype=False,
        )

    def test_merge_trees_two_unclassified_nodes(self):
        '''
        Same as `test_merge_trees_no_unclassified_nodes` except both reports
        have an unclassified node.
        '''
        tree1 = self._get_tree(
            'report-merging/two-unclassified/tree-1.report.txt'
        )
        tree2 = self._get_tree(
            'report-merging/two-unclassified/tree-2.report.txt'
        )
        obs_df = _dump_tree_to_report(*_merge_trees(tree1, tree2))

        exp_df = self._get_tree_df(
            'report-merging/two-unclassified/merged-tree.report.txt'
        )

        assert_frame_equal(
            obs_df.sort_values(by='perc_frags_covered').reset_index(drop=True),
            exp_df.sort_values(by='perc_frags_covered').reset_index(drop=True),
            check_dtype=False,
        )
