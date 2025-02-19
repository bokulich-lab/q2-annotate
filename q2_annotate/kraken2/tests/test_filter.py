# ----------------------------------------------------------------------------
# Copyright (c) 2022-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pandas as pd
from pandas.testing import assert_series_equal

from qiime2.plugin.testing import TestPluginBase
from q2_types.kraken2 import Kraken2ReportFormat, Kraken2ReportDirectoryFormat

from q2_annotate.kraken2.filter import (
    _report_df_to_tree,
    _trim_tree_dfs,
    _dump_tree_to_report,
    _write_report_dfs,
    filter_kraken2_reports,
)


class TestFilter(TestPluginBase):
    package = "q2_annotate.kraken2.tests"

    def setUp(self):
        super().setUp()

        curated_report_fp = self.get_data_path('filter/sample-1.report.txt')
        curated_report = Kraken2ReportFormat(curated_report_fp, mode='r')
        self.curated_report = curated_report.view(pd.DataFrame)

        real_report_fp = self.get_data_path('filter/SRR17001003.report.txt')
        real_report = Kraken2ReportFormat(real_report_fp, mode='r')
        self.real_report = real_report.view(pd.DataFrame)

        reports_fp = self.get_data_path('filter')
        self.reports = Kraken2ReportDirectoryFormat(reports_fp, mode='r')

    def test_curated_report_to_tree(self):
        '''
        Tests that a curated kraken2 report with no unclassified reads is
        properly parsed into a tree.
        '''
        root, unclassified_node = _report_df_to_tree(self.curated_report)

        # assert unclassified_node looks correct
        self.assertEqual(unclassified_node, None)

        # assert number of nodes correct
        self.assertEqual(root.count(), 20)

        # assert root is root
        self.assertTrue(root.is_root())
        self.assertEqual(root._kraken_data['name'], 'root')
        self.assertEqual(root._kraken_data['rank'], 'R')

        # assert num leaves is correct
        self.assertEqual(len(list(root.tips())), 9)

        # find node 1767
        nodes_list = list(root.find_by_func(
            lambda n: n._kraken_data['taxon_id'] == 1767
        ))

        node_1767 = nodes_list.pop(0)

        # only one node has taxon id 1767
        self.assertEqual(len(nodes_list), 0)

        # siblings match
        self.assertEqual(len(node_1767.siblings()), 5)
        node_1767_sibling_ids = [
            n._kraken_data['taxon_id'] for n in node_1767.siblings()
        ]
        self.assertEqual(
            set(node_1767_sibling_ids),
            set([1138383, 701042, 1764, 339268, 2775496])
        )

        # correct parent
        self.assertEqual(node_1767.parent._kraken_data['taxon_id'], 120793)

        # correct depth
        self.assertEqual(len(node_1767.ancestors()), 10)

        # assert node has proper kraken data
        exp = {
            'perc_frags_covered': 48.55,
            'n_frags_covered': 230249,
            'n_frags_assigned': 230249,
            'rank': 'S',
            'taxon_id': 1767,
            'name': 'Mycobacterium intracellulare'
        }
        obs = node_1767._kraken_data
        obs['name'] = obs['name'].strip()
        self.assertEqual(obs, exp)

    def test_real_report_to_tree(self):
        '''
        Tests that a real kraken2 report with unclassified reads is
        properly parsed into a tree and unclassified node.
        '''
        root, unclassified_node = _report_df_to_tree(self.real_report)

        # assert unclassified_node looks correct
        exp = {
            'perc_frags_covered': 4.09,
            'n_frags_covered': 332879,
            'n_frags_assigned': 332879,
            'n_read_minimizers': 0,
            'n_uniq_minimizers': 0,
            'rank': 'U',
            'taxon_id': 0,
            'name': 'unclassified'
        }
        obs = unclassified_node._kraken_data

        self.assertEqual(obs, exp)

        # assert number of nodes correct
        self.assertEqual(root.count(), 21237)

        # assert root is root
        self.assertTrue(root.is_root())
        self.assertEqual(root._kraken_data['name'], 'root')
        self.assertEqual(root._kraken_data['rank'], 'R')

        # assert num leaves less than num nodes
        self.assertLess(len(list(root.tips())), root.count())

        # find node 2732008
        nodes_list = list(root.find_by_func(
            lambda n: n._kraken_data['taxon_id'] == 2732008
        ))

        node_2732008 = nodes_list.pop(0)

        # only one node has taxon id 2732008
        self.assertEqual(len(nodes_list), 0)

        # correct parent
        self.assertEqual(node_2732008.parent._kraken_data['taxon_id'], 2732005)

        # correct children
        self.assertEqual(len(node_2732008.children), 2)
        node_2732008_children_ids = [
            n._kraken_data['taxon_id'] for n in node_2732008.children
        ]
        self.assertEqual(
            set(node_2732008_children_ids), set([2732529, 3044425])
        )

        # correct depth
        self.assertEqual(len(node_2732008.ancestors()), 4)

        # assert node has proper kraken data
        exp = {
            'perc_frags_covered': 0.00,
            'n_frags_covered': 3,
            'n_frags_assigned': 0,
            'n_read_minimizers': 28,
            'n_uniq_minimizers': 16,
            'rank': 'P',
            'taxon_id': 2732008,
            'name': 'Preplasmiviricota'
        }
        obs = node_2732008._kraken_data
        obs['name'] = obs['name'].strip()
        self.assertEqual(obs, exp)

    def test_trim_curated_tree(self):
        '''
        Tests that a curated report tree is properly trimmed.
        '''
        def find_by_name(name, root):
            nodes_list = list(root.find_by_func(
                lambda n: n._kraken_data['name'].strip() == name
            ))
            return nodes_list[0] if len(nodes_list) else None

        root, unclassified_node = _report_df_to_tree(self.curated_report)

        total_reads = root._kraken_data['n_frags_covered']

        filtered_root = _trim_tree_dfs(
            root.copy(deep=True),
            abundance_threshold=0.003,
            total_reads=total_reads
        )

        # number of nodes left is correct
        self.assertEqual(filtered_root.count(), 13)

        # some trimmed nodes are gone
        self.assertIsNotNone(find_by_name('Mycobacterium marseillense', root))
        self.assertIsNone(
            find_by_name('Mycobacterium marseillense', filtered_root)
        )

        self.assertIsNotNone(find_by_name('Mycobacterium gordonae', root))
        self.assertIsNone(
            find_by_name('Mycobacterium gordonae', filtered_root)
        )

        # some trimmed node's children should be gone
        self.assertIsNotNone(find_by_name('Mycobacterium tuberculosis', root))
        self.assertIsNone(
            find_by_name('Mycobacterium tuberculosis', filtered_root)
        )

        # a node with no n_frags_assigned whose children were all removed
        # is removed
        self.assertIsNotNone(
            find_by_name('Mycobacterium tuberculosis complex', root)
        )
        self.assertIsNone(
            find_by_name('Mycobacterium tuberculosis complex', filtered_root)
        )

        # trimmed node's ancestors have n_frags_covered updated
        myco_avium = find_by_name('Mycobacterium avium complex (MAC)', root)
        myco_avium_filtered = find_by_name(
            'Mycobacterium avium complex (MAC)', filtered_root
        )

        filtered_reads = 180

        self.assertEqual(
            myco_avium_filtered._kraken_data['n_frags_covered'],
            myco_avium._kraken_data['n_frags_covered'] - filtered_reads,
        )

        ancestors = list(myco_avium.ancestors())
        filtered_ancestors = list(myco_avium_filtered.ancestors())

        filtered_reads = 410

        for ancestor_index, _ in enumerate(ancestors):
            filtered_ancestor = filtered_ancestors[ancestor_index]
            ancestor = ancestors[ancestor_index]
            self.assertEqual(
                filtered_ancestor._kraken_data['n_frags_covered'],
                ancestor._kraken_data['n_frags_covered'] - filtered_reads,
            )

    def test_dump_tree_to_report_round_trip(self):
        '''
        Test that a report tree is properly dumped to a report dataframe.
        '''
        root, unclassified_node = _report_df_to_tree(self.curated_report)

        round_tripped_report = _dump_tree_to_report(root, unclassified_node)

        round_tripped_report.sort_values(
            'taxon_id', inplace=True, ascending=False
        )
        curated_report = self.curated_report.sort_values(
            'taxon_id', ascending=False
        )

        columns = [
            'n_frags_covered', 'n_frags_assigned', 'taxon_id', 'name', 'rank'
        ]
        for column in columns:
            assert_series_equal(
                round_tripped_report[column],
                curated_report[column],
                check_dtype=False,
                check_index=False,
            )

    def test_filter_kraken2_reports(self):
        '''
        '''
        filtered_reports = filter_kraken2_reports(
            self.reports, abundance_threshold=0.01
        )
