# ----------------------------------------------------------------------------
# Copyright (c) 2022-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from io import StringIO
import os.path
from pathlib import Path
import shutil
import tempfile
from unittest.mock import patch

import pandas as pd
from pandas.testing import assert_series_equal

import qiime2
from qiime2.plugin.testing import TestPluginBase
from q2_types.kraken2 import (
    Kraken2ReportFormat,
    Kraken2OutputFormat,
    Kraken2ReportDirectoryFormat,
    Kraken2OutputDirectoryFormat,
)

from q2_annotate.kraken2.filter import (
    _validate_parameters,
    _find_empty_reports,
    _create_filtered_results,
    _filter_kraken2_results_by_metadata,
    _validate_ids,
    _report_df_to_tree,
    _trim_tree_dfs,
    _dump_tree_to_report,
    _filter_kraken2_reports_by_abundance,
    _align_outputs_with_reports,
    _align_single_output_with_report,
)


class TestFilterKrakenReports(TestPluginBase):
    package = "q2_annotate.kraken2.tests"

    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        instance = cls()

        cls.report_mags = Kraken2ReportDirectoryFormat(
            instance.get_data_path("reports-mags"), "r"
        )
        cls.output_mags = Kraken2OutputDirectoryFormat(
            instance.get_data_path("outputs-mags"), "r"
        )
        cls.report_mags_unclassified_missing_frac = (
            Kraken2ReportDirectoryFormat(
                instance.get_data_path(
                    "reports-mags-unclassified-missing-frac"),
                "r"))
        cls.report_mags_unclassified_with_root = (
            Kraken2ReportDirectoryFormat(
                instance.get_data_path(
                    "reports-mags-unclassified-root"),
                "r"))

        cls.file_dict_report_mags = cls.report_mags.file_dict()
        cls.file_dict_output_mags = cls.output_mags.file_dict()

        cls.file_dict_report_unclassified = (
            cls.report_mags_unclassified_missing_frac.file_dict())
        cls.file_dict_report_unclassified_root = (
            cls.report_mags_unclassified_with_root.file_dict())

        cls.metadata_df = pd.read_csv(
            instance.get_data_path("metadata/metadata.tsv"),
            sep="\t",
            index_col="id"
        )

        cls.metadata1 = qiime2.Metadata(cls.metadata_df)

        cls.metadata_df.drop(
            "8894435a-c836-4c18-b475-8b38a9ab6c6b", inplace=True)
        cls.metadata2 = qiime2.Metadata(cls.metadata_df)

    def test_find_empty_reports(self):
        empty_reports = _find_empty_reports(
            file_dict={"": self.file_dict_report_mags}
        )
        self.assertEqual(
            empty_reports, {"8894435a-c836-4c18-b475-8b38a9ab6c6b"}
        )

    def test_find_empty_reports_missing_frac(self):
        empty_reports = _find_empty_reports(
            file_dict={"": self.file_dict_report_unclassified}
        )
        self.assertEqual(
            empty_reports, {"8894435a-c836-4c18-b475-8b38a9ab6c6b"}
        )

    def test_find_empty_reports_with_root(self):
        empty_reports = _find_empty_reports(
            file_dict={"": self.file_dict_report_unclassified_root}
        )
        self.assertEqual(
            empty_reports, {"8894435a-c836-4c18-b475-8b38a9ab6c6b"}
        )

    def test_create_filter_results_reports(self):
        result = _create_filtered_results(
            suffix="report",
            file_dict={"sample1": self.file_dict_report_unclassified},
            ids_to_keep={"8894435a-c836-4c18-b475-8b38a9ab6c6b"}
        )
        self.assertTrue(isinstance(result, Kraken2ReportDirectoryFormat))
        self.assertTrue(
            os.path.exists(
                os.path.join(
                    str(result),
                    "sample1",
                    "8894435a-c836-4c18-b475-8b38a9ab6c6b.report.txt"
                )
            )
        )

    def test_create_filter_results_outputs(self):
        result = _create_filtered_results(
            suffix="output",
            file_dict={"sample1": self.file_dict_output_mags},
            ids_to_keep={"8894435a-c836-4c18-b475-8b38a9ab6c6b"}
        )
        self.assertTrue(isinstance(result, Kraken2OutputDirectoryFormat))
        self.assertTrue(
            os.path.exists(
                os.path.join(
                    str(result),
                    "sample1",
                    "8894435a-c836-4c18-b475-8b38a9ab6c6b.output.txt"
                )
            )
        )

    def test_filter_kraken_reports_error(self):
        with self.assertRaisesRegex(
                ValueError, "No IDs remain after filtering."
        ):
            _filter_kraken2_results_by_metadata(
                reports=self.report_mags_unclassified_missing_frac,
                outputs=self.output_mags,
                metadata=self.metadata1,
                exclude_ids=True
            )

    def test_filter_kraken_reports_metadata(self):
        results = _filter_kraken2_results_by_metadata(
            reports=self.report_mags_unclassified_missing_frac,
            outputs=self.output_mags,
            metadata=self.metadata2,
        )
        for i, suffix in zip(range(2), ["report", "output"]):
            self.assertTrue(
                os.path.exists(
                    os.path.join(
                        str(results[i]),
                        f"3b72d1a7-ddb0-4dc7-ac36-080ceda04aaa.{suffix}.txt"
                    )
                )
            )

    def test_filter_kraken_reports_metadata_where(self):
        results = _filter_kraken2_results_by_metadata(
            reports=self.report_mags_unclassified_missing_frac,
            outputs=self.output_mags,
            metadata=self.metadata1,
            where="timepoint>10"
        )
        for i, suffix in zip(range(2), ["report", "output"]):
            self.assertTrue(
                os.path.exists(
                    os.path.join(
                        str(results[i]),
                        f"8894435a-c836-4c18-b475-8b38a9ab6c6b.{suffix}.txt"
                    )
                )
            )

    def test_filter_kraken_reports_where_print(self):
        with patch('sys.stdout', new_callable=StringIO) as mock_stdout:
            _filter_kraken2_results_by_metadata(
                reports=self.report_mags_unclassified_missing_frac,
                outputs=self.output_mags,
                metadata=self.metadata1,
                where="timepoint>30",
                exclude_ids=True
            )
            output = mock_stdout.getvalue().strip()

        self.assertIn("The filter query returned no IDs to filter out.", output)

    def test_filter_kraken_reports_empty(self):
        results = _filter_kraken2_results_by_metadata(
            reports=self.report_mags_unclassified_missing_frac,
            outputs=self.output_mags,
            remove_empty=True,
        )
        for i, suffix in zip(range(2), ["report", "output"]):
            self.assertTrue(
                os.path.exists(
                    os.path.join(
                        str(results[i]),
                        f"3b72d1a7-ddb0-4dc7-ac36-080ceda04aaa.{suffix}.txt"
                    )
                )
            )

    def test_missing_metadata_and_remove_empty(self):
        with self.assertRaisesRegex(
                ValueError, r'--m-metadata-file.*--p-remove-empty'
        ):
            _validate_parameters(metadata=None, remove_empty=False,
                                 where=None, exclude_ids=False)

    def test_where_without_metadata(self):
        with self.assertRaisesRegex(
                ValueError, r'--p-where.*--m-metadata-file'
        ):
            _validate_parameters(metadata=None, remove_empty=True,
                                 where=True, exclude_ids=False)

    def test_exclude_ids_without_metadata(self):
        with self.assertRaisesRegex(
                ValueError, r'--p-exclude-ids.*--m-metadata-file'
        ):
            _validate_parameters(metadata=None, remove_empty=True,
                                 where=None, exclude_ids=True)

    def test_validate_ids(self):
        ids = _validate_ids(
            file_dict_reports={"": self.file_dict_report_mags},
            file_dict_outputs={"": self.file_dict_output_mags},
        )
        self.assertTrue(ids == {"3b72d1a7-ddb0-4dc7-ac36-080ceda04aaa",
                                "8894435a-c836-4c18-b475-8b38a9ab6c6b"}
                        )

    def test_validate_ids_error(self):
        file_dict_reports = self.file_dict_report_mags
        file_dict_outputs = self.file_dict_output_mags
        del self.file_dict_report_mags["3b72d1a7-ddb0-4dc7-ac36-080ceda04aaa"]
        del self.file_dict_output_mags["8894435a-c836-4c18-b475-8b38a9ab6c6b"]
        msg = (r"reports\: \{'3b72d1a7-ddb0-4dc7-ac36-080ceda04aaa'\}.*\n.*outputs\: "
               r"\{'8894435a-c836-4c18-b475-8b38a9ab6c6b'\}.*")
        with self.assertRaisesRegex(ValueError, msg):
            _validate_ids(
                file_dict_reports={"": file_dict_reports},
                file_dict_outputs={"": file_dict_outputs},
            )


class TestAbundanceFilter(TestPluginBase):
    package = "q2_annotate.kraken2.tests"

    def setUp(self):
        super().setUp()

        curated_report_fp = self.get_data_path(
            'filter/reports/sample-1.report.txt'
        )
        curated_report = Kraken2ReportFormat(curated_report_fp, mode='r')
        self.curated_report = curated_report.view(pd.DataFrame)

        real_report_fp = self.get_data_path(
            'filter/reports/SRR17001003.report.txt'
        )
        real_report = Kraken2ReportFormat(real_report_fp, mode='r')
        self.real_report = real_report.view(pd.DataFrame)

        reports_fp = self.get_data_path('filter/reports')
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

    def test_curated_report_round_trip(self):
        '''
        Test that parsing and dumping a curated report results in no changes.
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

    def test_real_report_round_trip(self):
        '''
        Test that parsing and dumping a real report results in no changes.
        '''
        root, unclassified_node = _report_df_to_tree(self.real_report)

        round_tripped_report = _dump_tree_to_report(root, unclassified_node)

        round_tripped_report.sort_values(
            'taxon_id', inplace=True, ascending=False
        )
        real_report = self.real_report.sort_values('taxon_id', ascending=False)

        columns = [
            'n_frags_covered', 'n_frags_assigned', 'taxon_id', 'name', 'rank'
        ]
        for column in columns:
            assert_series_equal(
                round_tripped_report[column],
                real_report[column],
                check_dtype=False,
                check_index=False,
            )

    def test_filter_kraken2_reports_by_abundance(self):
        '''
        Test that the main `_filter_kraken2_reports_by_abundance` method runs,
        results in the same number of outputted formats as inputted ones.
        '''
        filtered_reports = _filter_kraken2_reports_by_abundance(
            self.reports, abundance_threshold=0.01
        )

        num_input_reports = len(list(
            self.reports.reports.iter_views(Kraken2ReportFormat)
        ))
        num_output_reports = len(list(
            filtered_reports.reports.iter_views(Kraken2ReportFormat)
        ))

        self.assertEqual(num_input_reports, num_output_reports)


class TestOutputReportAlignment(TestPluginBase):
    package = "q2_annotate.kraken2.tests"

    def setUp(self):
        super().setUp()

        report_fp = self.get_data_path(
            'filter/output-report-alignment/reports/sample-1.report.txt'
        )
        self.report = Kraken2ReportFormat(report_fp, mode='r')

        output_fp = self.get_data_path(
            'filter/output-report-alignment/outputs/sample-1.output.txt'
        )
        self.output = Kraken2OutputFormat(output_fp, mode='r')

        reports_fp = self.get_data_path(
            'filter/output-report-alignment/reports/'
        )
        self.reports = Kraken2ReportDirectoryFormat(reports_fp, mode='r')

        outputs_fp = self.get_data_path(
            'filter/output-report-alignment/outputs/'
        )
        self.outputs = Kraken2OutputDirectoryFormat(outputs_fp, mode='r')

    def test_align_single_output_with_report(self):
        '''
        Tests that a kraken2 report that has a subset of the taxon ids present
        in a kraken2 output results in a properly filtered output when the two
        are aligned.
        '''
        output_dir_fmt = Kraken2OutputDirectoryFormat()
        sample_id = 'sample-1'

        _align_single_output_with_report(
            self.output, self.report, output_dir_fmt, sample_id
        )

        _, aligned_output = list(
            output_dir_fmt.outputs.iter_views(Kraken2OutputFormat)
        )[0]

        report_df = self.report.view(pd.DataFrame)
        output_df = aligned_output.view(pd.DataFrame)

        # filtered taxon ids have been removed
        self.assertEqual(
            set(report_df['taxon_id']), set(output_df['taxon_id'])
        )

        # duplicate taxon id records still present in report
        self.assertGreater(len(output_df), len(report_df))

    def test_align_outputs_with_reports(self):
        '''
        Tests that reports and outputs are properly paired when performing
        alignment on entire report/output directory formats.
        '''
        aligned_outputs = _align_outputs_with_reports(
            self.outputs, self.reports
        )

        sample1_report_df = Kraken2ReportFormat(
            self.reports.path / 'sample-1.report.txt', mode='r'
        ).view(pd.DataFrame)
        sample2_report_df = Kraken2ReportFormat(
            self.reports.path / 'sample-2.report.txt', mode='r'
        ).view(pd.DataFrame)
        sample1_output_df = Kraken2OutputFormat(
            aligned_outputs.path / 'sample-1.output.txt', mode='r'
        ).view(pd.DataFrame)
        sample2_output_df = Kraken2OutputFormat(
            aligned_outputs.path / 'sample-2.output.txt', mode='r'
        ).view(pd.DataFrame)

        # reports, outputs paired properly
        self.assertEqual(
            set(sample1_report_df['taxon_id']),
            set(sample1_output_df['taxon_id'])
        )
        self.assertEqual(
            set(sample2_report_df['taxon_id']),
            set(sample2_output_df['taxon_id'])
        )


class TestFilterKraken2ResultsPipeline(TestPluginBase):
    package = "q2_annotate.kraken2.tests"

    @classmethod
    def setUpClass(cls):
        pm = qiime2.sdk.PluginManager()
        cls.filter_kraken2_results = pm.plugins['annotate'].pipelines[
            'filter_kraken2_results'
        ]

    def setUp(self):
        super().setUp()

        self.mag_reports = qiime2.Artifact.import_data(
            "FeatureData[Kraken2Report % Properties('mags')]",
            self.get_data_path('reports-mags')
        )
        self.mag_outputs = qiime2.Artifact.import_data(
            "FeatureData[Kraken2Output % Properties('mags')]",
            self.get_data_path('outputs-mags')
        )

    def test_empty_filter_with_mags(self):
        filtered_reports, filtered_outptus = self.filter_kraken2_results(
            self.mag_reports, self.mag_outputs, remove_empty=True
        )

    def test_abundance_filter_with_mags(self):
        pre_row_count = 0
        pre_reports = self.mag_reports.view(Kraken2ReportDirectoryFormat)
        for _, report in pre_reports.reports.iter_views(Kraken2ReportFormat):
            pre_row_count += report.view(pd.DataFrame).shape[0]

        filtered_reports, filtered_outptus = self.filter_kraken2_results(
            self.mag_reports, self.mag_outputs, abundance_threshold=0.05
        )

        post_row_count = 0
        post_reports = filtered_reports.view(Kraken2ReportDirectoryFormat)
        for _, report in post_reports.reports.iter_views(Kraken2ReportFormat):
            post_row_count += report.view(pd.DataFrame).shape[0]

        self.assertGreater(pre_row_count, post_row_count)

    def test_errors_on_no_filtering_parameters(self):
        with self.assertRaisesRegex(ValueError, "None of.*filtering"):
            filtered_reports, filtered_outptus = self.filter_kraken2_results(
                self.mag_reports, self.mag_outputs
            )

    def test_errors_on_mismatched_sample_ids(self):
        mag_reports = qiime2.Artifact.import_data(
            "FeatureData[Kraken2Report % Properties('mags')]",
            self.get_data_path('reports-mags')
        )
        with tempfile.TemporaryDirectory() as tempdir:
            filename = '3b72d1a7-ddb0-4dc7-ac36-080ceda04aaa.output.txt'
            shutil.copyfile(
                Path(self.get_data_path('outputs-mags')) / filename,
                Path(tempdir) / filename
            )

            mag_outputs = qiime2.Artifact.import_data(
                "FeatureData[Kraken2Output % Properties('mags')]", tempdir
            )

            with self.assertRaisesRegex(ValueError, "sample IDs"):
                self.filter_kraken2_results(
                    mag_reports, mag_outputs, remove_empty=True
                )
