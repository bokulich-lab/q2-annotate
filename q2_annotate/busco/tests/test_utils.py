# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json
import tempfile
import pandas as pd
from qiime2.plugin.testing import TestPluginBase
from q2_annotate.busco.utils import (
    _parse_busco_params, _collect_summaries, _parse_df_columns,
    _partition_dataframe, _get_feature_table, _calculate_summary_stats,
    _get_mag_lengths, _validate_lineage_dataset_input
)
from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt
from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_annotate.busco.types import BuscoDatabaseDirFmt
from q2_annotate.busco.busco import calculate_unbinned_percentage, count_contigs, filter_unbinned_for_partition, get_fasta_files_from_dir
from q2_types.per_sample_sequences import ContigSequencesDirFmt
#?
from unittest.mock import patch, ANY, call, MagicMock

from q2_types.per_sample_sequences._methods import partition_sample_data_mags
from q2_assembly.filter import filter_contigs
from qiime2 import Metadata
from pathlib import Path
import warnings
from qiime2 import Artifact
class TestBUSCOUtils(TestPluginBase):
    package = "q2_annotate.busco.tests"

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.mags = MultiMAGSequencesDirFmt(
            path=self.get_data_path('mags'),
            mode="r",
        )
        self.feature_data_mags = MAGSequencesDirFmt(
            path=self.get_data_path('mags/sample1'),
            mode="r",
        )
        self.df1 = pd.DataFrame({
            'sample_id': ['sample1'] * 6 + ['sample2'] * 4 + ['sample3'] * 5,
            'mag_id': [f'mag{i}' for i in range(1, 16)],
            'value': range(15)
        })
        self.df2 = pd.DataFrame({
            'sample_id': ['sample1'] * 6 + ['sample2'] * 6 + ['sample3'] * 3,
            'mag_id': [f'mag{i}' for i in range(1, 16)],
            'value': range(15)
        })
        self.df3 = pd.DataFrame({
            'mag_id': ['mag1', 'mag2', 'mag3'],
            'sample_id': ['sample1', 'sample2', 'sample3'],
            'dataset': ['dataset1', 'dataset2', 'dataset3'],
            'single': [1, 2, 3],
            'duplicated': [4, 5, 6],
            'fragmented': [7, 8, 9],
            'missing': [10, 11, 12],
            'complete': [13, 14, 15],
            'n_markers': [16, 17, 18],
            'contigs_n50': [19, 20, 21],
            'percent_gaps': [22, 23, 24],
            'scaffolds': [25, 26, 27],
            'length': [28, 29, 30]
        })
        self.df4 = pd.DataFrame({
            'id': ['mag1', 'mag2', 'mag3'],
            'percent_gaps': ['10%', '20%', '30%'],
            'single': ['1.0', '2.0', '3.0'],
            'duplicated': ['4.0', '5.0', '6.0'],
            'fragmented': ['7.0', '8.0', '9.0'],
            'missing': ['10.0', '11.0', '12.0'],
            'complete': ['13.0', '14.0', '15.0'],
            'n_markers': ['16', '17', '18']
        })
        self.df5 = pd.DataFrame({
            'index': [0, 1, 2],
            'mag_id': ['mag1', 'mag2', 'mag3'],
            'percent_gaps': [10.0, 20.0, 30.0],
            'single': [1.0, 2.0, 3.0],
            'duplicated': [4.0, 5.0, 6.0],
            'fragmented': [7.0, 8.0, 9.0],
            'missing': [10.0, 11.0, 12.0],
            'complete': [13.0, 14.0, 15.0],
            'n_markers': [16, 17, 18]
        })

    def test_parse_busco_params_1(self):
        observed = _parse_busco_params("auto_lineage", True)
        expected = ["--auto-lineage"]
        self.assertSetEqual(set(observed), set(expected))

    def test_parse_busco_params_2(self):
        observed = _parse_busco_params("evalue", 0.66)
        expected = ["--evalue", str(0.66)]
        self.assertSetEqual(set(observed), set(expected))

    def test_parse_busco_params_3(self):
        observed = _parse_busco_params("augustus", True)
        expected = ["--augustus"]
        self.assertSetEqual(set(observed), set(expected))

    def test_parse_busco_params_4(self):
        observed = _parse_busco_params("lineage_dataset", "bacteria-XYZ")
        expected = ["--lineage_dataset", "bacteria-XYZ"]
        self.assertSetEqual(set(observed), set(expected))

    def test_collect_summaries(self):
        with tempfile.TemporaryDirectory():
            paths = {}

            for i in range(1, 4):
                paths[f"sample{i}"] = self.get_data_path(
                    filename=f"batch_summary_sample{i}.txt"
                )

            obs = _collect_summaries(paths)
            exp = pd.read_csv(
                self.get_data_path(filename="all_batch_summaries.csv")
            )
            pd.set_option('display.max_columns', None)
            pd.testing.assert_frame_equal(obs, exp)

    def test_parse_df_columns(self):
        obs = _parse_df_columns(self.df4)
        exp = self.df5
        pd.testing.assert_frame_equal(obs, exp)

    def test_partition_dataframe_sample_data_max_rows_5(self):
        partitions = _partition_dataframe(self.df1, 5, True)
        self.assertEqual(len(partitions), 3)
        obs_shapes = [p.shape for p in partitions]
        exp_shapes = [(6, 3), (4, 3), (5, 3)]
        self.assertListEqual(obs_shapes, exp_shapes)

        partitions = _partition_dataframe(self.df2, 5, True)
        self.assertEqual(len(partitions), 3)
        obs_shapes = [p.shape for p in partitions]
        exp_shapes = [(6, 3), (6, 3), (3, 3)]
        self.assertListEqual(obs_shapes, exp_shapes)

    def test_partition_dataframe_sample_data_max_rows_10(self):
        partitions = _partition_dataframe(self.df1, 10, True)
        self.assertEqual(len(partitions), 2)
        obs_shapes = [p.shape for p in partitions]
        exp_shapes = [(10, 3), (5, 3)]
        self.assertListEqual(obs_shapes, exp_shapes)

        partitions = _partition_dataframe(self.df2, 10, True)
        self.assertEqual(len(partitions), 2)
        obs_shapes = [p.shape for p in partitions]
        exp_shapes = [(6, 3), (9, 3)]
        self.assertListEqual(obs_shapes, exp_shapes)

    def test_partition_dataframe_sample_data_max_rows_15(self):
        partitions = _partition_dataframe(self.df1, 15, True)
        self.assertEqual(len(partitions), 1)
        obs_shapes = [p.shape for p in partitions]
        exp_shapes = [(15, 3), ]
        self.assertListEqual(obs_shapes, exp_shapes)

        partitions = _partition_dataframe(self.df2, 15, True)
        self.assertEqual(len(partitions), 1)
        obs_shapes = [p.shape for p in partitions]
        exp_shapes = [(15, 3), ]
        self.assertListEqual(obs_shapes, exp_shapes)

    def test_partition_dataframe_feature_data_max_rows_5(self):
        n = 5
        df1 = self.df1.copy()
        df1 = df1.loc[df1["sample_id"] == "sample1"]
        partitions = _partition_dataframe(df1, n, False)
        self.assertEqual(len(partitions), 2)
        obs_shapes = [p.shape for p in partitions]
        exp_shapes = [(5, 3), (1, 3)]
        self.assertListEqual(obs_shapes, exp_shapes)

        df2 = self.df2.copy()
        df2 = df2.loc[df2["sample_id"] == "sample3"]
        partitions = _partition_dataframe(df2, n, False)
        self.assertEqual(len(partitions), 1)
        obs_shapes = [p.shape for p in partitions]
        exp_shapes = [(3, 3)]
        self.assertListEqual(obs_shapes, exp_shapes)

    def test_partition_dataframe_feature_data_max_rows_10(self):
        n = 10
        df1 = self.df1.copy()
        df1 = df1.loc[df1["sample_id"] == "sample1"]
        partitions = _partition_dataframe(df1, n, False)
        self.assertEqual(len(partitions), 1)
        obs_shapes = [p.shape for p in partitions]
        exp_shapes = [(6, 3)]
        self.assertListEqual(obs_shapes, exp_shapes)

        df2 = self.df2.copy()
        df2 = df2.loc[df2["sample_id"] == "sample2"]
        partitions = _partition_dataframe(df2, n, False)
        self.assertEqual(len(partitions), 1)
        obs_shapes = [p.shape for p in partitions]
        exp_shapes = [(6, 3)]
        self.assertListEqual(obs_shapes, exp_shapes)

    def test_partition_dataframe_feature_data_max_rows_15(self):
        n = 10
        df1 = self.df1.copy()
        df1 = df1.loc[df1["sample_id"] == "sample1"]
        partitions = _partition_dataframe(df1, n, False)
        self.assertEqual(len(partitions), 1)
        obs_shapes = [p.shape for p in partitions]
        exp_shapes = [(6, 3)]
        self.assertListEqual(obs_shapes, exp_shapes)

        df2 = self.df2.copy()
        df2 = df2.loc[df2["sample_id"] == "sample2"]
        partitions = _partition_dataframe(df2, n, False)
        self.assertEqual(len(partitions), 1)
        obs_shapes = [p.shape for p in partitions]
        exp_shapes = [(6, 3)]
        self.assertListEqual(obs_shapes, exp_shapes)

    def test_get_feature_table_sample_data(self):
        obs = json.loads(
            _get_feature_table(self.df3)
        )
        with open(
            self.get_data_path('feature_table_sample_data.json'), 'r'
        ) as f:
            exp = json.load(f)
        self.assertDictEqual(obs, exp)

    def test_get_feature_table_feature_data(self):
        df3 = self.df3.copy()
        df3 = df3.loc[df3["sample_id"] == "sample1"]
        obs = json.loads(
            _get_feature_table(df3)
        )
        with open(
            self.get_data_path('feature_table_feature_data.json'), 'r'
        ) as f:
            exp = json.load(f)
        self.assertDictEqual(obs, exp)

    def test_calculate_summary_stats(self):
        obs = _calculate_summary_stats(self.df3)
        exp = pd.DataFrame({
            "min": pd.Series({
                'single': 1,
                'duplicated': 4,
                'fragmented': 7,
                'missing': 10,
                'complete': 13
            }),
            "median": pd.Series({
                'single': 2.0,
                'duplicated': 5.0,
                'fragmented': 8.0,
                'missing': 11.0,
                'complete': 14.0
            }),
            "mean": pd.Series({
                'single': 2.0,
                'duplicated': 5.0,
                'fragmented': 8.0,
                'missing': 11.0,
                'complete': 14.0
            }),
            "max": pd.Series({
                'single': 3,
                'duplicated': 6,
                'fragmented': 9,
                'missing': 12,
                'complete': 15
            }),
            "count": pd.Series({
                'single': 3,
                'duplicated': 3,
                'fragmented': 3,
                'missing': 3,
                'complete': 3
            })
        }).T.to_json(orient='table')

        self.assertEqual(obs, exp)

    def test_get_mag_lengths_sample_data(self):
        obs = _get_mag_lengths(self.mags)
        exp = pd.Series(
            {
                '24dee6fe-9b84-45bb-8145-de7b092533a1': 1935,
                'ca7012fc-ba65-40c3-84f5-05aa478a7585': 3000,
                'fb0bc871-04f6-486b-a10e-8e0cb66f8de3': 2000,
                'd65a71fa-4279-4588-b937-0747ed5d604d': 3000,
                'db03f8b6-28e1-48c5-a47c-9c65f38f7357': 2000,
                'fa4d7420-d0a4-455a-b4d7-4fa66e54c9bf': 3000
            }, name="length"
        )
        pd.testing.assert_series_equal(obs, exp)

    def test_get_mag_lengths_feature_data(self):
        obs = _get_mag_lengths(self.feature_data_mags)
        exp = pd.Series(
            {
                '24dee6fe-9b84-45bb-8145-de7b092533a1': 1935,
                'ca7012fc-ba65-40c3-84f5-05aa478a7585': 3000,
                'fb0bc871-04f6-486b-a10e-8e0cb66f8de3': 2000,
            }, name="length"
        )
        pd.testing.assert_series_equal(obs, exp)

    def test_validate_lineage_dataset_input_valid(self):
        # Give path to valid database
        p = self.get_data_path("busco_db")
        busco_db = BuscoDatabaseDirFmt(path=p, mode="r")
        _validate_lineage_dataset_input(
            lineage_dataset="lineage_1",
            auto_lineage=False,
            auto_lineage_euk=False,
            auto_lineage_prok=False,
            busco_db=busco_db,
            kwargs={}
        )

    def test_validate_lineage_dataset_input_invalid(self):
        # Give path to valid database
        p = self.get_data_path("busco_db")
        busco_db = BuscoDatabaseDirFmt(path=p, mode="r")

        with self.assertRaisesRegex(
            ValueError,
            "is not present in input database."
        ):
            # Run busco
            _validate_lineage_dataset_input(
                lineage_dataset="lineage2",
                auto_lineage=False,
                auto_lineage_euk=False,
                auto_lineage_prok=False,
                busco_db=busco_db,
                kwargs={}
            )

    def test_validate_lineage_dataset_input_warning(self):
        # Give path to valid database
        p = self.get_data_path("busco_db")
        busco_db = BuscoDatabaseDirFmt(path=p, mode="r")
        kwargs = {
            "auto_lineage": True,
            "auto_lineage_euk": False,
            "auto_lineage_prok": False
        }
        with self.assertWarnsRegex(
            Warning,
            "`--p-auto-lineage` flags will be ignored."
        ):
            # Run busco
            _validate_lineage_dataset_input(
                lineage_dataset="lineage_1",
                auto_lineage=True,
                auto_lineage_euk=False,
                auto_lineage_prok=False,
                busco_db=busco_db,
                kwargs=kwargs
            )

        self.assertDictEqual(
            kwargs,
            {
                "auto_lineage": False,
                "auto_lineage_euk": False,
                "auto_lineage_prok": False
            }
        )
    def test_count_binned_contigs(self):
        # Load synthetic test data with known number of contigs
        # bins = MultiMAGSequencesDirFmt(
        #     path=self.get_data_path("mags"), mode='r'
        # )
        sample_path  = Path(self.get_data_path("mags")) / "sample1"
        # Extract all FASTA files from the directory
        fasta_files = get_fasta_files_from_dir(sample_path)
        count = count_contigs(fasta_files)
        #2+3+2 + 3+2+3 = 15
        self.assertEqual(count, 7)

    def test_count_unbinned_contigs(self):
        # Load synthetic test data with known number of contigs
        # bins = MultiMAGSequencesDirFmt(
        #     path=self.get_data_path("unbinned"), mode='r'
        # )
        sample_path = Path(self.get_data_path("unbinned")) / "sample1_contigs.fa"

        count = count_contigs([sample_path])
        #3+2=5
        self.assertEqual(count, 3)  


    def test_calculate_unbinned_percentage(self):

        mag_sample_path  = Path(self.get_data_path("mags")) / "sample1"
        # Extract all FASTA files from the directory
        # mag_sample_files = [fp for fp in mag_sample_path.glob("*") if fp.suffix in {".fa", ".fasta", ".fna"}]
        mag_sample_files = get_fasta_files_from_dir(mag_sample_path)
        # mags_count = count_contigs(mags)
        unbinned_sample_path = Path(self.get_data_path("unbinned")) / "sample1_contigs.fa"
        percentage, count = calculate_unbinned_percentage(mag_sample_files, [unbinned_sample_path])

        # Type checks
        self.assertIsInstance(percentage, float)
        self.assertIsInstance(count, int)

        expected_count = 3
        expected_percentage = (3 / (3 + 7)) * 100  

        self.assertEqual(count, expected_count)
        self.assertEqual(percentage, expected_percentage)
 
    def test_no_unbinned(self):
        # mags_count = count_contigs(mags)
        mag_sample_path  = Path(self.get_data_path("mags")) / "sample1"
        # Extract all files from the directory
        mag_sample_files = get_fasta_files_from_dir(mag_sample_path)
        # mags_count = count_contigs(mags)
        unbinned_sample_path = Path(self.get_data_path("unbinned_empty")) / "sample1_contigs.fa"

        percentage, count = calculate_unbinned_percentage(mag_sample_files, [unbinned_sample_path])
        # Type and range checks
        self.assertIsInstance(percentage, float)
        self.assertIsInstance(count, int)
        self.assertEqual(count, 0)
        self.assertEqual(percentage, 0.0)

    # def test_partial_unbinned(self):
        
    #     # mags_count = count_contigs(mags)
    #     percentage, count = calculate_unbinned_percentage(mags, unbinned)

    #     expected_count = 3
    #     expected_percentage = (3 / (3 + 15)) * 100  

    #     self.assertEqual(count, expected_count)
    #     self.assertEqual(percentage, expected_percentage)
    
    # def test_no_contigs_at_all(self):
    #     mags = MultiMAGSequencesDirFmt(
    #         path=self.get_data_path('mags_empty'),  # no .fasta files
    #         mode='r'
    #     )
    #     unbinned = ContigSequencesDirFmt(
    #         path=self.get_data_path('unbinned_empty'),  # empty .fasta files
    #         mode='r'
    #     )
    #     # mags_count = count_contigs(mags)
    #     percentage, count = calculate_unbinned_percentage(mags, unbinned)

    #     self.assertEqual(count, 0)
    #     self.assertEqual(percentage, 0.0)

    def test_only_unbinned(self):
        # mags_count = count_contigs(mags)
        mag_sample_path  = Path(self.get_data_path("mags_empty")) / "sample1"
        # Extract all files from the directory
        mag_sample_files = get_fasta_files_from_dir(mag_sample_path)
        # mags_count = count_contigs(mags)
        unbinned_sample_path = Path(self.get_data_path("unbinned")) / "sample1_contigs.fa"

        percentage, count = calculate_unbinned_percentage(mag_sample_files, [unbinned_sample_path])
        expected_count = 3
        self.assertEqual(count, expected_count)
        self.assertEqual(percentage, 100)

    
    def test_filtered_unbinned_matches_partition_1_sample(self):
        mag_fmt = MultiMAGSequencesDirFmt(
            path=self.get_data_path("partition_1_sample"),
            mode="r"
        )
        partitioned_mags = Artifact.import_data('SampleData[MAGs]', mag_fmt)

        unbinned = ContigSequencesDirFmt(
            path=self.get_data_path("unbinned"),
            mode="r"
        )

        expected_metadata = Metadata(pd.DataFrame(index=pd.Index(['sample1'], name="ID")))
        expected_where = "ID IN ('sample1')"
        # Mock _filter_contigs
        mock_filter_contigs = MagicMock(return_value=("filtered_result",))

        # Call function under test
        filter_unbinned_for_partition(unbinned, partitioned_mags, mock_filter_contigs)

        # Check arguments passed to the mock action
        mock_filter_contigs.assert_called_once_with(contigs=unbinned, metadata=expected_metadata, where=expected_where)

    def test_filtered_unbinned_matches_partition_2_samples(self):
        mag_fmt = MultiMAGSequencesDirFmt(
            path=self.get_data_path("partition_2_samples"),
            mode="r"
        )
        partitioned_mags = Artifact.import_data('SampleData[MAGs]', mag_fmt)
        # Load unbinned contigs
        unbinned = ContigSequencesDirFmt(
            path=self.get_data_path("unbinned"),
            mode="r"
        )
        expected_metadata = Metadata(pd.DataFrame(index=pd.Index(['sample1', 'sample2'], name="ID")))
        expected_where = "ID IN ('sample1', 'sample2')"
        # Mock _filter_contigs
        # mock_filter_contigs = MagicMock(return_value="filtered_result")
        mock_filter_contigs = MagicMock(return_value=("filtered_result",))

        # Call function under test
        filter_unbinned_for_partition(unbinned, partitioned_mags, mock_filter_contigs)
        
        # Check arguments passed to the mock action
        mock_filter_contigs.assert_called_once_with(contigs=unbinned, metadata=expected_metadata, where=expected_where)

