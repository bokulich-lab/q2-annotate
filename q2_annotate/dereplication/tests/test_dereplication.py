# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob

import os
import unittest

import pandas as pd
import qiime2
import skbio

from q2_annotate.dereplication.derep import (
    _find_similar_bins_fcluster,
    _get_bin_lengths,
    _remap_bins,
    _reassign_bins_to_samples,
    _write_unique_bins,
    _generate_pa_table,
    dereplicate_mags,
    _get_representatives,
)
from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt

from qiime2.plugin.testing import TestPluginBase


class TestDereplication(TestPluginBase):
    package = "q2_annotate.dereplication.tests"

    def setUp(self):
        super().setUp()
        self.bins = MultiMAGSequencesDirFmt(self.get_data_path("mags"), mode="r")
        self.dist_matrix_df = pd.read_csv(
            self.get_data_path("distance-matrix.tsv"), sep="\t", index_col=0
        )
        self.dist_matrix = skbio.DistanceMatrix(
            self.dist_matrix_df, ids=self.dist_matrix_df.index
        )
        self.clusters_99 = [
            ["24dee6fe-9b84-45bb-8145-de7b092533a1"],
            [
                "ca7012fc-ba65-40c3-84f5-05aa478a7585",
                "db03f8b6-28e1-48c5-a47c-9c65f38f7357",
                "fa4d7420-d0a4-455a-b4d7-4fa66e54c9bf",
            ],
            [
                "d65a71fa-4279-4588-b937-0747ed5d604d",
                "fb0bc871-04f6-486b-a10e-8e0cb66f8de3",
            ],
        ]
        # fmt: off
        self.bin_map = {
            "ca7012fc-ba65-40c3-84f5-05aa478a7585":
                "ca7012fc-ba65-40c3-84f5-05aa478a7585",
            "d65a71fa-4279-4588-b937-0747ed5d604d":
                "d65a71fa-4279-4588-b937-0747ed5d604d",
            "db03f8b6-28e1-48c5-a47c-9c65f38f7357":
                "ca7012fc-ba65-40c3-84f5-05aa478a7585",
            "fa4d7420-d0a4-455a-b4d7-4fa66e54c9bf":
                "ca7012fc-ba65-40c3-84f5-05aa478a7585",
            "fb0bc871-04f6-486b-a10e-8e0cb66f8de3":
                "d65a71fa-4279-4588-b937-0747ed5d604d",
            "24dee6fe-9b84-45bb-8145-de7b092533a1":
                "24dee6fe-9b84-45bb-8145-de7b092533a1",
        }
        # fmt: on
        self.dereplicated_bins = {
            "sample1": {
                "ca7012fc-ba65-40c3-84f5-05aa478a7585": 1,
                "d65a71fa-4279-4588-b937-0747ed5d604d": 1,
                "24dee6fe-9b84-45bb-8145-de7b092533a1": 1,
            },
            "sample2": {
                "ca7012fc-ba65-40c3-84f5-05aa478a7585": 2,
                "d65a71fa-4279-4588-b937-0747ed5d604d": 1,
                "24dee6fe-9b84-45bb-8145-de7b092533a1": 0,
            },
        }
        self.busco_results = pd.read_csv(
            filepath_or_buffer=self.get_data_path("busco_results.tsv"),
            sep="\t",
            index_col=0,
            dtype="str",
        )
        self.busco_results.index.name = "id"

    def test_find_clusters_fcluster_similar(self):
        obs = _find_similar_bins_fcluster(self.dist_matrix_df, 0.99)
        exp = self.clusters_99
        self.assertListEqual(exp, obs)

    def test_find_clusters_fcluster_different(self):
        obs = _find_similar_bins_fcluster(self.dist_matrix_df, 0.1)
        exp = [[x] for sublist in self.clusters_99 for x in sublist]
        self.assertListEqual(sorted(exp), sorted(obs))

    def test_bin_lengths(self):
        obs = _get_bin_lengths(self.bins)
        exp = pd.Series(
            [1935, 3000, 2000, 3000, 2000, 3000],
            name="length",
            index=[
                "24dee6fe-9b84-45bb-8145-de7b092533a1",
                "ca7012fc-ba65-40c3-84f5-05aa478a7585",
                "fb0bc871-04f6-486b-a10e-8e0cb66f8de3",
                "d65a71fa-4279-4588-b937-0747ed5d604d",
                "db03f8b6-28e1-48c5-a47c-9c65f38f7357",
                "fa4d7420-d0a4-455a-b4d7-4fa66e54c9bf",
            ],
        )
        pd.testing.assert_series_equal(exp, obs)

    def test_remap_bins(self):
        longest_bins = [
            "24dee6fe-9b84-45bb-8145-de7b092533a1",
            "ca7012fc-ba65-40c3-84f5-05aa478a7585",
            "d65a71fa-4279-4588-b937-0747ed5d604d",
        ]
        obs = _remap_bins(self.clusters_99, longest_bins, self.dist_matrix_df)
        exp = self.bin_map
        self.assertDictEqual(exp, obs)

    def test_write_unique_bins(self):
        obs = _write_unique_bins(self.bins, self.bin_map)
        exp = MAGSequencesDirFmt(self.get_data_path("mags-unique"), mode="r")
        self.assertIsInstance(obs, MAGSequencesDirFmt)

        obs_files = glob.glob(os.path.join(str(obs), "*.fasta"))
        exp_files = glob.glob(os.path.join(str(exp), "*.fasta"))
        self.assertSetEqual(
            set([os.path.basename(x) for x in exp_files]),
            set([os.path.basename(x) for x in obs_files]),
        )

    def test_reassign_bins_to_samples(self):
        obs = _reassign_bins_to_samples(
            self.bin_map, self.bins.manifest.view(pd.DataFrame)
        )
        exp = self.dereplicated_bins
        self.assertDictEqual(exp, obs)

    def test_generate_pa_table(self):
        obs = _generate_pa_table(self.dereplicated_bins)
        exp = pd.read_csv(self.get_data_path("pa-table.csv"), index_col=0)
        pd.testing.assert_frame_equal(exp, obs)

    def test_dereplicate_mags(self):
        mags = MultiMAGSequencesDirFmt(self.get_data_path("mags"), mode="r")

        obs_mags, obs_pa = dereplicate_mags(mags, self.dist_matrix, threshold=0.99)
        exp_mags = MAGSequencesDirFmt(self.get_data_path("mags-unique"), mode="r")
        exp_pa = pd.read_csv(self.get_data_path("pa-table.csv"), index_col=0)

        self.assertIsInstance(obs_mags, MAGSequencesDirFmt)

        # assert correct bins were formed
        obs_files = glob.glob(os.path.join(str(obs_mags), "*.fasta"))
        exp_files = glob.glob(os.path.join(str(exp_mags), "*.fasta"))
        self.assertSetEqual(
            set([os.path.basename(x) for x in exp_files]),
            set([os.path.basename(x) for x in obs_files]),
        )

        # assert correct PA table was generated
        pd.testing.assert_frame_equal(exp_pa, obs_pa)

    def test_get_representatives_metadata_max_value(self):
        obs = _get_representatives(
            mags=self.bins,
            metadata=qiime2.Metadata(self.busco_results),
            column="complete",
            bin_clusters=self.clusters_99,
            find_max=True,
        )
        exp = [
            "24dee6fe-9b84-45bb-8145-de7b092533a1",
            "fa4d7420-d0a4-455a-b4d7-4fa66e54c9bf",
            "d65a71fa-4279-4588-b937-0747ed5d604d",
        ]
        self.assertListEqual(exp, obs)

    def test_get_representatives_metadata_min_value(self):
        obs = _get_representatives(
            mags=self.bins,
            metadata=qiime2.Metadata(self.busco_results),
            column="complete",
            bin_clusters=self.clusters_99,
            find_max=False,
        )
        exp = [
            "24dee6fe-9b84-45bb-8145-de7b092533a1",
            "ca7012fc-ba65-40c3-84f5-05aa478a7585",
            "fb0bc871-04f6-486b-a10e-8e0cb66f8de3",
        ]
        self.assertListEqual(exp, obs)

    def test_get_representatives_length(self):
        obs = _get_representatives(
            mags=self.bins,
            metadata=None,
            column=None,
            bin_clusters=self.clusters_99,
            find_max=True,
        )
        exp = [
            "24dee6fe-9b84-45bb-8145-de7b092533a1",
            "ca7012fc-ba65-40c3-84f5-05aa478a7585",
            "d65a71fa-4279-4588-b937-0747ed5d604d",
        ]
        self.assertListEqual(exp, obs)

    def test_get_representatives_key_error(self):
        with self.assertRaisesRegex(KeyError, "version"):
            _get_representatives(
                mags=self.bins,
                metadata=qiime2.Metadata(self.busco_results),
                column="version",
                bin_clusters=self.clusters_99,
                find_max=True,
            )

    def test_get_representatives_value_error(self):
        with self.assertRaisesRegex(ValueError, "numerical"):
            _get_representatives(
                mags=self.bins,
                metadata=qiime2.Metadata(self.busco_results),
                column="dataset",
                bin_clusters=self.clusters_99,
                find_max=True,
            )


if __name__ == "__main__":
    unittest.main()
