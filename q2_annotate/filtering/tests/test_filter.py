# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import unittest

import pandas as pd
import qiime2

from qiime2.plugin.testing import TestPluginBase

from q2_annotate.busco.types import BUSCOResultsFormat
from q2_annotate.filtering.filter_mags import (
    _filter_ids,
    _filter_manifest,
    _mags_to_df,
    filter_derep_mags,
    filter_mags,
    _find_empty_mags,
    _validate_parameters,
)
from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt


class TestMAGFiltering(TestPluginBase):
    package = "q2_annotate.filtering.tests"

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.metadata_df = pd.DataFrame(
            {
                "sample_id": ["id1", "id1", "id2", "id3", "id3", "id3"],
                "complete": [20.0, 98.0, 68.5, 100.0, 98.0, 21.3],
            },
            index=pd.Index(["mag1", "mag2", "mag3", "mag4", "mag5", "mag6"], name="id"),
        )

    def setUp(self):
        super().setUp()
        self.manifest_df = pd.read_csv(self.get_data_path("MANIFEST"))
        self.manifest_df.set_index(["sample-id", "mag-id"], inplace=True)
        self.mag_data_dir = self.get_data_path("mags")
        self.mag_derep_data_dir = self.get_data_path("mags/sample2")
        self.mag_df = pd.DataFrame(
            {
                "sample_id": ["sample1", "sample1", "sample2", "sample2"],
                "mag_id": [
                    "24dee6fe-9b84-45bb-8145-de7b092533a1",
                    "28f9219f-83c5-42ea-a1d9-1e20df03c707",
                    "d65a71fa-4279-4588-b937-0747ed5d604d",
                    "db03f8b6-28e1-48c5-a47c-9c65f38f7357",
                ],
                "mag_fp": [
                    f"{self.mag_data_dir}/sample1/"
                    f"24dee6fe-9b84-45bb-8145-de7b092533a1.fasta",
                    f"{self.mag_data_dir}/sample1/"
                    f"28f9219f-83c5-42ea-a1d9-1e20df03c707.fasta",
                    f"{self.mag_data_dir}/sample2/"
                    f"d65a71fa-4279-4588-b937-0747ed5d604d.fasta",
                    f"{self.mag_data_dir}/sample2/"
                    f"db03f8b6-28e1-48c5-a47c-9c65f38f7357.fasta",
                ],
            }
        )
        self.metadata = qiime2.Metadata(self.metadata_df)

    def test_filter_ids_include(self):
        ids = {"mag1", "mag2", "mag3", "mag4", "mag5", "mag6"}
        obs = _filter_ids(ids, self.metadata, "complete>40", exclude_ids=False)
        exp = {"mag2", "mag3", "mag4", "mag5"}
        self.assertSetEqual(obs, exp)

    def test_filter_ids_include_none(self):
        ids = {"mag1", "mag2", "mag3", "mag4", "mag5", "mag6"}
        obs = _filter_ids(ids, self.metadata, "complete<10", exclude_ids=False)
        exp = ids
        self.assertSetEqual(obs, exp)

    def test_filter_ids_exclude(self):
        ids = {"mag1", "mag2", "mag3", "mag4", "mag5", "mag6"}
        obs = _filter_ids(ids, self.metadata, "complete>40", exclude_ids=True)
        exp = {"mag1", "mag6"}
        self.assertSetEqual(obs, exp)

    def test_filter_manifest_mags(self):
        ids = {"mag1", "mag2", "mag5"}
        obs = _filter_manifest(self.manifest_df, ids)
        exp = self.manifest_df[
            self.manifest_df.index.get_level_values("mag-id").isin(ids)
        ]
        pd.testing.assert_frame_equal(obs, exp)

    def test_mags_to_df_on_mag(self):
        mags = MultiMAGSequencesDirFmt(self.mag_data_dir, mode="r")
        obs = _mags_to_df(mags.sample_dict())
        exp = self.mag_df
        exp.set_index("mag_id", inplace=True)
        pd.testing.assert_frame_equal(obs, exp)

    def test_filter_derep_mags(self):
        mags = MAGSequencesDirFmt(self.mag_derep_data_dir, mode="r")
        metadata = BUSCOResultsFormat(
            self.get_data_path("metadata-derep.tsv"), mode="r"
        ).view(qiime2.Metadata)

        obs = filter_derep_mags(mags, metadata, where="complete>10")
        obs_features = obs.feature_dict()
        exp_features = ["db03f8b6-28e1-48c5-a47c-9c65f38f7357"]
        self.assertListEqual(list(obs_features.keys()), exp_features)

    def test_filter_derep_mags_no_mags_remain(self):
        mags = MAGSequencesDirFmt(self.mag_derep_data_dir, mode="r")
        metadata = BUSCOResultsFormat(
            self.get_data_path("metadata-derep.tsv"), mode="r"
        ).view(qiime2.Metadata)

        with self.assertRaisesRegex(ValueError, "No MAGs remain"):
            filter_derep_mags(mags, metadata, where="complete<1")

    def test_filter_mags_remove_empty(self):
        mags = MultiMAGSequencesDirFmt(self.mag_data_dir, mode="r")
        obs = filter_mags(mags, remove_empty=True)
        obs_dict = obs.sample_dict()
        self.assertEqual(obs_dict.keys(), {"sample1", "sample2"})
        self.assertEqual(
            obs_dict["sample1"].keys(), {"24dee6fe-9b84-45bb-8145-de7b092533a1"}
        )
        self.assertEqual(
            obs_dict["sample2"].keys(),
            {
                "d65a71fa-4279-4588-b937-0747ed5d604d",
                "db03f8b6-28e1-48c5-a47c-9c65f38f7357",
            },
        )

    def test_find_empty_mags(self):
        mags = MultiMAGSequencesDirFmt(self.mag_data_dir, mode="r")
        mags_df = _mags_to_df(mags.sample_dict())
        empty_mags = _find_empty_mags(mags_df)
        self.assertEqual(empty_mags, {"28f9219f-83c5-42ea-a1d9-1e20df03c707"})

    def test_filter_derep_mags_remove_empty(self):
        mags = MAGSequencesDirFmt(self.get_data_path("mags/sample1"), mode="r")
        obs = filter_derep_mags(mags, remove_empty=True)
        obs_dict = obs.feature_dict()
        self.assertEqual(obs_dict.keys(), {"24dee6fe-9b84-45bb-8145-de7b092533a1"})

    def test_filter_mags_features(self):
        mags = MultiMAGSequencesDirFmt(self.mag_data_dir, mode="r")
        metadata = BUSCOResultsFormat(
            self.get_data_path("metadata-derep.tsv"), mode="r"
        ).view(qiime2.Metadata)

        obs = filter_mags(mags, metadata, where="length<1000000", on="mag")
        obs_samples = obs.sample_dict()
        exp_samples = ["sample2"]
        self.assertListEqual(list(obs_samples.keys()), exp_samples)

        obs_features = [
            mag_id
            for feature_dict in obs_samples.values()
            for mag_id in feature_dict.keys()
        ]
        exp_features = [
            "d65a71fa-4279-4588-b937-0747ed5d604d",
            "db03f8b6-28e1-48c5-a47c-9c65f38f7357",
        ]
        self.assertListEqual(obs_features, exp_features)

    def test_filter_mags_samples(self):
        mags = MultiMAGSequencesDirFmt(self.mag_data_dir, mode="r")
        metadata = qiime2.Metadata(
            pd.read_csv(
                self.get_data_path("metadata-sample.tsv"), sep="\t", index_col=0
            )
        )

        obs = filter_mags(mags, metadata, where="metric<5", on="sample")
        obs_samples = obs.sample_dict()
        exp_samples = ["sample2"]
        self.assertListEqual(list(obs_samples.keys()), exp_samples)

        obs_feature_count = len(obs.sample_dict()["sample2"])
        exp_feature_count = 2
        self.assertEqual(obs_feature_count, exp_feature_count)

    def test_filter_mags_samples_empty(self):
        mags = MultiMAGSequencesDirFmt(self.mag_data_dir, mode="r")
        metadata = qiime2.Metadata(
            pd.read_csv(
                self.get_data_path("metadata-sample.tsv"), sep="\t", index_col=0
            )
        )

        obs = filter_mags(
            mags, metadata, where="metric>5", on="sample", remove_empty=True
        )
        obs_samples = obs.sample_dict()
        exp_samples = ["sample1"]
        self.assertListEqual(list(obs_samples.keys()), exp_samples)

        obs_feature_count = len(obs.sample_dict()["sample1"])
        exp_feature_count = 1
        self.assertEqual(obs_feature_count, exp_feature_count)

    def test_filter_mags_no_mags_remain(self):
        mags = MultiMAGSequencesDirFmt(self.mag_data_dir, mode="r")
        metadata = qiime2.Metadata(
            pd.read_csv(
                self.get_data_path("metadata-sample.tsv"), sep="\t", index_col=0
            )
        )
        with self.assertRaisesRegex(ValueError, "No MAGs remain after filtering"):
            filter_mags(
                mags, metadata, where="metric>50", on="sample", remove_empty=True
            )

    def test_validate_parameters_no_parameter_error(self):
        with self.assertRaisesRegex(ValueError, "At least one of the"):
            _validate_parameters(None, None)


if __name__ == "__main__":
    unittest.main()
