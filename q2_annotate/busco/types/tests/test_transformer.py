# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pandas as pd
import qiime2
from qiime2.plugin.testing import TestPluginBase
from q2_annotate.busco.types import BUSCOResultsFormat
from q2_annotate.busco.types._transformer import _read_dataframe


class TestBUSCOTransformers(TestPluginBase):
    package = "q2_annotate.busco.types.tests"

    def setUp(self):
        super().setUp()
        self.fp = self.get_data_path("busco_results.tsv")

    def test_read_dataframe(self):
        obs = _read_dataframe(self.fp)

        self.assertIsInstance(obs, pd.DataFrame)
        self.assertEqual(obs.shape, (3, 16))
        self.assertEqual(obs.index.name, "id")

    def test_result_to_dataframe_transformer(self):
        transformer = self.get_transformer(BUSCOResultsFormat, pd.DataFrame)
        obs = transformer(BUSCOResultsFormat(self.fp, mode="r"))
        exp = pd.read_csv(self.fp, sep="\t", header=0, index_col=0, dtype="str")
        exp.index.name = "id"

        pd.testing.assert_frame_equal(obs, exp)

    def test_dataframe_to_result_transformer(self):
        transformer = self.get_transformer(pd.DataFrame, BUSCOResultsFormat)
        df = pd.read_csv(self.fp, sep="\t", header=0, index_col=False, dtype="str")
        df.index.name = "id"
        obs = transformer(df)

        obs.validate()
        self.assertIsInstance(obs, BUSCOResultsFormat)

    def test_result_to_metadata_transformer(self):
        transformer = self.get_transformer(BUSCOResultsFormat, qiime2.Metadata)
        obs = transformer(BUSCOResultsFormat(self.fp, mode="r"))

        df = pd.read_csv(self.fp, sep="\t", header=0, index_col=0, dtype="str")
        df.index.name = "id"
        for col in [
            "complete",
            "single",
            "duplicated",
            "fragmented",
            "missing",
            "n_markers",
            "scaffold_n50",
            "contigs_n50",
            "scaffolds",
            "length",
            "completeness",
            "contamination",
        ]:
            df[col] = pd.to_numeric(df[col])
        exp = qiime2.Metadata(df)

        self.assertEqual(obs, exp)
