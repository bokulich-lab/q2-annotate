# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json
import biom
import numpy as np
import pandas as pd
import scipy.sparse as sp
from qiime2.plugin.testing import TestPluginBase

from q2_annotate.tfa import (
    estimate_tfa,
    tfa_by_function,
    tfa_by_taxon,
)


class TestTFA(TestPluginBase):
    package = "q2_annotate.tfa.tests"

    def setUp(self):
        super().setUp()
        self.abundance = pd.DataFrame(
            {
                "C1": [100, 0, 0],
                "C2": [50, 0, 0],
                "C3": [0, 120, 0],
                "C4": [0, 40, 0],
                "C5": [0, 0, 90],
                "C6": [0, 0, 60],
            },
            index=["S1", "S2", "S3"],
        )
        self.inventory = pd.DataFrame(
            {
                "bla_TEM": [2, 1, 2, 1, 2, 1],
                "vanA": [1, 0, 1, 0, 1, 0],
                "mecA": [0, 3, 0, 3, 0, 3],
            },
            index=["C1", "C2", "C3", "C4", "C5", "C6"],
        )
        self.taxonomy = pd.DataFrame(
            {"Taxon": ["Escherichia coli", "Klebsiella pneumoniae"]},
            index=["T1", "T2"],
        )
        self.mapping = {"T1": ["C1", "C3", "C5"], "T2": ["C2", "C4", "C6"]}

    def estimate(self, abundance=None, inventory=None, taxonomy=None, mapping=None):
        abundance = self.abundance if abundance is None else abundance
        inventory = self.inventory if inventory is None else inventory
        taxonomy = self.taxonomy if taxonomy is None else taxonomy
        mapping = self.mapping if mapping is None else mapping
        abundance, inventory = abundance.T, inventory.T
        return estimate_tfa(
            biom.Table(
                abundance.values,
                observation_ids=abundance.index,
                sample_ids=abundance.columns,
            ),
            biom.Table(
                inventory.values,
                observation_ids=inventory.index,
                sample_ids=inventory.columns,
            ),
            taxonomy,
            mapping,
        )

    def test_sample_resolved_values_and_sparse_pairs(self):
        observed = self.estimate()
        self.assertEqual(list(observed.ids(axis="sample")), ["S1", "S2", "S3"])
        self.assertEqual(observed.shape, (4, 3))
        self.assertTrue(sp.issparse(observed.matrix_data))
        self.assertEqual(observed.matrix_data.nnz, 12)
        expected = {
            ("T1", "bla_TEM"): [200, 240, 180],
            ("T1", "vanA"): [100, 120, 90],
            ("T2", "bla_TEM"): [50, 40, 60],
            ("T2", "mecA"): [150, 120, 180],
        }
        for pair_id, values in zip(
            observed.ids(axis="observation"), observed.matrix_data.toarray()
        ):
            np.testing.assert_array_equal(values, expected[tuple(json.loads(pair_id))])

    def test_taxon_and_function_projections(self):
        observed = self.estimate()
        by_taxon = tfa_by_taxon(observed, ["T1"])["T1"]
        pd.testing.assert_frame_equal(
            by_taxon.to_dataframe(dense=True),
            pd.DataFrame(
                [[200, 240, 180], [100, 120, 90]],
                index=["bla_TEM", "vanA"],
                columns=["S1", "S2", "S3"],
            ),
            check_dtype=False,
        )
        by_function = tfa_by_function(observed, ["bla_TEM"])["bla_TEM"]
        pd.testing.assert_frame_equal(
            by_function.to_dataframe(dense=True),
            pd.DataFrame(
                [[200, 240, 180], [50, 40, 60]],
                index=["T1", "T2"],
                columns=["S1", "S2", "S3"],
            ),
            check_dtype=False,
        )
        self.assertTrue(sp.issparse(by_taxon.matrix_data))
        with self.assertRaisesRegex(ValueError, "No TFA rows"):
            tfa_by_taxon(observed, ["unknown"])

    def test_multiple_selections_only_return_requested_ids(self):
        observed = self.estimate()
        by_taxon = tfa_by_taxon(observed, ["T2", "T1", "T2"])
        self.assertEqual(list(by_taxon), ["T2", "T1"])
        self.assertEqual(
            list(by_taxon["T2"].ids(axis="observation")), ["bla_TEM", "mecA"]
        )
        by_function = tfa_by_function(observed, ["vanA", "mecA"])
        self.assertEqual(list(by_function), ["vanA", "mecA"])
        self.assertEqual(list(by_function["vanA"].ids(axis="observation")), ["T1"])
        with self.assertRaisesRegex(ValueError, "at least one taxon"):
            tfa_by_taxon(observed, [])
        with self.assertRaisesRegex(ValueError, "unknown"):
            tfa_by_function(observed, ["vanA", "unknown"])

    def test_no_common_contigs_retains_samples(self):
        observed = self.estimate(mapping={"T1": ["not-a-contig"]})
        self.assertEqual(observed.shape, (0, 3))
        self.assertEqual(list(observed.ids(axis="sample")), ["S1", "S2", "S3"])

    def test_ignores_extra_contigs(self):
        abundance = self.abundance.assign(C7=[1000, 1000, 1000])
        taxonomy = pd.DataFrame(
            {"Taxon": ["E. coli", "K. pneumoniae", "Other"]},
            index=["T1", "T2", "T3"],
        )
        mapping = {**self.mapping, "T3": ["C7", "C8"]}
        observed = self.estimate(
            abundance=abundance, taxonomy=taxonomy, mapping=mapping
        )
        self.assertEqual(observed.shape, (4, 3))
        self.assertEqual(observed.matrix_data.nnz, 12)

    def test_zero_abundance_keeps_observed_pairs(self):
        observed = self.estimate(abundance=self.abundance * 0)
        self.assertEqual(observed.shape, (4, 3))
        self.assertEqual(observed.matrix_data.nnz, 0)

    def test_extreme_abundances(self):
        abundance = self.abundance * 0.0
        abundance.loc["S1", "C1"] = 1e-9
        abundance.loc["S1", "C2"] = 1e12
        observed = self.estimate(abundance=abundance)
        np.testing.assert_allclose(
            tfa_by_taxon(observed, ["T1"])["T1"].matrix_data.toarray(),
            [[2e-9, 0, 0], [1e-9, 0, 0]],
            rtol=1e-12,
            atol=0,
        )
        np.testing.assert_allclose(
            tfa_by_taxon(observed, ["T2"])["T2"].matrix_data.toarray(),
            [[1e12, 0, 0], [3e12, 0, 0]],
            rtol=1e-12,
            atol=0,
        )

    def test_pair_ids_round_trip_ambiguous_characters(self):
        inventory = self.inventory.rename(columns={"bla_TEM": 'a,"b'})
        taxonomy = self.taxonomy.rename(index={"T1": 't,"1'})
        mapping = {'t,"1': self.mapping["T1"], "T2": self.mapping["T2"]}
        observed = self.estimate(
            inventory=inventory, taxonomy=taxonomy, mapping=mapping
        )
        self.assertIn(
            'a,"b',
            next(iter(tfa_by_taxon(observed, ['t,"1']).values())).ids(
                axis="observation"
            ),
        )
        self.assertIn(
            't,"1',
            next(iter(tfa_by_function(observed, ['a,"b']).values())).ids(
                axis="observation"
            ),
        )

    def test_invalid_pair_id_rejected(self):
        malformed = biom.Table(
            np.array([[1]]), observation_ids=["not-a-pair"], sample_ids=["S1"]
        )
        with self.assertRaisesRegex(ValueError, "Invalid TFA feature ID"):
            tfa_by_taxon(malformed, ["T1"])

    def test_ambiguous_taxon_assignment_rejected(self):
        mapping = {**self.mapping, "T2": ["C1", "C2", "C4", "C6"]}
        with self.assertRaisesRegex(ValueError, "multiple taxa"):
            self.estimate(mapping=mapping)
