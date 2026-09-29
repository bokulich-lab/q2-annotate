# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import biom
import numpy as np
import pandas as pd
import qiime2
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_table import FeatureTable, Frequency

from q2_annotate.plugin_setup import plugin
from q2_annotate.tfa.tfa import _collection_key


class TestTFAAction(TestPluginBase):
    package = "q2_annotate.tfa.tests"

    def test_action_registration(self):
        self.assertIn("estimate_tfa", plugin.methods)
        self.assertNotIn("_estimate_tfa_table", plugin.methods)
        self.assertIn("tfa_by_taxon", plugin.methods)
        self.assertIn("tfa_by_function", plugin.methods)
        self.assertIn("tfa_collapse_samples", plugin.methods)
        self.assertNotIn("estimate_tfa", plugin.pipelines)
        self.assertNotIn("TFA", plugin.type_fragments)

    def _inputs(self, inventory_contigs=("C1", "C2", "C3")):
        abundance = qiime2.Artifact.import_data(
            "FeatureTable[Frequency]",
            biom.Table(
                np.array([[2, 0], [1, 3], [0, 4]]),
                observation_ids=["C1", "C2", "C3"],
                sample_ids=["S1", "S2"],
            ),
        )
        inventory = qiime2.Artifact.import_data(
            "FeatureTable[Frequency]",
            biom.Table(
                np.array([[2, 1, 0], [0, 1, 3]])[:, : len(inventory_contigs)],
                observation_ids=["geneA", "gene B"],
                sample_ids=list(inventory_contigs),
            ),
        )
        taxonomy = qiime2.Artifact.import_data(
            "FeatureData[Taxonomy]",
            pd.DataFrame(
                {"Taxon": ["k__Bacteria; p__A", "k__Bacteria; p__B"]},
                index=pd.Index(["T1", "T2"], name="Feature ID"),
            ),
        )
        mapping = qiime2.Artifact.import_data(
            "FeatureMap[TaxonomyToContigs]",
            {"T1": ["C1", "C2"], "T2": ["C3"]},
        )
        return dict(
            contig_abundance=abundance,
            feature_inventory=inventory,
            taxonomy=taxonomy,
            taxon_to_contig_map=mapping,
        )

    def test_estimate_returns_frequency_and_gene_taxonomy(self):
        result = plugin.methods["estimate_tfa"](**self._inputs())
        self.assertEqual(len(result), 2)
        self.assertEqual(
            str(result.tfa.type), "FeatureTable[Frequency % Properties('tfa')]"
        )
        self.assertLessEqual(result.tfa.type, FeatureTable[Frequency])
        self.assertEqual(
            str(result.gene_taxonomy.type),
            "FeatureData[Taxonomy % Properties('tfa')]",
        )
        self.assertEqual(
            list(result.tfa.view(biom.Table).ids(axis="sample")),
            ["S1", "S2"],
        )
        self.assertEqual(
            list(result.tfa.view(biom.Table).ids(axis="observation")),
            list(result.gene_taxonomy.view(pd.DataFrame).index),
        )

    def test_no_contig_overlap_reports_empty_result(self):
        """Reject estimation when no pair has a nonzero sample load."""
        with self.assertRaisesRegex(
            ValueError, "No taxon/gene pairs have nonzero loads"
        ):
            plugin.methods["estimate_tfa"](**self._inputs(("C4",)))

    def test_tfa_actions_require_property(self):
        feature_load, gene_taxonomy = plugin.methods["estimate_tfa"](**self._inputs())
        unmarked = qiime2.Artifact.import_data(
            "FeatureTable[Frequency]", feature_load.view(biom.Table)
        )
        for action, parameters in (
            ("tfa_by_taxon", {"taxon_ids": ["T1"]}),
            ("tfa_by_function", {"function_ids": ["geneA"]}),
            ("tfa_collapse_samples", {}),
        ):
            with self.assertRaisesRegex(TypeError, "Properties.*tfa"):
                plugin.methods[action](
                    feature_load=unmarked, gene_taxonomy=gene_taxonomy, **parameters
                )

    def test_public_selectors_return_requested_sparse_tables(self):
        feature_load, gene_taxonomy = plugin.methods["estimate_tfa"](**self._inputs())
        by_taxon = plugin.methods["tfa_by_taxon"](
            feature_load=feature_load,
            gene_taxonomy=gene_taxonomy,
            taxon_ids=["T2", "T1"],
        ).tables
        self.assertEqual(set(by_taxon.keys()), {"T1", "T2"})
        by_taxon = by_taxon["T1"].view(biom.Table)
        self.assertEqual(list(by_taxon.ids(axis="observation")), ["geneA", "gene B"])
        np.testing.assert_array_equal(by_taxon.matrix_data.toarray(), [[5, 3], [1, 3]])
        by_function = plugin.methods["tfa_by_function"](
            feature_load=feature_load,
            gene_taxonomy=gene_taxonomy,
            function_ids=["gene B"],
        ).tables
        self.assertEqual(set(by_function.keys()), {_collection_key("gene B")})
        by_function = by_function[_collection_key("gene B")].view(biom.Table)
        self.assertEqual(list(by_function.ids(axis="observation")), ["T1", "T2"])
        np.testing.assert_array_equal(
            by_function.matrix_data.toarray(), [[1, 3], [0, 12]]
        )
        self.assertEqual(by_function.matrix_data.nnz, 3)

    def test_collapse_samples_returns_taxon_by_function_artifact(self):
        feature_load, gene_taxonomy = plugin.methods["estimate_tfa"](**self._inputs())
        result = plugin.methods["tfa_collapse_samples"](
            feature_load=feature_load, gene_taxonomy=gene_taxonomy
        )
        self.assertEqual(str(result.table.type), "FeatureTable[Frequency]")
        table = result.table.view(biom.Table)
        self.assertEqual(list(table.ids(axis="observation")), ["T1", "T2"])
        self.assertEqual(list(table.ids(axis="sample")), ["geneA", "gene B"])
        np.testing.assert_array_equal(table.matrix_data.toarray(), [[8, 4], [0, 12]])
        mean = plugin.methods["tfa_collapse_samples"](
            feature_load=feature_load, gene_taxonomy=gene_taxonomy, aggregation="mean"
        ).table.view(biom.Table)
        np.testing.assert_array_equal(mean.matrix_data.toarray(), [[4, 2], [0, 6]])
        median = plugin.methods["tfa_collapse_samples"](
            feature_load=feature_load, gene_taxonomy=gene_taxonomy, aggregation="median"
        ).table.view(biom.Table)
        np.testing.assert_array_equal(median.matrix_data.toarray(), [[4, 2], [0, 6]])

    def test_collection_keys_are_safe_and_distinct(self):
        self.assertEqual(_collection_key("T1"), "T1")
        self.assertNotEqual(_collection_key("id_T1"), "id_T1")
        self.assertNotEqual(_collection_key("gene B"), _collection_key("gene_B"))
        self.assertLessEqual(len(_collection_key("x" * 300)), 200)
