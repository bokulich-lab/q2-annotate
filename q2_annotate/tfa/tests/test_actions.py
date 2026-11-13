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
        self.assertNotIn("estimate_tfa", plugin.methods)
        self.assertNotIn("_estimate_tfa_table", plugin.methods)
        self.assertIn("tfa_by_taxon", plugin.methods)
        self.assertIn("tfa_by_function", plugin.methods)
        self.assertIn("tfa_collapse_samples", plugin.methods)
        self.assertIn("estimate_tfa", plugin.pipelines)
        self.assertIn("_visualize_tfa", plugin.visualizers)
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

    def test_estimate_returns_frequency_gene_taxonomy_and_visualization(self):
        result = self.plugin.pipelines["estimate_tfa"](**self._inputs())
        self.assertEqual(len(result), 3)
        self.assertIsInstance(result.visualization, qiime2.Visualization)
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
            self.plugin.pipelines["estimate_tfa"](**self._inputs(("C4",)))

    def test_tfa_actions_require_property(self):
        result = self.plugin.pipelines["estimate_tfa"](**self._inputs())
        feature_load, gene_taxonomy = result.tfa, result.gene_taxonomy
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

    def test_by_taxon_returns_requested_sparse_tables(self):
        result = self.plugin.pipelines["estimate_tfa"](**self._inputs())
        feature_load, gene_taxonomy = result.tfa, result.gene_taxonomy
        by_taxon = plugin.methods["tfa_by_taxon"](
            feature_load=feature_load,
            gene_taxonomy=gene_taxonomy,
            taxon_ids=["T2", "T1"],
        ).tables
        self.assertEqual(set(by_taxon.keys()), {"T1", "T2"})
        by_taxon = by_taxon["T1"].view(biom.Table)
        self.assertEqual(list(by_taxon.ids(axis="observation")), ["geneA", "gene B"])
        np.testing.assert_array_equal(by_taxon.matrix_data.toarray(), [[5, 3], [1, 3]])

    def test_by_function_returns_requested_sparse_tables(self):
        result = self.plugin.pipelines["estimate_tfa"](**self._inputs())
        feature_load, gene_taxonomy = result.tfa, result.gene_taxonomy
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
        result = self.plugin.pipelines["estimate_tfa"](**self._inputs())
        feature_load, gene_taxonomy = result.tfa, result.gene_taxonomy
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

    def test_follow_up_actions_preserve_fractional_pipeline_loads_with_metadata(self):
        """Use metadata-aware pipeline outputs directly in all follow-up actions."""
        inputs = self._inputs()
        abundance = inputs["contig_abundance"].view(biom.Table)
        inputs["contig_abundance"] = qiime2.Artifact.import_data(
            "FeatureTable[Frequency]",
            biom.Table(
                abundance.matrix_data / 4,
                abundance.ids(axis="observation"),
                abundance.ids(axis="sample"),
            ),
        )
        metadata = qiime2.Metadata.load(self.get_data_path("sample-metadata.tsv"))
        result = self.plugin.pipelines["estimate_tfa"](**inputs, metadata=metadata)
        follow_up_inputs = dict(
            feature_load=result.tfa, gene_taxonomy=result.gene_taxonomy
        )
        by_taxon = self.plugin.methods["tfa_by_taxon"](
            **follow_up_inputs, taxon_ids=["T1"]
        ).tables["T1"].view(biom.Table)
        np.testing.assert_array_equal(
            by_taxon.matrix_data.toarray(), [[1.25, 0.75], [0.25, 0.75]]
        )
        by_function = self.plugin.methods["tfa_by_function"](
            **follow_up_inputs, function_ids=["gene B"]
        ).tables[_collection_key("gene B")].view(biom.Table)
        np.testing.assert_array_equal(
            by_function.matrix_data.toarray(), [[0.25, 0.75], [0, 3]]
        )
        for aggregation, expected in (
            ("sum", [[2, 1], [0, 3]]),
            ("mean", [[1, 0.5], [0, 1.5]]),
            ("median", [[1, 0.5], [0, 1.5]]),
        ):
            with self.subTest(aggregation=aggregation):
                collapsed = self.plugin.methods["tfa_collapse_samples"](
                    **follow_up_inputs, aggregation=aggregation
                ).table.view(biom.Table)
                np.testing.assert_array_equal(collapsed.matrix_data.toarray(), expected)

    def test_collection_keys_are_safe_and_distinct(self):
        self.assertEqual(_collection_key("T1"), "T1")
        self.assertNotEqual(_collection_key("id_T1"), "id_T1")
        self.assertNotEqual(_collection_key("gene B"), _collection_key("gene_B"))
        self.assertLessEqual(len(_collection_key("x" * 300)), 200)
