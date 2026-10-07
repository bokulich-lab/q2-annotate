# -----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
import json
import re
from pathlib import Path

import biom
import pandas as pd
import qiime2
import scipy.sparse as sp
from qiime2.plugin.testing import TestPluginBase

from q2_annotate.plugin_setup import plugin
from q2_annotate.tfa.tfa import _gene_taxonomy_id, estimate_tfa
from q2_annotate.tfa.visualizer import explore_tfa


class TestTFAVisualizer(TestPluginBase):
    package = "q2_annotate.tfa.tests"

    def _table(self):
        data = pd.read_csv(
            self.get_data_path("tfa-table.tsv"),
            sep="\t",
            dtype={"taxon_id": str, "function_id": str},
        )
        samples = list(data.columns[2:])
        return biom.Table(
            sp.csr_matrix(data[samples].to_numpy(dtype=float)),
            observation_ids=[
                _gene_taxonomy_id(taxon, function)
                for taxon, function in zip(data.taxon_id, data.function_id)
            ],
            sample_ids=samples,
        )

    def _mapping(self, filename="gene-taxonomy.tsv"):
        return pd.read_csv(
            self.get_data_path(filename),
            sep="\t",
            index_col=0,
            dtype=str,
            keep_default_na=False,
        )

    def _payload(self, html):
        match = re.search(
            r'<script id="tfa-data" type="application/json">(.*?)</script>',
            html,
            re.DOTALL,
        )
        self.assertIsNotNone(match)
        return json.loads(match.group(1))

    def test_sparse_payload_and_metadata_groups(self):
        metadata = qiime2.Metadata.load(self.get_data_path("sample-metadata.tsv"))
        explore_tfa(
            self.temp_dir.name, self._table(), self._mapping(), metadata=metadata
        )
        output = Path(self.temp_dir.name)
        html = (output / "index.html").read_text(encoding="utf-8")
        payload = self._payload(html)
        self.assertEqual(payload["samples"], ["S1", "S2", "S3"])
        self.assertEqual([pair["total"] for pair in payload["pairs"]], [8, 4.25, 16])
        self.assertEqual(
            [pair["mean"] for pair in payload["pairs"]], [8 / 3, 4.25 / 3, 16 / 3]
        )
        self.assertEqual([pair["median"] for pair in payload["pairs"]], [3, 1.25, 4])
        self.assertEqual(len(payload["loads"]), 6)
        self.assertEqual(payload["loads"][0], [0, 0, 5.0])
        self.assertEqual(
            payload["groups"],
            {
                "group": ["case", "control", "case"],
                "batch": ["B1", "B2", "B1"],
            },
        )
        self.assertIn("Explore taxon", html)
        self.assertTrue((output / "explore.js").is_file())
        self.assertTrue((output / "style.css").is_file())
        self.assertTrue((output / "taxonomy.js").is_file())
        self.assertIn('id="level"', html)
        self.assertEqual(
            payload["pairs"][0]["lineage"], ["Bacteria", "Escherichia coli"]
        )

    def test_estimated_zero_load_pairs_are_absent_from_explorer(self):
        """Show only taxon-gene pairs with a load in the explorer."""
        abundance = biom.Table(
            sp.csr_matrix([[2, 0], [0, 0]]),
            observation_ids=["C1", "C2"],
            sample_ids=["S1", "S2"],
        )
        inventory = biom.Table(
            sp.csr_matrix([[1, 0], [0, 1]]),
            observation_ids=["geneA", "geneB"],
            sample_ids=["C1", "C2"],
        )
        taxonomy = pd.DataFrame({"Taxon": ["Taxon A", "Taxon B"]}, index=["T1", "T2"])
        table, mapping = estimate_tfa(
            abundance, inventory, taxonomy, {"T1": ["C1"], "T2": ["C2"]}
        )
        explore_tfa(self.temp_dir.name, table, mapping)
        html = (Path(self.temp_dir.name) / "index.html").read_text(encoding="utf-8")
        pairs = self._payload(html)["pairs"]
        self.assertEqual(len(pairs), 1)
        self.assertEqual((pairs[0]["taxon_id"], pairs[0]["function"]), ("T1", "geneA"))

    def test_optional_metadata_and_html_escaping(self):
        table = biom.Table(
            sp.csr_matrix([[2]]),
            observation_ids=[_gene_taxonomy_id("T<script>", "gene")],
            sample_ids=["S1"],
        )
        explore_tfa(self.temp_dir.name, table, self._mapping("gene-taxonomy-html.tsv"))
        html = (Path(self.temp_dir.name) / "index.html").read_text(encoding="utf-8")
        self.assertNotIn("T<script>", html)
        self.assertEqual(self._payload(html)["pairs"][0]["taxon"], "T<script>")
        self.assertEqual(self._payload(html)["groups"], {})

    def test_even_sample_medians_include_sparse_zeros(self):
        table = self._table().filter(["S1", "S2"], axis="sample", inplace=False)
        explore_tfa(self.temp_dir.name, table, self._mapping())
        html = (Path(self.temp_dir.name) / "index.html").read_text(encoding="utf-8")
        self.assertEqual(
            [pair["median"] for pair in self._payload(html)["pairs"]], [4, 2.125, 6]
        )

    def test_mostly_absent_pairs_have_zero_median(self):
        original = self._table()
        table = biom.Table(
            sp.hstack([original.matrix_data, sp.csr_matrix((3, 3))]),
            observation_ids=original.ids(axis="observation"),
            sample_ids=["S1", "S2", "S3", "S4", "S5", "S6"],
        )
        explore_tfa(self.temp_dir.name, table, self._mapping())
        html = (Path(self.temp_dir.name) / "index.html").read_text(encoding="utf-8")
        pairs = self._payload(html)["pairs"]
        self.assertEqual([pair["median"] for pair in pairs], [0, 0, 0])
        self.assertEqual([pair["mean"] for pair in pairs], [8 / 6, 4.25 / 6, 16 / 6])

    def test_samples_missing_from_metadata_remain_grouped(self):
        metadata = qiime2.Metadata.load(
            self.get_data_path("sample-metadata-partial.tsv")
        )
        explore_tfa(
            self.temp_dir.name, self._table(), self._mapping(), metadata=metadata
        )
        html = (Path(self.temp_dir.name) / "index.html").read_text(encoding="utf-8")
        self.assertEqual(
            self._payload(html)["groups"]["group"],
            ["case", "control", "(missing)"],
        )

    def test_action_registration_and_execution(self):
        self.assertIn("explore_tfa", plugin.visualizers)
        artifact = qiime2.Artifact.import_data(
            "FeatureTable[Frequency % Properties('tfa')]", self._table()
        )
        mapping = qiime2.Artifact.import_data(
            "FeatureData[Taxonomy % Properties('tfa')]", self._mapping()
        )
        metadata = qiime2.Metadata.load(self.get_data_path("sample-metadata.tsv"))
        result = plugin.visualizers["explore_tfa"](
            feature_load=artifact, gene_taxonomy=mapping, metadata=metadata
        )
        self.assertIsInstance(result.visualization, qiime2.Visualization)
        without_metadata = plugin.visualizers["explore_tfa"](
            feature_load=artifact, gene_taxonomy=mapping
        )
        self.assertIsInstance(without_metadata.visualization, qiime2.Visualization)

    def test_mapping_supplies_full_and_short_taxonomy_labels(self):
        """Use mapping lineages for selectors and terminal names for heatmap axes."""
        mapping = self._mapping()
        mapping.loc[mapping["Taxon ID"] == "T1", "Taxon"] = (
            "d__Bacteria; g__Bacteroides; s__fragilis"
        )
        explore_tfa(self.temp_dir.name, self._table(), mapping)
        html = (Path(self.temp_dir.name) / "index.html").read_text(encoding="utf-8")
        pairs = self._payload(html)["pairs"]
        self.assertEqual(pairs[0]["taxon_id"], "T1")
        self.assertEqual(pairs[0]["taxon"], "d__Bacteria; g__Bacteroides; s__fragilis")
        self.assertEqual(pairs[0]["taxon_short"], "s__fragilis")
        self.assertEqual(pairs[2]["taxon"], "Bacteria; Klebsiella pneumoniae")
        self.assertEqual(pairs[2]["taxon_short"], "Klebsiella pneumoniae")

    def test_action_requires_tfa_property(self):
        table = qiime2.Artifact.import_data("FeatureTable[Frequency]", self._table())
        mapping = qiime2.Artifact.import_data(
            "FeatureData[Taxonomy % Properties('tfa')]", self._mapping()
        )
        with self.assertRaisesRegex(TypeError, "Properties.*tfa"):
            plugin.visualizers["explore_tfa"](feature_load=table, gene_taxonomy=mapping)

    def test_action_requires_tfa_taxonomy_property(self):
        """Require the tfa property on the gene taxonomy input."""
        table = qiime2.Artifact.import_data(
            "FeatureTable[Frequency % Properties('tfa')]", self._table()
        )
        ordinary = qiime2.Artifact.import_data("FeatureData[Taxonomy]", self._mapping())
        with self.assertRaisesRegex(TypeError, "Properties.*tfa"):
            plugin.visualizers["explore_tfa"](
                feature_load=table, gene_taxonomy=ordinary
            )

    def test_mapping_labels_and_alignment_for_filtered_table(self):
        mapping = self._mapping().iloc[::-1]
        table = self._table().filter(
            [_gene_taxonomy_id("T1", "geneA")], axis="observation", inplace=False
        )
        explore_tfa(self.temp_dir.name, table, mapping)
        html = (Path(self.temp_dir.name) / "index.html").read_text(encoding="utf-8")
        pairs = self._payload(html)["pairs"]
        self.assertEqual(len(pairs), 1)
        self.assertEqual(pairs[0]["taxon"], "Bacteria; Escherichia coli")
        self.assertEqual(pairs[0]["function"], "geneA")
        self.assertEqual(pairs[0]["total"], 8)

    def test_missing_mapping_reports_feature_id(self):
        with self.assertRaisesRegex(
            ValueError, "Missing gene taxonomy mappings.*T1\\|gene%20B"
        ):
            explore_tfa(self.temp_dir.name, self._table(), self._mapping().iloc[:1])
