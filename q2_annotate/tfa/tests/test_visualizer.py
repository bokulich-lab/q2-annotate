# -----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
import json
from pathlib import Path

import biom
import pandas as pd
import qiime2
import scipy.sparse as sp
from qiime2.plugin.testing import TestPluginBase

from q2_annotate.plugin_setup import plugin
from q2_annotate.tfa.tfa import _gene_taxonomy_id, _estimate_tfa
from q2_annotate.tfa.visualizer import _visualize_tfa


class TestTFAVisualizer(TestPluginBase):
    package = "q2_annotate.tfa.tests"

    def setUp(self):
        super().setUp()
        self.expected_payload = json.loads(
            Path(self.get_data_path("expected-explorer-data.json")).read_text()
        )

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

    def _payload(self):
        return json.loads(
            (Path(self.temp_dir.name) / "data.json").read_text(encoding="utf-8")
        )

    def test_sparse_payload_and_metadata_groups(self):
        metadata = qiime2.Metadata.load(self.get_data_path("sample-metadata.tsv"))
        _visualize_tfa(
            self.temp_dir.name, self._table(), self._mapping(), metadata=metadata
        )
        output = Path(self.temp_dir.name)
        html = (output / "index.html").read_text(encoding="utf-8")
        self.expected_payload["groups"] = {
            "group": ["case", "control", "case"],
            "batch": ["B1", "B2", "B1"],
        }
        self.assertEqual(self._payload(), self.expected_payload)
        self.assertIn('aria-label="Explorer controls"', html)
        self.assertNotIn("<header>", html)
        self.assertNotIn('id="tfa-data"', html)
        self.assertNotIn("__TFA_DATA__", html)
        self.assertEqual(
            {file.name for file in output.iterdir()},
            {
                "index.html",
                "data.json",
                "explore.js",
                "plots.json",
                "style.css",
                "taxonomy.js",
            },
        )
        self.assertIn('id="level"', html)
        self.assertEqual(
            set(json.loads((output / "plots.json").read_text())),
            {"heatmap", "boxplot", "samples"},
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
        table, mapping = _estimate_tfa(
            abundance, inventory, taxonomy, {"T1": ["C1"], "T2": ["C2"]}
        )
        _visualize_tfa(self.temp_dir.name, table, mapping)
        expected = {
            "samples": ["S1", "S2"],
            "pairs": [
                {
                    "taxon": "Taxon A",
                    "taxon_short": "Taxon A",
                    "taxon_id": "T1",
                    "lineage": ["Taxon A"],
                    "function": "geneA",
                    "total": 2,
                    "mean": 1,
                    "median": 1,
                }
            ],
            "loads": [[0, 0, 2]],
            "groups": {},
        }
        self.assertEqual(self._payload(), expected)

    def test_optional_metadata_and_labels_are_separate_from_html(self):
        table = biom.Table(
            sp.csr_matrix([[2]]),
            observation_ids=[_gene_taxonomy_id("T<script>", "gene")],
            sample_ids=["S1"],
        )
        _visualize_tfa(
            self.temp_dir.name, table, self._mapping("gene-taxonomy-html.tsv")
        )
        html = (Path(self.temp_dir.name) / "index.html").read_text(encoding="utf-8")
        self.assertNotIn("T<script>", html)
        expected = {
            "samples": ["S1"],
            "pairs": [
                {
                    "taxon": "T<script>",
                    "taxon_short": "T<script>",
                    "taxon_id": "T<script>",
                    "lineage": ["T<script>"],
                    "function": "gene",
                    "total": 2,
                    "mean": 2,
                    "median": 2,
                }
            ],
            "loads": [[0, 0, 2]],
            "groups": {},
        }
        self.assertEqual(self._payload(), expected)

    def test_even_sample_medians_include_sparse_zeros(self):
        table = self._table().filter(["S1", "S2"], axis="sample", inplace=False)
        _visualize_tfa(self.temp_dir.name, table, self._mapping())
        self.expected_payload["samples"] = ["S1", "S2"]
        self.expected_payload["loads"] = self.expected_payload["loads"][:-1]
        for pair, summary in zip(
            self.expected_payload["pairs"],
            [(8, 4, 4), (4.25, 2.125, 2.125), (12, 6, 6)],
        ):
            pair.update(zip(("total", "mean", "median"), summary))
        self.assertEqual(self._payload(), self.expected_payload)

    def test_mostly_absent_pairs_have_zero_median(self):
        original = self._table()
        table = biom.Table(
            sp.hstack([original.matrix_data, sp.csr_matrix((3, 3))]),
            observation_ids=original.ids(axis="observation"),
            sample_ids=["S1", "S2", "S3", "S4", "S5", "S6"],
        )
        _visualize_tfa(self.temp_dir.name, table, self._mapping())
        self.expected_payload["samples"] = ["S1", "S2", "S3", "S4", "S5", "S6"]
        for pair, mean in zip(
            self.expected_payload["pairs"], [8 / 6, 4.25 / 6, 16 / 6]
        ):
            pair.update(mean=mean, median=0)
        self.assertEqual(self._payload(), self.expected_payload)

    def test_samples_missing_from_metadata_remain_grouped(self):
        metadata = qiime2.Metadata.load(
            self.get_data_path("sample-metadata-partial.tsv")
        )
        _visualize_tfa(
            self.temp_dir.name, self._table(), self._mapping(), metadata=metadata
        )
        self.expected_payload["groups"] = {"group": ["case", "control", "(missing)"]}
        self.assertEqual(self._payload(), self.expected_payload)

    def test_action_registration_and_execution(self):
        self.assertIn("_visualize_tfa", plugin.visualizers)
        artifact = qiime2.Artifact.import_data(
            "FeatureTable[Frequency % Properties('tfa')]", self._table()
        )
        mapping = qiime2.Artifact.import_data(
            "FeatureData[Taxonomy % Properties('tfa')]", self._mapping()
        )
        metadata = qiime2.Metadata.load(self.get_data_path("sample-metadata.tsv"))
        result = plugin.visualizers["_visualize_tfa"](
            feature_load=artifact, gene_taxonomy=mapping, metadata=metadata
        )
        self.assertIsInstance(result.visualization, qiime2.Visualization)
        without_metadata = plugin.visualizers["_visualize_tfa"](
            feature_load=artifact, gene_taxonomy=mapping
        )
        self.assertIsInstance(without_metadata.visualization, qiime2.Visualization)

    def test_mapping_supplies_full_and_short_taxonomy_labels(self):
        """Use mapping lineages for selectors and terminal names for heatmap axes."""
        mapping = self._mapping()
        mapping.loc[mapping["Taxon ID"] == "T1", "Taxon"] = (
            "d__Bacteria; g__Bacteroides; s__fragilis"
        )
        _visualize_tfa(self.temp_dir.name, self._table(), mapping)
        for pair in self.expected_payload["pairs"][:2]:
            pair.update(
                taxon="d__Bacteria; g__Bacteroides; s__fragilis",
                taxon_short="s__fragilis",
                lineage=["d__Bacteria", "g__Bacteroides", "s__fragilis"],
            )
        self.assertEqual(self._payload(), self.expected_payload)

    def test_mapping_labels_and_alignment_for_filtered_table(self):
        mapping = self._mapping().iloc[::-1]
        table = self._table().filter(
            [_gene_taxonomy_id("T1", "geneA")], axis="observation", inplace=False
        )
        _visualize_tfa(self.temp_dir.name, table, mapping)
        self.expected_payload["pairs"] = self.expected_payload["pairs"][:1]
        self.expected_payload["loads"] = [[0, 0, 5], [0, 1, 3]]
        self.assertEqual(self._payload(), self.expected_payload)

    def test_missing_mapping_reports_feature_id(self):
        with self.assertRaisesRegex(
            ValueError, "Missing gene taxonomy mappings.*T1\\|gene%20B"
        ):
            _visualize_tfa(self.temp_dir.name, self._table(), self._mapping().iloc[:1])
