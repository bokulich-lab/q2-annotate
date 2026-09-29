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

from ..tfa import (
    _align_contig_matrices,
    _gene_taxonomy_id,
    _map_contigs_to_taxa,
    estimate_tfa,
    tfa_by_taxon,
    tfa_by_function,
    tfa_collapse_samples,
)


class TestTFA(TestPluginBase):
    package = "q2_annotate.tfa.tests"

    def setUp(self):
        super().setUp()
        self.abundance = pd.read_csv(
            self.get_data_path("contig-abundance.tsv"), sep="\t", index_col=0
        ).T
        self.inventory = pd.read_csv(
            self.get_data_path("gene-inventory.tsv"), sep="\t", index_col=0
        ).T
        self.taxonomy = pd.DataFrame(
            {"Taxon": ["Escherichia coli", "Klebsiella pneumoniae"]},
            index=pd.Index(["T1", "T2"], name="Feature ID"),
        )
        self.mapping = json.load(open(self.get_data_path("taxon-to-contigs.json")))
        self.expected_tfa = pd.read_csv(
            self.get_data_path("expected-tfa-table.tsv"), sep="\t", index_col=0
        ).astype(float)
        self.expected_gene_taxonomy = pd.read_csv(
            self.get_data_path("expected-gene-taxonomy.tsv"),
            sep="\t",
            index_col=0,
            dtype=str,
            keep_default_na=False,
        )

    def _estimate(self, abundance=None, inventory=None, taxonomy=None, mapping=None):
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
        """Estimate per-sample loads while storing only observed taxon-gene pairs."""
        observed, gene_taxonomy = self._estimate()
        self.assertTrue(sp.issparse(observed.matrix_data))
        self.assertEqual(
            observed.matrix_data.nnz, np.count_nonzero(self.expected_tfa.values)
        )
        pd.testing.assert_frame_equal(
            observed.to_dataframe(dense=True).rename_axis("Feature ID"),
            self.expected_tfa,
        )
        pd.testing.assert_frame_equal(gene_taxonomy, self.expected_gene_taxonomy)
        self.assertEqual(
            list(observed.ids(axis="observation")), list(gene_taxonomy.index)
        )

    def test_no_common_contigs_rejected(self):
        """Report when no mapped contigs overlap the input tables."""
        with self.assertRaisesRegex(
            ValueError, "No taxon/gene pairs have nonzero loads"
        ):
            self._estimate(mapping={"T1": ["not-a-contig"]})

    def test_ignores_extra_contigs(self):
        """Ignore mapped contigs absent from the gene inventory."""
        abundance = self.abundance.assign(C7=[1000, 0, 0])
        taxonomy = pd.DataFrame(
            {"Taxon": ["Escherichia coli", "Klebsiella pneumoniae", "Other"]},
            index=pd.Index(["T1", "T2", "T3"]),
        )
        mapping = {**self.mapping, "T3": ["C7", "C8"]}
        observed, gene_taxonomy = self._estimate(
            abundance=abundance, taxonomy=taxonomy, mapping=mapping
        )
        pd.testing.assert_frame_equal(
            observed.to_dataframe(dense=True).rename_axis("Feature ID"),
            self.expected_tfa,
        )
        pd.testing.assert_frame_equal(gene_taxonomy, self.expected_gene_taxonomy)

    def test_zero_abundance_rejected(self):
        """Report when all pairs have zero loads and no taxonomy can be returned."""
        with self.assertRaisesRegex(
            ValueError, "No taxon/gene pairs have nonzero loads"
        ):
            self._estimate(abundance=self.abundance * 0)

    def test_omits_pairs_with_zero_load_across_samples(self):
        """Keep loaded pairs while dropping an entirely absent taxon."""
        abundance = self.abundance.copy()
        abundance[["C2", "C4", "C6"]] = 0
        observed, gene_taxonomy = self._estimate(abundance=abundance)
        self.assertEqual(
            list(observed.ids(axis="observation")), ["T1|bla_TEM", "T1|vanA"]
        )
        self.assertEqual(list(gene_taxonomy.index), ["T1|bla_TEM", "T1|vanA"])

    def test_extreme_abundances(self):
        """Preserve very small and very large abundance values."""
        abundance = self.abundance * 0.0
        abundance.loc["S1", "C1"] = 1e-9
        abundance.loc["S1", "C2"] = 1e12
        observed, gene_taxonomy = self._estimate(abundance=abundance)
        loads = {
            tuple(gene_taxonomy.loc[pair_id, ["Taxon ID", "Gene ID"]]): values
            for pair_id, values in zip(
                observed.ids(axis="observation"), observed.matrix_data.toarray()
            )
        }
        np.testing.assert_allclose(
            loads[("T1", "bla_TEM")],
            [2e-9, 0, 0],
            rtol=1e-12,
            atol=0,
        )
        np.testing.assert_allclose(
            loads[("T2", "mecA")],
            [3e12, 0, 0],
            rtol=1e-12,
            atol=0,
        )

    def test_mapping_preserves_ambiguous_characters(self):
        """Preserve punctuation in taxon and gene IDs through estimation."""
        inventory = self.inventory.rename(columns={"bla_TEM": 'a,"b'})
        taxonomy = self.taxonomy.rename(index={"T1": 't,"1'})
        mapping = {'t,"1': self.mapping["T1"], "T2": self.mapping["T2"]}
        observed, gene_taxonomy = self._estimate(
            inventory=inventory, taxonomy=taxonomy, mapping=mapping
        )
        pairs = set(
            gene_taxonomy[["Taxon ID", "Gene ID"]].itertuples(index=False, name=None)
        )
        self.assertIn(('t,"1', 'a,"b'), pairs)

    def test_ambiguous_taxon_assignment_rejected(self):
        """Reject a contig assigned to more than one taxon."""
        mapping = {**self.mapping, "T2": ["C1", "C2", "C4", "C6"]}
        with self.assertRaisesRegex(ValueError, "multiple taxa"):
            self._estimate(mapping=mapping)

    def test_contig_mapping_skips_absent_taxonomy(self):
        """Skip mapped taxa missing from the taxonomy table."""
        self.assertEqual(
            _map_contigs_to_taxa(self.taxonomy.loc[["T1"]], self.mapping),
            {"C1": "T1", "C3": "T1", "C5": "T1"},
        )

    def test_align_contig_matrices_preserves_axes_and_values(self):
        """Align contigs across matrices without changing sample or gene values."""
        abundance = self.abundance.T.iloc[::-1]
        inventory = self.inventory.T.iloc[:, ::-1]
        contigs, abundances, genes = _align_contig_matrices(
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
            {"C1": "T1", "C3": "T1", "C5": "T1"},
        )
        self.assertEqual(contigs, ["C1", "C3", "C5"])
        self.assertIsInstance(abundances, sp.csr_matrix)
        self.assertIsInstance(genes, sp.csr_matrix)
        np.testing.assert_array_equal(
            abundances.toarray(), self.abundance[contigs].T.values
        )
        np.testing.assert_array_equal(
            genes.toarray(), self.inventory.loc[contigs].values
        )

    def test_stable_feature_ids_and_mapping(self):
        """Keep feature IDs and loads stable when input axes are reordered."""
        table, mapping = self._estimate()
        reordered, reordered_mapping = self._estimate(
            inventory=self.inventory[self.inventory.columns[::-1]],
            abundance=self.abundance[self.abundance.columns[::-1]],
        )
        pd.testing.assert_frame_equal(
            mapping.sort_index(), reordered_mapping.sort_index()
        )
        self.assertEqual(list(table.ids(axis="observation")), list(mapping.index))
        self.assertEqual(set(mapping.index), set(self.expected_gene_taxonomy.index))
        pd.testing.assert_frame_equal(
            table.to_dataframe(dense=True).sort_index(),
            reordered.to_dataframe(dense=True).sort_index(),
        )

    def test_readable_ids_escape_ambiguous_components(self):
        """Build unique, readable feature IDs from taxon and gene IDs."""
        cases = pd.read_csv(self.get_data_path("pair-ids.tsv"), sep="\t", dtype=str)
        observed = [
            _gene_taxonomy_id(taxon, gene)
            for taxon, gene in zip(cases["Taxon ID"], cases["Gene ID"])
        ]
        self.assertEqual(observed, list(cases["Feature ID"]))
        self.assertEqual(len(observed), len(set(observed)))

    def test_taxon_and_function_projections(self):
        observed, gene_taxonomy = self._estimate()
        by_taxon = tfa_by_taxon(observed, gene_taxonomy, ["T1"])["T1"]
        pd.testing.assert_frame_equal(
            by_taxon.to_dataframe(dense=True),
            pd.DataFrame(
                [[200, 240, 180], [100, 120, 90]],
                index=["bla_TEM", "vanA"],
                columns=["S1", "S2", "S3"],
            ),
            check_dtype=False,
        )
        by_function = tfa_by_function(observed, gene_taxonomy, ["bla_TEM"])["bla_TEM"]
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
            tfa_by_taxon(observed, gene_taxonomy, ["unknown"])

    def test_collapse_samples_sum_mean_median(self):
        observed, gene_taxonomy = self._estimate()
        expected = np.array([[620, 310, 0], [150, 0, 450]])
        for aggregation, values in (
            ("sum", expected),
            ("mean", expected / 3),
            ("median", np.array([[200, 100, 0], [50, 0, 150]])),
        ):
            collapsed = tfa_collapse_samples(observed, gene_taxonomy, aggregation)
            self.assertTrue(sp.issparse(collapsed.matrix_data))
            self.assertEqual(list(collapsed.ids(axis="observation")), ["T1", "T2"])
            self.assertEqual(
                list(collapsed.ids(axis="sample")), ["bla_TEM", "vanA", "mecA"]
            )
            np.testing.assert_allclose(collapsed.matrix_data.toarray(), values)
        np.testing.assert_array_equal(
            tfa_collapse_samples(observed, gene_taxonomy).matrix_data.toarray(),
            expected,
        )

    def test_collapse_median_includes_implicit_zeros(self):
        abundance = self.abundance.copy()
        abundance.loc[["S2", "S3"], :] = 0
        observed, gene_taxonomy = self._estimate(abundance=abundance)
        collapsed = tfa_collapse_samples(observed, gene_taxonomy, "median")
        self.assertEqual(collapsed.matrix_data.nnz, 0)
        self.assertEqual(collapsed.shape, (2, 3))
        with self.assertRaisesRegex(ValueError, "Unsupported sample aggregation"):
            tfa_collapse_samples(observed, gene_taxonomy, "maximum")

    def test_collapse_empty_tfa(self):
        table, gene_taxonomy = self._estimate()
        empty = table.filter([], axis="observation", inplace=False)
        gene_taxonomy = gene_taxonomy.iloc[:0]
        self.assertEqual(tfa_collapse_samples(empty, gene_taxonomy).shape, (0, 0))
        self.assertEqual(
            tfa_collapse_samples(empty, gene_taxonomy, "median").shape, (0, 0)
        )

    def test_collapse_mean_requires_samples(self):
        original, gene_taxonomy = self._estimate()
        empty_samples = biom.Table(
            sp.csr_matrix((original.shape[0], 0)),
            observation_ids=original.ids(axis="observation"),
            sample_ids=[],
        )
        with self.assertRaisesRegex(ValueError, "without samples"):
            tfa_collapse_samples(empty_samples, gene_taxonomy, "mean")

    def test_multiple_selections_only_return_requested_ids(self):
        observed, gene_taxonomy = self._estimate()
        by_taxon = tfa_by_taxon(observed, gene_taxonomy, ["T2", "T1", "T2"])
        self.assertEqual(list(by_taxon), ["T2", "T1"])
        self.assertEqual(
            list(by_taxon["T2"].ids(axis="observation")), ["bla_TEM", "mecA"]
        )
        by_function = tfa_by_function(observed, gene_taxonomy, ["vanA", "mecA"])
        self.assertEqual(list(by_function), ["vanA", "mecA"])
        self.assertEqual(list(by_function["vanA"].ids(axis="observation")), ["T1"])
        with self.assertRaisesRegex(ValueError, "at least one taxon"):
            tfa_by_taxon(observed, gene_taxonomy, [])
        with self.assertRaisesRegex(ValueError, "unknown"):
            tfa_by_function(observed, gene_taxonomy, ["vanA", "unknown"])

    def test_missing_feature_mapping_rejected(self):
        table, gene_taxonomy = self._estimate()
        missing = gene_taxonomy.iloc[1:]
        with self.assertRaisesRegex(ValueError, "Missing gene taxonomy mappings"):
            tfa_by_taxon(table, missing, ["T1"])
        with self.assertRaisesRegex(ValueError, "Missing gene taxonomy mappings"):
            tfa_collapse_samples(table, missing)

    def test_mapping_rows_are_joined_by_feature_id(self):
        table, gene_taxonomy = self._estimate()
        expected = tfa_by_taxon(table, gene_taxonomy, ["T1"])["T1"]
        observed = tfa_by_taxon(table, gene_taxonomy.iloc[::-1], ["T1"])["T1"]
        np.testing.assert_array_equal(
            observed.matrix_data.toarray(), expected.matrix_data.toarray()
        )
