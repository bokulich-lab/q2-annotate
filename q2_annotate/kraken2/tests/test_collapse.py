# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import tempfile
import shutil
import json
from unittest.mock import patch, MagicMock, ANY

import biom
import numpy as np
import pandas as pd
from qiime2 import Artifact
from qiime2.plugin.testing import TestPluginBase
from q2_types.kraken2 import (
    Kraken2OutputDirectoryFormat,
    Kraken2ReportDirectoryFormat,
)

from q2_annotate.kraken2.collapse import (
    _build_contig_map,
    _average_by_count,
    _df_to_json_per_sample,
    _table_to_json,
    _extract_mean_abundances,
    _visualize_collapsed_contigs,
    collapse_contigs,
    map_taxonomy_to_contigs,
)


class TestBuildContigMap(TestPluginBase):
    package = "q2_annotate.kraken2.tests"

    def test_build_contig_map_single_output(self):
        """Test building contig map from a single output file."""
        outputs_dir = self.get_data_path("collapse/single-output")
        outputs_format = Kraken2OutputDirectoryFormat(outputs_dir, mode="r")

        obs = _build_contig_map(outputs_format)

        exp = {
            "taxon1": ["contig1", "contig2"],
            "taxon2": ["contig3"],
        }

        self.assertEqual(obs, exp)

    def test_build_contig_map_multiple_outputs(self):
        """Test building contig map from multiple output files."""
        outputs_dir = self.get_data_path("collapse/multiple-outputs")
        outputs_format = Kraken2OutputDirectoryFormat(outputs_dir, mode="r")

        obs = _build_contig_map(outputs_format)

        exp = {
            "taxon1": ["contig1", "contig3"],
            "taxon2": ["contig2", "contig4"],
        }

        self.assertEqual(obs, exp)

    def test_build_contig_map_unclassified_only(self):
        """Test building contig map with only unclassified sequences."""
        outputs_dir = self.get_data_path("collapse/unclassified-only")
        outputs_format = Kraken2OutputDirectoryFormat(outputs_dir, mode="r")

        obs = _build_contig_map(outputs_format)

        # Unclassified sequences map to taxon "0"
        self.assertEqual(obs, {"0": ["contig1"]})


class TestAverageByCount(TestPluginBase):
    package = "q2_annotate.kraken2.tests"

    def setUp(self):
        super().setUp()
        # Create a collapsed biom table (taxa x samples)
        collapsed_data = np.array([[30.0, 50.0], [30.0, 40.0]])
        collapsed_obs_ids = ["taxon1", "taxon2"]
        sample_ids = ["sample1", "sample2"]
        self.collapsed_table = biom.Table(collapsed_data, collapsed_obs_ids, sample_ids)

        # Create original biom table (contigs x samples)
        # taxon1 has contig1 (10) and contig2 (20) in sample1, contig1 (30) and contig2 (20) in sample2
        # taxon2 has contig3 (30) in sample1, contig3 (40) in sample2
        original_data = np.array([[10.0, 30.0], [20.0, 20.0], [30.0, 40.0]])
        original_obs_ids = ["contig1", "contig2", "contig3"]
        self.original_table = biom.Table(original_data, original_obs_ids, sample_ids)

        # Create reverse contig map (contig -> taxon)
        self.contig_map_rev = {
            "contig1": "taxon1",
            "contig2": "taxon1",
            "contig3": "taxon2",
        }

    def test_average_by_count(self):
        """Test averaging abundances by contig count per sample."""
        obs = _average_by_count(self.collapsed_table, self.original_table, self.contig_map_rev)

        # Expected: divide by count of contigs per taxon per sample
        # taxon1 in sample1: 30 / 2 = 15.0 (contig1=10, contig2=20)
        # taxon1 in sample2: 50 / 2 = 25.0 (contig1=30, contig2=20)
        # taxon2 in sample1: 30 / 1 = 30.0 (contig3=30)
        # taxon2 in sample2: 40 / 1 = 40.0 (contig3=40)
        exp_data = np.array([[15.0, 25.0], [30.0, 40.0]])
        exp_table = biom.Table(exp_data, ["taxon1", "taxon2"], ["sample1", "sample2"])

        obs_df = obs.to_dataframe(dense=True)
        exp_df = exp_table.to_dataframe(dense=True)

        pd.testing.assert_frame_equal(obs_df, exp_df)

    def test_average_by_count_single_contig(self):
        """Test averaging when each taxon has one contig per sample."""
        # Create collapsed table
        collapsed_data = np.array([[10.0, 20.0], [30.0, 40.0]])
        collapsed_table = biom.Table(collapsed_data, ["taxon1", "taxon2"], ["sample1", "sample2"])

        # Create original table where each taxon has one contig
        original_data = np.array([[10.0, 20.0], [30.0, 40.0]])
        original_table = biom.Table(original_data, ["contig1", "contig2"], ["sample1", "sample2"])

        contig_map_rev = {
            "contig1": "taxon1",
            "contig2": "taxon2",
        }

        obs = _average_by_count(collapsed_table, original_table, contig_map_rev)

        # Should be unchanged since dividing by 1
        obs_df = obs.to_dataframe(dense=True)
        exp_df = collapsed_table.to_dataframe(dense=True)

        pd.testing.assert_frame_equal(obs_df, exp_df)

    def test_average_by_count_per_sample_variation(self):
        """Test that averaging uses per-sample contig counts, not global counts."""
        # Create collapsed table
        collapsed_data = np.array([[30.0, 10.0]])
        collapsed_table = biom.Table(collapsed_data, ["taxon1"], ["sample1", "sample2"])

        # Create original table where taxon1 has 2 contigs in sample1, 1 in sample2
        # sample1: contig1=10, contig2=20 (both present)
        # sample2: contig1=10, contig2=0 (only contig1 present)
        original_data = np.array([[10.0, 10.0], [20.0, 0.0]])
        original_table = biom.Table(original_data, ["contig1", "contig2"], ["sample1", "sample2"])

        contig_map_rev = {
            "contig1": "taxon1",
            "contig2": "taxon1",
        }

        obs = _average_by_count(collapsed_table, original_table, contig_map_rev)

        # Expected:
        # sample1: 30 / 2 = 15.0 (2 contigs with non-zero values)
        # sample2: 10 / 1 = 10.0 (1 contig with non-zero value)
        exp_data = np.array([[15.0, 10.0]])
        exp_table = biom.Table(exp_data, ["taxon1"], ["sample1", "sample2"])

        obs_df = obs.to_dataframe(dense=True)
        exp_df = exp_table.to_dataframe(dense=True)

        pd.testing.assert_frame_equal(obs_df, exp_df)


class TestDfToJsonPerSample(TestPluginBase):
    package = "q2_annotate.kraken2.tests"

    def setUp(self):
        super().setUp()
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_df_to_json_per_sample(self):
        """Test writing DataFrame to separate JSON files per sample."""
        # Create DataFrame with array column
        df = pd.DataFrame(
            {
                "taxon": ["taxon1", "taxon2", "taxon1"],
                "sample": ["sample1", "sample1", "sample2"],
                "abundances": [[1.0, 2.0, 3.0], [4.0, 5.0], [6.0]],
            }
        )

        # Create data directory
        os.makedirs(os.path.join(self.temp_dir, "data"), exist_ok=True)

        _df_to_json_per_sample(df, self.temp_dir)

        # Verify JSON files were created
        sample1_file = os.path.join(self.temp_dir, "data", "sample1.json")
        sample2_file = os.path.join(self.temp_dir, "data", "sample2.json")
        self.assertTrue(os.path.exists(sample1_file))
        self.assertTrue(os.path.exists(sample2_file))

        # Read back and verify sample1
        with open(sample1_file, "r") as f:
            sample1_data = json.load(f)

        # Sample1 should have exploded abundances for taxon1 and taxon2
        expected_sample1 = [
            {"taxon": "taxon1", "abundance": 1.0},
            {"taxon": "taxon1", "abundance": 2.0},
            {"taxon": "taxon1", "abundance": 3.0},
            {"taxon": "taxon2", "abundance": 4.0},
            {"taxon": "taxon2", "abundance": 5.0},
        ]
        self.assertEqual(len(sample1_data), 5)
        # Sort for comparison
        sample1_data_sorted = sorted(
            sample1_data, key=lambda x: (x["taxon"], x["abundance"])
        )
        expected_sample1_sorted = sorted(
            expected_sample1, key=lambda x: (x["taxon"], x["abundance"])
        )
        self.assertEqual(sample1_data_sorted, expected_sample1_sorted)

        # Read back and verify sample2
        with open(sample2_file, "r") as f:
            sample2_data = json.load(f)

        expected_sample2 = [
            {"taxon": "taxon1", "abundance": 6.0},
        ]
        self.assertEqual(sample2_data, expected_sample2)


class TestTableToJson(TestPluginBase):
    package = "q2_annotate.kraken2.tests"

    def setUp(self):
        super().setUp()
        self.temp_dir = tempfile.mkdtemp()

        # Create a sparse biom table
        # Each contig is unique to one sample:
        # Contig1 -> taxon1 in sample1, Contig2 -> taxon1 in sample1,
        # Contig3 -> taxon2 in sample2
        data = np.array([[10.0, 0.0], [20.0, 0.0], [0.0, 30.0]])
        obs_ids = ["contig1", "contig2", "contig3"]
        sample_ids = ["sample1", "sample2"]

        self.table = biom.Table(data, obs_ids, sample_ids)

        # Reverse contig map: contig -> taxon
        self.contig_map_rev = {
            "contig1": "taxon1",
            "contig2": "taxon1",
            "contig3": "taxon2",
        }

        # Taxonomy mapping
        self.taxonomy = pd.Series(
            {
                "taxon1": (
                    "d__Bacteria;p__Firmicutes;c__Bacilli;"
                    "g__Bacillus;s__Bacillus_subtilis"
                ),
                "taxon2": (
                    "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;"
                    "g__Escherichia;s__Escherichia_coli"
                ),
            }
        )

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_table_to_json_with_taxonomy(self):
        """Test saving table data as JSON files with taxonomy."""
        # Create data directory (normally done by _visualize_collapsed_contigs)
        os.makedirs(os.path.join(self.temp_dir, "data"), exist_ok=True)

        _table_to_json(self.table, self.contig_map_rev, self.taxonomy, self.temp_dir)

        # Verify JSON files were created
        sample1_file = os.path.join(self.temp_dir, "data", "sample1.json")
        sample2_file = os.path.join(self.temp_dir, "data", "sample2.json")
        self.assertTrue(os.path.exists(sample1_file))
        self.assertTrue(os.path.exists(sample2_file))

        # Read back and verify sample1
        with open(sample1_file, "r") as f:
            sample1_data = json.load(f)

        # Sample1 should have taxon1 with abundances [10.0, 20.0]
        expected_sample1 = [
            {"taxon": self.taxonomy["taxon1"], "abundance": 10.0},
            {"taxon": self.taxonomy["taxon1"], "abundance": 20.0},
        ]
        sample1_sorted = sorted(sample1_data, key=lambda x: x["abundance"])
        expected_sample1_sorted = sorted(expected_sample1, key=lambda x: x["abundance"])
        self.assertEqual(sample1_sorted, expected_sample1_sorted)

        # Read back and verify sample2
        with open(sample2_file, "r") as f:
            sample2_data = json.load(f)

        expected_sample2 = [
            {"taxon": self.taxonomy["taxon2"], "abundance": 30.0},
        ]
        self.assertEqual(sample2_data, expected_sample2)

    def test_table_to_json_without_taxonomy(self):
        """Test saving table data as JSON files without taxonomy."""
        # Create data directory (normally done by _visualize_collapsed_contigs)
        os.makedirs(os.path.join(self.temp_dir, "data"), exist_ok=True)

        _table_to_json(self.table, self.contig_map_rev, None, self.temp_dir)

        # Verify JSON files were created
        sample1_file = os.path.join(self.temp_dir, "data", "sample1.json")
        sample2_file = os.path.join(self.temp_dir, "data", "sample2.json")
        self.assertTrue(os.path.exists(sample1_file))
        self.assertTrue(os.path.exists(sample2_file))

        # Read back and verify sample1 (using taxon IDs)
        with open(sample1_file, "r") as f:
            sample1_data = json.load(f)

        expected_sample1 = [
            {"taxon": "taxon1", "abundance": 10.0},
            {"taxon": "taxon1", "abundance": 20.0},
        ]
        sample1_sorted = sorted(sample1_data, key=lambda x: x["abundance"])
        expected_sample1_sorted = sorted(expected_sample1, key=lambda x: x["abundance"])
        self.assertEqual(sample1_sorted, expected_sample1_sorted)

        # Read back and verify sample2
        with open(sample2_file, "r") as f:
            sample2_data = json.load(f)

        expected_sample2 = [
            {"taxon": "taxon2", "abundance": 30.0},
        ]
        self.assertEqual(sample2_data, expected_sample2)

    def test_table_to_json_missing_contig(self):
        """Test saving when contig is not in map (should use '0')."""
        # Create data directory (normally done by _visualize_collapsed_contigs)
        os.makedirs(os.path.join(self.temp_dir, "data"), exist_ok=True)

        contig_map_rev = {
            "contig1": "taxon1",
            # contig2 missing - will map to "0"
            "contig3": "taxon2",
        }

        _table_to_json(self.table, contig_map_rev, None, self.temp_dir)

        # Verify JSON files were created
        sample1_file = os.path.join(self.temp_dir, "data", "sample1.json")
        sample2_file = os.path.join(self.temp_dir, "data", "sample2.json")
        self.assertTrue(os.path.exists(sample1_file))
        self.assertTrue(os.path.exists(sample2_file))

        # Read back and verify sample1 (contig2 maps to "0")
        with open(sample1_file, "r") as f:
            sample1_data = json.load(f)

        # Sample1 should have taxon1 (contig1=10.0) and "0" (contig2=20.0)
        expected_sample1 = [
            {"taxon": "taxon1", "abundance": 10.0},
            {"taxon": "0", "abundance": 20.0},
        ]
        sample1_sorted = sorted(
            sample1_data, key=lambda x: (x["taxon"], x["abundance"])
        )
        expected_sample1_sorted = sorted(
            expected_sample1, key=lambda x: (x["taxon"], x["abundance"])
        )
        self.assertEqual(sample1_sorted, expected_sample1_sorted)

        # Read back and verify sample2
        with open(sample2_file, "r") as f:
            sample2_data = json.load(f)

        expected_sample2 = [
            {"taxon": "taxon2", "abundance": 30.0},
        ]
        self.assertEqual(sample2_data, expected_sample2)


class TestExtractMeanAbundances(TestPluginBase):
    package = "q2_annotate.kraken2.tests"

    def test_extract_mean_abundances_with_taxonomy(self):
        """Test extracting mean abundances with taxonomy mapping."""
        # Create collapsed table
        data = np.array([[15.0, 20.0], [30.0, 0.0]])
        obs_ids = ["taxon1", "taxon2"]
        sample_ids = ["sample1", "sample2"]
        collapsed_table = biom.Table(data, obs_ids, sample_ids)

        # Create taxonomy mapping
        taxonomy = pd.Series(
            {
                "taxon1": (
                    "d__Bacteria;p__Firmicutes;c__Bacilli;"
                    "g__Bacillus;s__Bacillus_subtilis"
                ),
                "taxon2": (
                    "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;"
                    "g__Escherichia;s__Escherichia_coli"
                ),
            }
        )

        result = _extract_mean_abundances(collapsed_table, taxonomy)

        expected = {
            taxonomy["taxon1"]: {"sample1": 15.0, "sample2": 20.0},
            taxonomy["taxon2"]: {"sample1": 30.0},
        }
        self.assertEqual(result, expected)

    def test_extract_mean_abundances_without_taxonomy(self):
        """Test extracting mean abundances without taxonomy (uses taxon IDs)."""
        # Create collapsed table
        data = np.array([[15.0], [30.0]])
        obs_ids = ["taxon1", "taxon2"]
        sample_ids = ["sample1"]
        collapsed_table = biom.Table(data, obs_ids, sample_ids)

        result = _extract_mean_abundances(collapsed_table, None)

        expected = {
            "taxon1": {"sample1": 15.0},
            "taxon2": {"sample1": 30.0},
        }
        self.assertEqual(result, expected)

    def test_extract_mean_abundances_excludes_zero_values(self):
        """Test that zero values are excluded from mean abundances."""
        # Create collapsed table with zero values
        data = np.array([[10.0, 0.0], [0.0, 20.0], [5.0, 5.0]])
        obs_ids = ["taxon1", "taxon2", "taxon3"]
        sample_ids = ["sample1", "sample2"]
        collapsed_table = biom.Table(data, obs_ids, sample_ids)

        result = _extract_mean_abundances(collapsed_table, None)

        expected = {
            "taxon1": {"sample1": 10.0},
            "taxon2": {"sample2": 20.0},
            "taxon3": {"sample1": 5.0, "sample2": 5.0},
        }
        self.assertEqual(result, expected)

    def test_extract_mean_abundances_with_missing_taxonomy_mapping(self):
        """Test that missing taxonomy mappings fall back to taxon ID."""
        # Create collapsed table
        data = np.array([[15.0]])
        obs_ids = ["taxon1"]
        sample_ids = ["sample1"]
        collapsed_table = biom.Table(data, obs_ids, sample_ids)

        # Create taxonomy that doesn't include taxon1
        taxonomy = pd.Series(
            {
                "taxon2": "d__Bacteria;p__Firmicutes",
            }
        )

        result = _extract_mean_abundances(collapsed_table, taxonomy)

        # Should fall back to taxon ID since taxon1 not in taxonomy
        expected = {"taxon1": {"sample1": 15.0}}
        self.assertEqual(result, expected)

    def test_extract_mean_abundances_empty_table(self):
        """Test extracting mean abundances from empty table."""
        # Create empty collapsed table
        data = np.array([]).reshape(0, 0)
        obs_ids = []
        sample_ids = []
        collapsed_table = biom.Table(data, obs_ids, sample_ids)

        result = _extract_mean_abundances(collapsed_table, None)

        # Should return empty dict
        self.assertEqual(result, {})


class TestVisualizeCollapsedContigs(TestPluginBase):
    package = "q2_annotate.kraken2.tests"

    def setUp(self):
        super().setUp()
        self.temp_dir = tempfile.mkdtemp()

        # Create a simple biom table
        data = np.array([[10.0], [20.0]])
        obs_ids = ["contig1", "contig2"]
        sample_ids = ["sample1"]

        self.table = biom.Table(data, obs_ids, sample_ids)

        self.contig_map = {
            "taxon1": ["contig1", "contig2"],
        }

        # Create collapsed table (averaged abundances per taxon)
        # taxon1 has contig1 and contig2, so average = (10.0 + 20.0) / 2 = 15.0
        collapsed_data = np.array([[15.0]])
        collapsed_obs_ids = ["taxon1"]
        self.collapsed_table = biom.Table(collapsed_data, collapsed_obs_ids, sample_ids)

        taxon1_str = (
            "d__Bacteria;p__Firmicutes;c__Bacilli;" "g__Bacillus;s__Bacillus_subtilis"
        )
        self.taxonomy = pd.Series(
            {
                "taxon1": taxon1_str,
            }
        )

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    @patch("q2_annotate.kraken2.collapse.q2templates.render")
    @patch("q2_annotate.kraken2.collapse.shutil.copytree")
    def test_visualize_collapsed_contigs_with_taxonomy(
        self, mock_copytree, mock_render
    ):
        """Test visualization generation with taxonomy."""
        _visualize_collapsed_contigs(
            self.temp_dir,
            self.table,
            self.collapsed_table,
            self.contig_map,
            self.taxonomy,
        )

        # Verify data directory was created
        self.assertTrue(os.path.exists(os.path.join(self.temp_dir, "data")))

        # Verify JSON file was created
        json_file = os.path.join(self.temp_dir, "data", "sample1.json")
        self.assertTrue(os.path.exists(json_file))

        # Verify JS/CSS files were copied
        self.assertEqual(mock_copytree.call_count, 2)

        # Verify render was called with correct context
        mock_render.assert_called_once()
        call_args = mock_render.call_args
        context = call_args[1]["context"]

        # Build expected context
        expected_context = {
            "samples": json.dumps(["sample1"]),
            "mean_abundances": json.dumps({self.taxonomy["taxon1"]: {"sample1": 15.0}}),
        }

        # Verify samples content
        context_no_vega = context.copy()
        del context_no_vega["vega_abundance_histogram_spec"]
        self.assertDictEqual(context_no_vega, expected_context)

        # Verify Vega spec is a valid dict with expected structure
        vega_spec = json.loads(context["vega_abundance_histogram_spec"])
        self.assertIsInstance(vega_spec, dict)
        self.assertIn("$schema", vega_spec)
        self.assertIn("data", vega_spec)
        self.assertEqual(vega_spec["data"]["name"], "source")

    @patch("q2_annotate.kraken2.collapse.q2templates.render")
    @patch("q2_annotate.kraken2.collapse.shutil.copytree")
    def test_visualize_collapsed_contigs_without_taxonomy(
        self, mock_copytree, mock_render
    ):
        """Test visualization generation without taxonomy."""
        _visualize_collapsed_contigs(
            self.temp_dir, self.table, self.collapsed_table, self.contig_map, None
        )

        # Verify data directory was created
        self.assertTrue(os.path.exists(os.path.join(self.temp_dir, "data")))

        # Verify JSON file was created
        json_file = os.path.join(self.temp_dir, "data", "sample1.json")
        self.assertTrue(os.path.exists(json_file))

        # Verify render was called with correct context
        mock_render.assert_called_once()
        call_args = mock_render.call_args
        context = call_args[1]["context"]

        # Build expected context (without taxonomy, uses taxon IDs)
        expected_context = {
            "samples": json.dumps(["sample1"]),
            "mean_abundances": json.dumps({"taxon1": {"sample1": 15.0}}),
        }

        # Verify samples content
        context_no_vega = context.copy()
        del context_no_vega["vega_abundance_histogram_spec"]
        self.assertDictEqual(context_no_vega, expected_context)

        # Verify Vega spec is a valid dict with expected structure
        vega_spec = json.loads(context["vega_abundance_histogram_spec"])
        self.assertIsInstance(vega_spec, dict)
        self.assertIn("$schema", vega_spec)
        self.assertIn("data", vega_spec)
        self.assertEqual(vega_spec["data"]["name"], "source")


class TestCollapseContigs(TestPluginBase):
    package = "q2_annotate.kraken2.tests"

    def setUp(self):
        super().setUp()

        # Create a simple biom table
        data = np.array([[10.0], [20.0], [30.0]])
        obs_ids = ["contig1", "contig2", "contig3"]
        sample_ids = ["sample1"]

        biom_table = biom.Table(data, obs_ids, sample_ids)

        contig_map = {
            "taxon1": ["contig1", "contig2"],
            "taxon2": ["contig3"],
        }

        # Create real QIIME 2 artifacts
        self.table_artifact = Artifact.import_data(
            "FeatureTable[Frequency]", biom_table
        )
        self.contig_map_artifact = Artifact.import_data(
            "FeatureMap[TaxonomyToContigs]", contig_map
        )

        # Create taxonomy artifact
        taxonomy = pd.Series(
            {
                "taxon1": (
                    "d__Bacteria;p__Firmicutes;c__Bacilli;"
                    "g__Bacillus;s__Bacillus_subtilis"
                ),
                "taxon2": (
                    "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;"
                    "g__Escherichia;s__Escherichia_coli"
                ),
            },
            name="Taxon",
        )
        taxonomy.index.name = "Feature ID"
        self.taxonomy_artifact = Artifact.import_data("FeatureData[Taxonomy]", taxonomy)

    @patch("q2_annotate.kraken2.collapse._visualize_collapsed_contigs")
    def test_collapse_contigs_basic(self, mock_visualize):
        """Test basic contig collapsing."""
        mock_ctx = MagicMock()
        mock_visualize_action = MagicMock(return_value=(MagicMock(),))
        mock_ctx.get_action.return_value = mock_visualize_action

        mock_artifact = MagicMock()
        mock_ctx.make_artifact.return_value = mock_artifact

        result_table, result_viz = collapse_contigs(
            mock_ctx, self.contig_map_artifact, self.table_artifact, None
        )

        # Verify get_action was called with correct arguments
        mock_ctx.get_action.assert_called_once_with(
            "annotate", "_visualize_collapsed_contigs"
        )

        # Verify visualization action was called with correct artifacts
        mock_visualize_action.assert_called_once()
        call_args = mock_visualize_action.call_args[0]
        self.assertEqual(call_args[0], self.table_artifact)
        # Second argument should be the collapsed_table artifact
        # (created by make_artifact)
        # Since make_artifact returns a mock, we just verify it's the mock
        self.assertEqual(call_args[1], mock_artifact)
        self.assertEqual(call_args[2], self.contig_map_artifact)
        self.assertIsNone(call_args[3])

        # Verify make_artifact was called with correct arguments
        mock_ctx.make_artifact.assert_called_once()
        make_artifact_call_args = mock_ctx.make_artifact.call_args[0]
        self.assertEqual(make_artifact_call_args[0], "FeatureTable[Frequency]")
        # Verify the second argument is a biom.Table
        self.assertIsInstance(make_artifact_call_args[1], biom.Table)

        # Verify the collapsed table has correct structure
        collapsed_table = make_artifact_call_args[1]
        collapsed_df = collapsed_table.to_dataframe(dense=True)

        # Build expected DataFrame
        expected_data = np.array([[15.0], [30.0]])
        expected_df = pd.DataFrame(
            expected_data,
            index=["taxon1", "taxon2"],
            columns=["sample1"],
        )

        pd.testing.assert_frame_equal(collapsed_df, expected_df)

        self.assertEqual(result_table, mock_artifact)

    @patch("q2_annotate.kraken2.collapse._visualize_collapsed_contigs")
    def test_collapse_contigs_with_taxonomy(self, mock_visualize):
        """Test contig collapsing with taxonomy."""
        mock_ctx = MagicMock()
        mock_visualize_action = MagicMock(return_value=(MagicMock(),))
        mock_ctx.get_action.return_value = mock_visualize_action

        mock_artifact = MagicMock()
        mock_ctx.make_artifact.return_value = mock_artifact

        result_table, result_viz = collapse_contigs(
            mock_ctx,
            self.contig_map_artifact,
            self.table_artifact,
            self.taxonomy_artifact,
        )

        # Verify get_action was called with correct arguments
        mock_ctx.get_action.assert_called_once_with(
            "annotate", "_visualize_collapsed_contigs"
        )

        # Verify visualization action was called with correct artifacts
        # including taxonomy
        mock_visualize_action.assert_called_once()
        call_args = mock_visualize_action.call_args[0]
        self.assertEqual(call_args[0], self.table_artifact)
        # Second argument should be the collapsed_table artifact
        # (created by make_artifact)
        # Since make_artifact returns a mock, we just verify it's the mock
        self.assertEqual(call_args[1], mock_artifact)
        self.assertEqual(call_args[2], self.contig_map_artifact)
        self.assertEqual(call_args[3], self.taxonomy_artifact)

        # Verify make_artifact was called with correct arguments
        mock_ctx.make_artifact.assert_called_once()
        make_artifact_call_args = mock_ctx.make_artifact.call_args[0]
        self.assertEqual(make_artifact_call_args[0], "FeatureTable[Frequency]")
        # Verify the second argument is a biom.Table
        self.assertIsInstance(make_artifact_call_args[1], biom.Table)

        # Verify the collapsed table has correct structure
        collapsed_table = make_artifact_call_args[1]
        collapsed_df = collapsed_table.to_dataframe(dense=True)

        # Build expected DataFrame
        expected_data = np.array([[15.0], [30.0]])
        expected_df = pd.DataFrame(
            expected_data,
            index=["taxon1", "taxon2"],
            columns=["sample1"],
        )

        pd.testing.assert_frame_equal(collapsed_df, expected_df)

        self.assertEqual(result_table, mock_artifact)


class TestMapTaxonomyToContigs(TestPluginBase):
    package = "q2_annotate.kraken2.tests"

    def setUp(self):
        super().setUp()
        self.reports = Kraken2ReportDirectoryFormat(
            self.get_data_path("reports-mags"), "r"
        )
        self.outputs = Kraken2OutputDirectoryFormat(
            self.get_data_path("outputs-mags"), "r"
        )

        self.reports_artifact = Artifact.import_data(
            'SampleData[Kraken2Report % Properties("contigs")]', self.reports
        )
        self.outputs_artifact = Artifact.import_data(
            'SampleData[Kraken2Output % Properties("contigs")]', self.outputs
        )

    @patch("q2_annotate.kraken2.collapse._build_contig_map")
    def test_map_taxonomy_to_contigs_basic(self, mock_build_map):
        """Test mapping taxonomy to contigs."""
        mock_ctx = MagicMock()

        # Create real taxonomy artifact for kraken2_to_features return value
        taxonomy = pd.Series(
            {
                "1912795": "d__Bacteria;p__Firmicutes",
                "1583098": "d__Bacteria;p__Proteobacteria",
            },
            name="Taxon",
        )
        taxonomy.index.name = "Feature ID"
        taxonomy_artifact = Artifact.import_data("FeatureData[Taxonomy]", taxonomy)

        # Mock kraken2_to_features action
        mock_to_features_action = MagicMock()
        mock_to_features_action.return_value = (None, taxonomy_artifact)
        mock_ctx.get_action.return_value = mock_to_features_action

        # Mock contig map
        mock_contig_map = {
            "1912795": ["contig1", "contig2"],
            "1583098": ["contig3"],
        }
        mock_build_map.return_value = mock_contig_map

        # Create real artifacts for return values
        taxonomy_result = pd.Series(
            {
                "0": "d__Unclassified",
                "1912795": "d__Bacteria;p__Firmicutes",
                "1583098": "d__Bacteria;p__Proteobacteria",
            },
            name="Taxon",
        )
        taxonomy_result.index.name = "Feature ID"
        taxonomy_result_artifact = Artifact.import_data(
            "FeatureData[Taxonomy]", taxonomy_result
        )

        map_result_dict = {
            "1912795": ["contig1", "contig2"],
            "1583098": ["contig3"],
        }
        map_result_artifact = Artifact.import_data(
            "FeatureMap[TaxonomyToContigs]", map_result_dict
        )

        mock_ctx.make_artifact.side_effect = [
            taxonomy_result_artifact,
            map_result_artifact,
        ]

        result_map, result_taxonomy = map_taxonomy_to_contigs(
            mock_ctx,
            self.reports_artifact,
            self.outputs_artifact,
            coverage_threshold=0.1,
        )

        # Verify get_action was called with correct arguments
        mock_ctx.get_action.assert_called_once_with("annotate", "kraken2_to_features")

        # Verify kraken2_to_features action was called with correct arguments
        mock_to_features_action.assert_called_once()
        to_features_call_args = mock_to_features_action.call_args[0]
        self.assertEqual(to_features_call_args[0], self.reports_artifact)
        self.assertEqual(to_features_call_args[1], 0.1)

        # Verify _build_contig_map was called with correct arguments
        mock_build_map.assert_called_once()
        build_map_call_args = mock_build_map.call_args[0]
        self.assertEqual(
            build_map_call_args[0],
            ANY,
        )

        # Verify artifacts were created with correct arguments
        self.assertEqual(mock_ctx.make_artifact.call_count, 2)

        # Verify taxonomy artifact creation
        taxonomy_call_args = mock_ctx.make_artifact.call_args_list[0]
        self.assertEqual(taxonomy_call_args[0][0], "FeatureData[Taxonomy]")
        taxonomy = taxonomy_call_args[0][1]
        self.assertIsInstance(taxonomy, pd.Series)
        self.assertIn("0", taxonomy.index)
        self.assertEqual(taxonomy["0"], "d__Unclassified")

        # Verify map artifact creation
        map_call_args = mock_ctx.make_artifact.call_args_list[1]
        self.assertEqual(map_call_args[0][0], "FeatureMap[TaxonomyToContigs]")
        self.assertIsInstance(map_call_args[0][1], dict)

    @patch("q2_annotate.kraken2.collapse._build_contig_map")
    def test_map_taxonomy_to_contigs_adds_unclassified(self, mock_build_map):
        """Test that unclassified taxon is added."""
        mock_ctx = MagicMock()

        # Create real taxonomy artifact for kraken2_to_features return value
        taxonomy = pd.Series(
            {
                "1912795": "d__Bacteria;p__Firmicutes",
            },
            name="Taxon",
        )
        taxonomy.index.name = "Feature ID"
        taxonomy_artifact = Artifact.import_data("FeatureData[Taxonomy]", taxonomy)

        mock_to_features_action = MagicMock()
        mock_to_features_action.return_value = (None, taxonomy_artifact)
        mock_ctx.get_action.return_value = mock_to_features_action

        mock_contig_map = {
            "1912795": ["contig1"],
            "999": ["contig2"],  # Missing from taxonomy
        }
        mock_build_map.return_value = mock_contig_map

        # Create real artifacts for return values
        taxonomy_result = pd.Series(
            {
                "0": "d__Unclassified",
                "1912795": "d__Bacteria;p__Firmicutes",
            },
            name="Taxon",
        )
        taxonomy_result.index.name = "Feature ID"
        taxonomy_result_artifact = Artifact.import_data(
            "FeatureData[Taxonomy]", taxonomy_result
        )

        map_result = {
            "1912795": ["contig1"],
            "0": ["contig2"],  # Missing taxon moved to "0"
        }
        map_result_artifact = Artifact.import_data(
            "FeatureMap[TaxonomyToContigs]", map_result
        )

        mock_ctx.make_artifact.side_effect = [
            taxonomy_result_artifact,
            map_result_artifact,
        ]

        result_map, result_taxonomy = map_taxonomy_to_contigs(
            mock_ctx,
            self.reports_artifact,
            self.outputs_artifact,
            coverage_threshold=0.1,
        )

        # Verify get_action was called with correct arguments
        mock_ctx.get_action.assert_called_once_with("annotate", "kraken2_to_features")

        # Verify kraken2_to_features action was called with correct arguments
        mock_to_features_action.assert_called_once()
        to_features_call_args = mock_to_features_action.call_args[0]
        self.assertEqual(to_features_call_args[0], self.reports_artifact)
        self.assertEqual(to_features_call_args[1], 0.1)

        # Verify _build_contig_map was called with correct arguments
        mock_build_map.assert_called_once()
        build_map_call_args = mock_build_map.call_args[0]
        self.assertEqual(
            build_map_call_args[0],
            ANY,
        )

        # Verify unclassified was added to taxonomy
        taxonomy_call_args = mock_ctx.make_artifact.call_args_list[0]
        self.assertEqual(taxonomy_call_args[0][0], "FeatureData[Taxonomy]")
        taxonomy = taxonomy_call_args[0][1]
        self.assertIsInstance(taxonomy, pd.Series)
        self.assertIn("0", taxonomy.index)
        self.assertEqual(taxonomy["0"], "d__Unclassified")

    @patch("q2_annotate.kraken2.collapse._build_contig_map")
    def test_map_taxonomy_to_contigs_handles_missing_taxa(self, mock_build_map):
        """Test handling of taxa missing from taxonomy (below threshold)."""
        mock_ctx = MagicMock()

        # Create real taxonomy artifact for kraken2_to_features return value
        taxonomy = pd.Series(
            {
                "1912795": "d__Bacteria;p__Firmicutes",
            },
            name="Taxon",
        )
        taxonomy.index.name = "Feature ID"
        taxonomy_artifact = Artifact.import_data("FeatureData[Taxonomy]", taxonomy)

        mock_to_features_action = MagicMock()
        mock_to_features_action.return_value = (None, taxonomy_artifact)
        mock_ctx.get_action.return_value = mock_to_features_action

        # Contig map has taxa not in taxonomy (below threshold)
        # Note: "0" must exist in the map for the function to work
        mock_contig_map = {
            "0": [],  # Unclassified - must exist
            "1912795": ["contig1"],
            "999": ["contig2"],  # Missing from taxonomy
        }
        mock_build_map.return_value = mock_contig_map

        # Create real artifacts for return values
        taxonomy_result = pd.Series(
            {
                "0": "d__Unclassified",
                "1912795": "d__Bacteria;p__Firmicutes",
            },
            name="Taxon",
        )
        taxonomy_result.index.name = "Feature ID"
        taxonomy_result_artifact = Artifact.import_data(
            "FeatureData[Taxonomy]", taxonomy_result
        )

        map_result = {
            "0": ["contig2"],  # Missing taxon moved to "0"
            "1912795": ["contig1"],
        }
        map_result_artifact = Artifact.import_data(
            "FeatureMap[TaxonomyToContigs]", map_result
        )

        mock_ctx.make_artifact.side_effect = [
            taxonomy_result_artifact,
            map_result_artifact,
        ]

        result_map, result_taxonomy = map_taxonomy_to_contigs(
            mock_ctx,
            self.reports_artifact,
            self.outputs_artifact,
            coverage_threshold=0.1,
        )

        # Verify get_action was called with correct arguments
        mock_ctx.get_action.assert_called_once_with("annotate", "kraken2_to_features")

        # Verify kraken2_to_features action was called with correct arguments
        mock_to_features_action.assert_called_once()
        to_features_call_args = mock_to_features_action.call_args[0]
        self.assertEqual(to_features_call_args[0], self.reports_artifact)
        self.assertEqual(to_features_call_args[1], 0.1)

        # Verify _build_contig_map was called with correct arguments
        mock_build_map.assert_called_once()
        build_map_call_args = mock_build_map.call_args[0]
        self.assertEqual(
            build_map_call_args[0],
            ANY,
        )

        # Verify missing taxa were moved to "0"
        # The function modifies mock_contig_map in place
        # After the function, "999" should be deleted and its contigs moved to "0"
        self.assertNotIn("999", mock_contig_map)
        self.assertIn("0", mock_contig_map)
        self.assertIn("contig2", mock_contig_map["0"])

        # Verify artifacts were created with correct arguments
        self.assertEqual(mock_ctx.make_artifact.call_count, 2)

        # Verify taxonomy artifact creation
        taxonomy_call_args = mock_ctx.make_artifact.call_args_list[0]
        self.assertEqual(taxonomy_call_args[0][0], "FeatureData[Taxonomy]")

        # Verify map artifact creation
        map_call_args = mock_ctx.make_artifact.call_args_list[1]
        self.assertEqual(map_call_args[0][0], "FeatureMap[TaxonomyToContigs]")
        self.assertIsInstance(map_call_args[0][1], dict)
