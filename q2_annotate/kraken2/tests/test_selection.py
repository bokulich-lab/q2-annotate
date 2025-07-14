# ----------------------------------------------------------------------------
# Copyright (c) 2023-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import shutil
import tempfile
import unittest

import pandas as pd
import pandas.testing
from pandas._testing import assert_frame_equal
import skbio
from q2_annotate.kraken2 import kraken2_to_features
from q2_annotate.kraken2.select import (
    _kraken_to_ncbi_tree,
    _find_lcas,
    kraken2_to_mag_features,
)
from qiime2.plugin.testing import TestPluginBase

from q2_types.kraken2 import (
    Kraken2ReportDirectoryFormat,
    Kraken2OutputDirectoryFormat,
)


class MockTempDir(tempfile.TemporaryDirectory):
    pass


class TestKrakenSelect(TestPluginBase):
    package = "q2_annotate.kraken2.tests"

    def setUp(self):
        super().setUp()
        self.temp_dir = tempfile.mkdtemp()

        fp = self.get_data_path("kraken2-reports-select/kraken2_reads_table.json")
        self.kraken2_reads_table = pd.read_json(fp)
        self.kraken2_reads_table.columns = [
            str(x) for x in self.kraken2_reads_table.columns
        ]

        fp = self.get_data_path(
            "kraken2-reports-select/kraken2_reads_table_filtered.json"
        )
        self.kraken2_reads_table_filtered = pd.read_json(fp)
        self.kraken2_reads_table_filtered.columns = [
            str(x) for x in self.kraken2_reads_table_filtered.columns
        ]

        fp = self.get_data_path("kraken2-reports-select/kraken2_taxonomy.tsv")
        self.kraken_taxonomy = pd.read_csv(
            fp, sep="\t", header=0, dtype={"Feature ID": "object", "Taxon": str}
        )
        self.kraken_taxonomy.set_index("Feature ID", inplace=True)

        fp = self.get_data_path("kraken2-reports-select/kraken2_taxonomy_filtered.tsv")
        self.kraken_taxonomy_filtered = pd.read_csv(
            fp, sep="\t", header=0, dtype={"Feature ID": "object", "Taxon": str}
        )
        self.kraken_taxonomy_filtered.set_index("Feature ID", inplace=True)

        self.taxa_mag1 = pd.read_csv(self.get_data_path("mag-taxa-1.csv"), index_col=0)
        self.taxa_mag2 = pd.read_csv(self.get_data_path("mag-taxa-2.csv"), index_col=0)
        self.taxa_mag3 = pd.read_csv(self.get_data_path("mag-taxa-3.csv"), index_col=0)
        self.taxa_mag4 = pd.read_csv(self.get_data_path("mag-taxa-4.csv"), index_col=0)
        self.taxa_mag5 = pd.read_csv(self.get_data_path("mag-taxa-5.csv"), index_col=0)

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_kraken2_to_features_infra_clade(self):
        reports = Kraken2ReportDirectoryFormat(self.get_data_path("infra-clade/"), "r")
        obs_table, obs_taxonomy = kraken2_to_features(reports, coverage_threshold=0.0)
        # Check that expected taxons(tip of the tree) are in the taxonomy/table
        assert "237084" in obs_taxonomy.index
        assert "237084" in obs_table.columns
        assert "5199" in obs_taxonomy.index
        assert "5199" in obs_table.columns

    def test_kraken2_to_features_duplicated_genus(self):
        reports = Kraken2ReportDirectoryFormat(
            self.get_data_path("duplicated-genus/"), "r"
        )
        obs_table, obs_taxonomy = kraken2_to_features(reports, coverage_threshold=0.0)
        # 51977 is the taxon code for pilophorus which shares a genus name with
        # taxon code 237084
        assert "51977" in obs_taxonomy.index
        assert "51977" in obs_table.columns
        assert "237084" in obs_taxonomy.index
        assert "237084" in obs_table.columns

        # Taxon codes 51977 and 5199 are sister genera and 5199 needs to be in
        # the first evaluated tree and 51977 must be absent
        # to induce this duplicated genus name issue in the second tree
        assert "5199" in obs_taxonomy.index
        assert "5199" in obs_table.columns

    def test_kraken2_to_features_coverage_threshold(self):
        reports = Kraken2ReportDirectoryFormat(
            self.get_data_path("kraken2-reports-select/samples"), "r"
        )
        obs_table, obs_taxonomy = kraken2_to_features(reports, coverage_threshold=0.1)

        assert_frame_equal(obs_table, self.kraken2_reads_table_filtered)
        assert_frame_equal(obs_taxonomy, self.kraken_taxonomy_filtered)

    def test_kraken2_to_features_no_coverage_threshold(self):
        reports = Kraken2ReportDirectoryFormat(
            self.get_data_path("kraken2-reports-select/samples"), "r"
        )
        obs_table, obs_taxonomy = kraken2_to_features(reports, 0.0)

        assert_frame_equal(obs_table, self.kraken2_reads_table)
        assert_frame_equal(obs_taxonomy, self.kraken_taxonomy)

    def test_kraken_to_ncbi_tree(self):
        reports = Kraken2ReportDirectoryFormat(
            self.get_data_path("kraken2-to-ncbi-tree/example1"), "r"
        )
        reports = list(reports.reports.iter_views(pd.DataFrame))
        report_df = reports[0][1]

        exp_tree_fp = self.get_data_path("kraken2-to-ncbi-tree/exp-tree1.txt")
        exp_tree = skbio.TreeNode.read(exp_tree_fp)

        obs_tree = _kraken_to_ncbi_tree(report_df)

        # skbio.TreeNode doesn't define an equality operator, so performing
        # this test on newick strings (there is only one ordering, so unless
        # skbio changes its iteration order, this should be stable enough)
        self.assertEqual(str(obs_tree), str(exp_tree))

    def test_kraken2_to_mag_features(self):
        reports = Kraken2ReportDirectoryFormat(self.get_data_path("reports-mags"), "r")
        hits = Kraken2OutputDirectoryFormat(self.get_data_path("outputs-mags"), "r")
        obs = kraken2_to_mag_features(reports, hits, 0.0)
        exp = pd.DataFrame.from_dict(
            {
                "3b72d1a7-ddb0-4dc7-ac36-080ceda04aaa": [
                    "d__Bacteria;k__Bacteria;p__Firmicutes;c__Bacilli"
                ],
                "8894435a-c836-4c18-b475-8b38a9ab6c6b": ["d__Unclassified"],
            },
            orient="index",
        )
        exp.columns = ["Taxon"]
        exp.index.name = "Feature ID"
        pandas.testing.assert_frame_equal(obs, exp)

    def test_kraken2_to_mag_features_incorrect_fraction_unclassified(self):
        reports = Kraken2ReportDirectoryFormat(
            self.get_data_path("reports-mags-unclassified-wrong-frac"), "r"
        )
        hits = Kraken2OutputDirectoryFormat(self.get_data_path("outputs-mags"), "r")
        with self.assertRaisesRegex(
            ValueError,
            "fraction for MAG '8894435a-c836-4c18-b475-8b38a9ab6c6b' "
            "is not .* 99.01%",
        ):
            kraken2_to_mag_features(reports, hits, 0.0)

    def test_kraken2_to_mag_features_unclassified_root(self):
        reports = Kraken2ReportDirectoryFormat(
            self.get_data_path("reports-mags-unclassified-root"), "r"
        )
        hits = Kraken2OutputDirectoryFormat(self.get_data_path("outputs-mags"), "r")
        obs_taxonomy = kraken2_to_mag_features(reports, hits, 0.0)
        exp_taxonomy = pd.DataFrame(
            {
                "Taxon": [
                    "d__Bacteria;k__Bacteria;p__Firmicutes;c__Bacilli",
                    "d__Unclassified",
                ]
            },
            index=pd.Index(
                [
                    "3b72d1a7-ddb0-4dc7-ac36-080ceda04aaa",
                    "8894435a-c836-4c18-b475-8b38a9ab6c6b",
                ],
                name="Feature ID",
            ),
        )

        pandas.testing.assert_frame_equal(obs_taxonomy, exp_taxonomy)

    def test_kraken2_to_mag_features_unclassified_root_low_coverage(self):
        reports = Kraken2ReportDirectoryFormat(
            self.get_data_path("reports-mags-unclassified-root"), "r"
        )
        hits = Kraken2OutputDirectoryFormat(self.get_data_path("outputs-mags"), "r")
        obs_taxonomy = kraken2_to_mag_features(reports, hits, 1.0)
        exp_taxonomy = pd.DataFrame(
            {
                "Taxon": [
                    "d__Bacteria;k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;"
                    "f__Staphylococcaceae;g__Staphylococcus;s__Staphylococcus aureus",
                    "d__Unclassified",
                ]
            },
            index=pd.Index(
                [
                    "3b72d1a7-ddb0-4dc7-ac36-080ceda04aaa",
                    "8894435a-c836-4c18-b475-8b38a9ab6c6b",
                ],
                name="Feature ID",
            ),
        )

        pandas.testing.assert_frame_equal(obs_taxonomy, exp_taxonomy)

    def test_kraken2_to_mag_features_unclassified_root_only(self):
        reports = Kraken2ReportDirectoryFormat(
            self.get_data_path("reports-mags-unclassified-missing-frac"), "r"
        )
        hits = Kraken2OutputDirectoryFormat(self.get_data_path("outputs-mags"), "r")
        obs_taxonomy = kraken2_to_mag_features(reports, hits, 0.0)
        exp_taxonomy = pd.DataFrame(
            {
                "Taxon": [
                    "d__Bacteria;k__Bacteria;p__Firmicutes;c__Bacilli",
                    "d__Unclassified",
                ]
            },
            index=pd.Index(
                [
                    "3b72d1a7-ddb0-4dc7-ac36-080ceda04aaa",
                    "8894435a-c836-4c18-b475-8b38a9ab6c6b",
                ],
                name="Feature ID",
            ),
        )

        pandas.testing.assert_frame_equal(obs_taxonomy, exp_taxonomy)

    def test_kraken2_to_mag_features_unclassified_no_add_up(self):
        reports = Kraken2ReportDirectoryFormat(
            self.get_data_path("reports-mags-unclassified-no-add-up"), "r"
        )
        hits = Kraken2OutputDirectoryFormat(self.get_data_path("outputs-mags"), "r")

        with self.assertRaisesRegex(
            ValueError,
            "fraction for MAG '8894435a-c836-4c18-b475-8b38a9ab6c6b' " "is not 100.0",
        ):
            kraken2_to_mag_features(reports, hits, 0.1)

    def test_kraken2_to_mag_features_unclassified_all(self):
        reports = Kraken2ReportDirectoryFormat(
            self.get_data_path("reports-mags-unclassified-all"), "r"
        )
        hits = Kraken2OutputDirectoryFormat(self.get_data_path("outputs-mags"), "r")

        with self.assertRaisesRegex(
            ValueError, "The resulting feature table was empty"
        ):
            kraken2_to_mag_features(reports, hits, 0.1)

    def test_find_lcas_mode_lca(self):
        taxa = [self.taxa_mag1, self.taxa_mag2, self.taxa_mag3, self.taxa_mag4]
        obs = _find_lcas(taxa, mode="lca")
        exp = pd.DataFrame.from_dict(
            {
                "0e514d88-16c4-4273-a1df-1a360eb2c823": [
                    "d__Bacteria;k__Bacteria;p__Actinomycetota;c__Actinomycetes;"
                    "o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium"
                ],
                "3acec411-b0d0-4441-b936-5b8b571fa328": [
                    "d__Bacteria;k__Bacteria;p__Actinomycetota;c__Actinomycetes;"
                    "o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium"
                ],
                "3af39f47-e90a-46f2-9b5f-b236ae6551f0": [
                    "d__Bacteria;k__Bacteria;p__Actinomycetota;c__Actinomycetes;"
                    "o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium"
                ],
                "fed92059-3222-4573-b0ec-726c49fbfabb": [
                    "d__Bacteria;k__Bacteria;p__Actinomycetota;c__Actinomycetes;"
                    "o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium;"
                    "s__Mycobacterium florentinum"
                ],
            },
            orient="index",
        )
        exp.columns = ["Taxon"]
        exp.index.name = "Feature ID"
        pandas.testing.assert_frame_equal(obs, exp)

    def test_find_lcas_mode_lca_unclassified(self):
        taxa = [self.taxa_mag1, self.taxa_mag2, self.taxa_mag5]
        obs = _find_lcas(taxa, mode="lca")
        exp = pd.DataFrame.from_dict(
            {
                "0e514d88-16c4-4273-a1df-1a360eb2c823": [
                    "d__Bacteria;k__Bacteria;p__Actinomycetota;c__Actinomycetes;"
                    "o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium"
                ],
                "3acec411-b0d0-4441-b936-5b8b571fa328": [
                    "d__Bacteria;k__Bacteria;p__Actinomycetota;c__Actinomycetes;"
                    "o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium"
                ],
                "fed92059-3222-4573-b0ec-726c49fbfabb": ["d__Unclassified"],
            },
            orient="index",
        )
        exp.columns = ["Taxon"]
        exp.index.name = "Feature ID"
        pandas.testing.assert_frame_equal(obs, exp)

    # def test_find_lcas_mode_majority(self):
    #     taxa = [
    #       self.taxa_mag1, self.taxa_mag2, self.taxa_mag3, self.taxa_mag4
    #     ]
    #     obs = _find_lcas(taxa, mode='majority')
    #     exp = pd.DataFrame.from_dict({
    #         '0e514d88-16c4-4273-a1df-1a360eb2c823': [
    #             'd__Bacteria;k__Bacteria;p__Actinomycetota;c__Actinomycetes;'
    #             'o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium,'
    #             's__Mycobacterium avium;'
    #             'ssp__Mycobacterium avium subsp. hominissuis'
    #         ],
    #         '3acec411-b0d0-4441-b936-5b8b571fa328': [
    #             'd__Bacteria;k__Bacteria;p__Actinomycetota;c__Actinomycetes;'
    #             'o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium;'
    #             's__Mycobacterium florentinum'
    #         ],
    #         '3af39f47-e90a-46f2-9b5f-b236ae6551f0': [
    #             'd__Bacteria;k__Bacteria;p__Actinomycetota;c__Actinomycetes;'
    #             'o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium;'
    #             's__Mycobacterium avium;'
    #             'ssp__Mycobacterium avium subsp. hominissuis'
    #         ],
    #         'fed92059-3222-4573-b0ec-726c49fbfabb': [
    #             'd__Bacteria;k__Bacteria;p__Actinomycetota;c__Actinomycetes;'
    #             'o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium;'
    #             's__Mycobacterium florentinum'
    #         ]
    #     }, orient='index')
    #     exp.columns = ['Taxon']
    #     exp.index.name = 'Feature ID'
    #     pandas.testing.assert_frame_equal(obs, exp)
    #
    # def test_find_lcas_mode_majority_terminal_none(self):
    #     taxa = [self.taxa_mag2, self.taxa_mag4]
    #     obs = _find_lcas(taxa, mode='majority')
    #     exp = pd.DataFrame.from_dict({
    #         '3acec411-b0d0-4441-b936-5b8b571fa328': [
    #             'd__Bacteria;k__Bacteria;p__Actinomycetota;c__Actinomycetes;'
    #             'o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium;'
    #             's__Mycobacterium florentinum'
    #         ],
    #         'fed92059-3222-4573-b0ec-726c49fbfabb': [
    #             'd__Bacteria;k__Bacteria;p__Actinomycetota;c__Actinomycetes;'
    #             'o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium;'
    #             's__Mycobacterium florentinum'
    #         ]
    #     }, orient='index')
    #     exp.columns = ['Taxon']
    #     exp.index.name = 'Feature ID'
    #     pandas.testing.assert_frame_equal(obs, exp)


class TestKrakenSelectEdgeCases(unittest.TestCase):

    def make_dirfmt(self, string, coverage=False):
        """
        This is a trivial format so that the indentation and rank can be
        easily inspected by the reader.
        """
        rows = []
        for idx, line in enumerate(string.split("\n")):
            line = line.strip()
            if not line:
                continue
            tax_id = idx  # we can pretend
            if coverage:
                cover, rank, taxonomy_fragment = line.split(";")
            else:
                rank, taxonomy_fragment = line.split(";")
                cover = 0.5
            rows.append(
                dict(
                    perc_frags_covered=cover,
                    pad1=0,
                    pad2=0,
                    rank=rank,
                    taxon_id=tax_id,
                    name=taxonomy_fragment,
                )
            )
        df = pd.DataFrame(rows)

        dirfmt = Kraken2ReportDirectoryFormat()
        with open(dirfmt.path / "sample1.report.txt", "w") as fh:
            df.to_csv(fh, header=False, index=False, sep="\t")

        return dirfmt

    def make_exp(self, rows):
        taxonomy = pd.DataFrame(rows, columns=["Feature ID", "Taxon"])
        taxonomy["Feature ID"] = taxonomy["Feature ID"].apply(str)
        taxonomy = taxonomy.set_index("Feature ID")

        table = pd.DataFrame({str(idx): True for idx, _ in rows}, index=["sample1"])
        return taxonomy, table

    def test_kraken_to_ncbi_tree_simple(self):
        dirfmt = self.make_dirfmt(
            """
        R;root
        D;  a
        D;  b
        D;  c
        """
        )
        exp_tax, exp_table = self.make_exp(
            [
                (2, "d__a"),
                (3, "d__b"),
                (4, "d__c"),
            ]
        )

        table, taxonomy = kraken2_to_features(dirfmt)

        pandas.testing.assert_frame_equal(exp_table, table)
        pandas.testing.assert_frame_equal(exp_tax, taxonomy)

    def test_kraken_to_ncbi_tree_no_tricks(self):
        dirfmt = self.make_dirfmt(
            """
        R;root
        D;  a
        K;    a.a
        K;    a.b
        D;  b
        K;    b.a
        D;  c
        """
        )
        exp_tax, exp_table = self.make_exp(
            [
                (3, "d__a;k__a.a"),
                (4, "d__a;k__a.b"),
                (6, "d__b;k__b.a"),
                (7, "d__c"),
            ]
        )

        table, taxonomy = kraken2_to_features(dirfmt)

        pandas.testing.assert_frame_equal(exp_table, table)
        pandas.testing.assert_frame_equal(exp_tax, taxonomy)

    def test_kraken_to_ncbi_tree_no_nested_end(self):
        dirfmt = self.make_dirfmt(
            """
        R;root
        D;  a
        K;    a.a
        K;    a.b
        D;  b
        K;    b.a
        P;      b.a.a
        """
        )
        exp_tax, exp_table = self.make_exp(
            [
                (3, "d__a;k__a.a"),
                (4, "d__a;k__a.b"),
                (7, "d__b;k__b.a;p__b.a.a"),
            ]
        )

        table, taxonomy = kraken2_to_features(dirfmt)

        pandas.testing.assert_frame_equal(exp_table, table)
        pandas.testing.assert_frame_equal(exp_tax, taxonomy)

    def test_kraken_to_ncbi_tree_infraclade_ranks(self):
        dirfmt = self.make_dirfmt(
            """
         R;root
         D;  a
        D1;    a.s
         K;    a.b
        K1;      a.b.s
         P;        a.b.s.a
         D;  b
         K;    b.a
         P;      b.a.a
        P1;        b.a.a.a
        """
        )
        exp_tax, exp_table = self.make_exp(
            [
                (6, "d__a;k__a.b;p__a.b.s.a"),
                (9, "d__b;k__b.a;p__b.a.a"),
            ]
        )

        table, taxonomy = kraken2_to_features(dirfmt)

        pandas.testing.assert_frame_equal(exp_table, table)
        pandas.testing.assert_frame_equal(exp_tax, taxonomy)

    def test_kraken_to_ncbi_tree_subspecies_and_skip(self):
        dirfmt = self.make_dirfmt(
            """
         R;root
         K;  a
        K1;    a.s
         P;    a.b
        P1;      a.b.s
         C;        a.b.s.a
         S;          a.b.s.a.a
        S1;            a.b.s.a.a.a
        """
        )

        exp_tax, exp_table = self.make_exp(
            [
                (
                    8,
                    "d__containing k__a;k__a;p__a.b;c__a.b.s.a;"
                    "o__containing s__a.b.s.a.a;f__containing s__a.b.s.a.a;"
                    "g__containing s__a.b.s.a.a;s__a.b.s.a.a;s1__a.b.s.a.a.a",
                ),
            ]
        )

        table, taxonomy = kraken2_to_features(dirfmt)

        pandas.testing.assert_frame_equal(exp_table, table)
        pandas.testing.assert_frame_equal(exp_tax, taxonomy)

    def test_kraken_to_ncbi_tree_with_coverage_filter(self):
        dirfmt = self.make_dirfmt(
            """
        .5;R;root
        .4;K;  a
        .4;P;    a.a
        .0;P;    a.b
        .2;K;  b
        .0;P;    b.a
        .0;K;  c
        """,
            coverage=True,
        )
        exp_tax, exp_table = self.make_exp(
            [
                (3, "d__containing k__a;k__a;p__a.a"),
                (5, "d__containing k__b;k__b"),
            ]
        )

        table, taxonomy = kraken2_to_features(dirfmt)

        pandas.testing.assert_frame_equal(exp_table, table)
        pandas.testing.assert_frame_equal(exp_tax, taxonomy)

    def test_kraken_to_ncbi_tree_kingdom_promotion(self):
        dirfmt = self.make_dirfmt(
            """
        R;root
        D;  Bacteria
        P;    A
        D;  Archaea
        P;    B
        """
        )
        exp_tax, exp_table = self.make_exp(
            [
                (3, "d__Bacteria;k__Bacteria;p__A"),
                (5, "d__Archaea;k__Archaea;p__B"),
            ]
        )

        table, taxonomy = kraken2_to_features(dirfmt)

        pandas.testing.assert_frame_equal(exp_table, table)
        pandas.testing.assert_frame_equal(exp_tax, taxonomy)
