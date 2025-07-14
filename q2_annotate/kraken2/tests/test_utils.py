# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import unittest

import pandas as pd
from qiime2.plugin.testing import TestPluginBase

from q2_annotate.kraken2.utils import _process_kraken2_arg, _find_lca, _join_ranks


class TestKraken2Utils(TestPluginBase):
    package = "q2_annotate.kraken2.tests"

    def setUp(self):
        super().setUp()
        self.m_flor = [
            [
                "Bacteria",
                "Bacteria",
                "Actinomycetota",
                "Actinomycetes",
                "Mycobacteriales",
                "Mycobacteriaceae",
                "Mycobacterium",
                "Mycobacterium florentinum",
            ]
        ] * 5
        self.m_methano = [
            [
                "Bacteria",
                "Bacteria",
                "Actinomycetota",
                "Actinomycetes",
                "Mycobacteriales",
                "Mycobacteriaceae",
                "Mycobacterium",
                "Candidatus Mycobacterium methanotrophicum",
            ]
        ] * 2
        self.c_glu_r = [
            [
                "Bacteria",
                "Bacteria",
                "Actinomycetota",
                "Actinomycetes",
                "Mycobacteriales",
                "Corynebacteriaceae",
                "Corynebacterium",
                "Corynebacterium glutamicum",
                "Corynebacterium glutamicum ATCC R",
            ]
        ] * 2

        self.taxa = pd.Series(
            [
                *self.m_flor,
                *self.m_methano,
                [
                    "Bacteria",
                    "Bacteria",
                    "Actinomycetota",
                    "Actinomycetes",
                    "Mycobacteriales",
                    "Mycobacteriaceae",
                    "Mycobacterium",
                    "Mycobacterium avium",
                    "Mycobacterium avium subsp. hominissuis",
                ],
                [
                    "Bacteria",
                    "Bacteria",
                    "Actinomycetota",
                    "Actinomycetes",
                    "Mycobacteriales",
                    "Mycobacteriaceae",
                    "Mycobacterium",
                    "Mycobacterium heckeshornense",
                ],
                [
                    "Bacteria",
                    "Bacteria",
                    "Actinomycetota",
                    "Actinomycetes",
                    "Mycobacteriales",
                    "Mycobacteriaceae",
                    "Mycobacterium",
                    "Mycobacterium heckeshornense",
                ],
            ],
            name="Taxon",
        )

        self.taxa_mixed = pd.Series(
            [
                *self.m_methano,
                *self.c_glu_r,
                [
                    "Bacteria",
                    "Bacteria",
                    "Actinomycetota",
                    "Actinomycetes",
                    "Mycobacteriales",
                    "Corynebacteriaceae",
                    "Corynebacterium",
                    "Corynebacterium glutamicum",
                    "Corynebacterium glutamicum ATCC 13032",
                ],
            ],
            name="Taxon",
        )

    def test_process_kraken2_arg_bool(self):
        obs = _process_kraken2_arg("quick", True)
        exp = ["--quick"]
        self.assertListEqual(obs, exp)

    def test_process_kraken2_arg_number(self):
        obs = _process_kraken2_arg("threads", 3)
        exp = ["--threads", "3"]
        self.assertListEqual(obs, exp)

    def test_process_kraken2_arg_string(self):
        obs = _process_kraken2_arg("db", "/some/where/test_db")
        exp = ["--db", "/some/where/test_db"]
        self.assertListEqual(obs, exp)

    def test_process_kraken2_arg_unknown(self):
        with self.assertRaisesRegex(
            NotImplementedError,
            "Parsing arguments of type \"<class 'list'>\" " "is not supported.",
        ):
            _process_kraken2_arg("fake_param", [1, 2])

    def test_find_lca_1(self):
        obs = list(_find_lca(self.taxa))
        exp = [
            "Bacteria",
            "Bacteria",
            "Actinomycetota",
            "Actinomycetes",
            "Mycobacteriales",
            "Mycobacteriaceae",
            "Mycobacterium",
        ]
        self.assertListEqual(obs, exp)

    def test_find_lca_2(self):
        obs = list(_find_lca(self.taxa_mixed))
        exp = [
            "Bacteria",
            "Bacteria",
            "Actinomycetota",
            "Actinomycetes",
            "Mycobacteriales",
        ]
        self.assertListEqual(obs, exp)

    # def test_find_lca_majority_1(self):
    #     obs = list(_find_lca_majority(self.taxa))
    #     exp = [
    #         'Bacteria', 'Bacteria', 'Actinomycetota', 'Actinomycetes',
    #         'Mycobacteriales', 'Mycobacteriaceae', 'Mycobacterium',
    #         'Mycobacterium florentinum', None
    #     ]
    #     self.assertListEqual(obs, exp)
    #
    # def test_find_lca_majority_2(self):
    #     obs = list(_find_lca_majority(self.taxa_mixed))
    #     exp = [
    #         'Bacteria', 'Bacteria', 'Actinomycetota', 'Actinomycetes',
    #         'Mycobacteriales', 'Corynebacteriaceae', 'Corynebacterium',
    #         'Corynebacterium glutamicum',
    #     ]
    #     self.assertListEqual(obs, exp)

    def test_join_ranks_full(self):
        ranks = ["d__", "k__", "p__", "c__", "o__", "f__", "g__", "s__", "ssp__"]
        obs = _join_ranks(self.c_glu_r[0], ranks)
        exp = (
            "d__Bacteria;k__Bacteria;p__Actinomycetota;c__Actinomycetes;"
            "o__Mycobacteriales;f__Corynebacteriaceae;g__Corynebacterium;"
            "s__Corynebacterium glutamicum;"
            "ssp__Corynebacterium glutamicum ATCC R"
        )
        self.assertEqual(obs, exp)

    def test_join_ranks_no_ssp(self):
        ranks = ["d__", "k__", "p__", "c__", "o__", "f__", "g__", "s__"]
        obs = _join_ranks(self.c_glu_r[0], ranks)
        exp = (
            "d__Bacteria;k__Bacteria;p__Actinomycetota;c__Actinomycetes;"
            "o__Mycobacteriales;f__Corynebacteriaceae;g__Corynebacterium;"
            "s__Corynebacterium glutamicum"
        )
        self.assertEqual(obs, exp)

    def test_join_ranks_with_terminal_nones(self):
        ranks = [
            "d__",
            "k__",
            "p__",
        ]
        obs = _join_ranks(["Bacteria", "Bacteria", "Actinomycetota", None], ranks)
        exp = "d__Bacteria;k__Bacteria;p__Actinomycetota"
        self.assertEqual(obs, exp)

    def test_join_ranks_with_other_nones(self):
        ranks = ["d__", "k__", "p__", "c__"]
        obs = _join_ranks(["Bacteria", "Bacteria", None, "Actinomycetes"], ranks)
        exp = "d__Bacteria;k__Bacteria;p__;c__Actinomycetes"
        self.assertEqual(obs, exp)

    def test_join_ranks_with_nones_all_over(self):
        ranks = ["d__", "k__", "p__", "c__"]
        obs = _join_ranks(["Bacteria", "Bacteria", None, "Actinomycetes", None], ranks)
        exp = "d__Bacteria;k__Bacteria;p__;c__Actinomycetes"
        self.assertEqual(obs, exp)


if __name__ == "__main__":
    unittest.main()
