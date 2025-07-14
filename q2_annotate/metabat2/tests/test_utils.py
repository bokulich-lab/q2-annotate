# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
from qiime2.plugin.testing import TestPluginBase

from q2_annotate.metabat2.utils import _process_metabat2_arg


class TestMetabat2Utils(TestPluginBase):
    package = "q2_annotate.metabat2.tests"

    def test_process_metabat_arg_tnf(self):
        obs = _process_metabat2_arg("p_tnf", 12)
        exp = ["--pTNF", "12"]
        self.assertListEqual(obs, exp)

    def test_process_metabat_arg_min_cv(self):
        obs = _process_metabat2_arg("min_cv", 10)
        exp = ["--minCV", "10"]
        self.assertListEqual(obs, exp)

    def test_process_metabat_arg_min_cv_sum(self):
        obs = _process_metabat2_arg("min_cv_sum", 50)
        exp = ["--minCVSum", "50"]
        self.assertListEqual(obs, exp)

    def test_process_metabat_arg_other(self):
        obs = _process_metabat2_arg("max_edges", 200)
        exp = ["--maxEdges", "200"]
        self.assertListEqual(obs, exp)

    def test_process_metabat_arg_bool(self):
        obs = _process_metabat2_arg("no_add", True)
        exp = ["--noAdd"]
        self.assertListEqual(obs, exp)


if __name__ == "__main__":
    unittest.main()
