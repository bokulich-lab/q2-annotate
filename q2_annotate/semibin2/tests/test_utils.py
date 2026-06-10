# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
from qiime2.plugin.testing import TestPluginBase

from q2_annotate.semibin2.utils import _process_semibin2_arg


class TestMetabat2Utils(TestPluginBase):
    package = "q2_annotate.semibin2.tests"

    def test_process_semibin2_arg_training_type_semi(self):
        obs = _process_semibin2_arg("training_type", "semi")
        exp = ["--semi-supervised"]
        self.assertListEqual(obs, exp)

    def test_process_semibin2_arg_training_type_self(self):
        obs = _process_semibin2_arg("training_type", "self")
        exp = ["--self-supervised"]
        self.assertListEqual(obs, exp)

    def test_process_semibin2_arg_option(self):
        obs = _process_semibin2_arg("epochs", 16)
        exp = ["--epochs", "16"]
        self.assertListEqual(obs, exp)

    def test_process_semibin2_arg_bool(self):
        obs = _process_semibin2_arg("no_recluster", True)
        exp = ["--no-recluster"]
        self.assertListEqual(obs, exp)


if __name__ == "__main__":
    unittest.main()
