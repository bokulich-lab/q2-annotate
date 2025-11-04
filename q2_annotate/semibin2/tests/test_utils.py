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


class TestSemiBin2Utils(TestPluginBase):
    package = "q2_annotate.semibin2.tests"

    def test_process_semibin2_arg_standard_params(self):
        """Test processing of standard SemiBin2 parameters."""
        # Test threads parameter
        obs = _process_semibin2_arg("threads", 4)
        exp = ["--threads", "4"]
        self.assertListEqual(exp, obs)

        # Test min_len parameter
        obs = _process_semibin2_arg("min_len", 1000)
        exp = ["--min-len", "1000"]
        self.assertListEqual(exp, obs)

        # Test batch_size parameter
        obs = _process_semibin2_arg("batch_size", 32)
        exp = ["--batch-size", "32"]
        self.assertListEqual(exp, obs)

        # Test epochs parameter
        obs = _process_semibin2_arg("epochs", 20)
        exp = ["--epochs", "20"]
        self.assertListEqual(exp, obs)

        # Test random_seed parameter
        obs = _process_semibin2_arg("random_seed", 42)
        exp = ["--random-seed", "42"]
        self.assertListEqual(exp, obs)

        # Test sequencing_type parameter
        obs = _process_semibin2_arg("sequencing_type", "human_gut")
        exp = ["--data-type", "human_gut"]
        self.assertListEqual(exp, obs)

    def test_process_semibin2_arg_boolean_params(self):
        """Test processing of boolean parameters."""
        # Test boolean True - should return only the flag
        obs = _process_semibin2_arg("some_flag", True)
        exp = ["--some-flag"]
        self.assertListEqual(exp, obs)

        # Test boolean False - should return flag with value
        obs = _process_semibin2_arg("some_flag", False)
        exp = ["--some-flag", "False"]
        self.assertListEqual(exp, obs)

    def test_process_semibin2_arg_unmapped_params(self):
        """Test processing of parameters not in the mapping."""
        # Test parameter not in mapping - should convert underscores to hyphens
        obs = _process_semibin2_arg("custom_param", "value")
        exp = ["--custom-param", "value"]
        self.assertListEqual(exp, obs)

        # Test parameter with multiple underscores
        obs = _process_semibin2_arg("some_complex_param_name", 123)
        exp = ["--some-complex-param-name", "123"]
        self.assertListEqual(exp, obs)

    def test_process_semibin2_arg_numeric_params(self):
        """Test processing of numeric parameters."""
        # Test integer
        obs = _process_semibin2_arg("threads", 8)
        exp = ["--threads", "8"]
        self.assertListEqual(exp, obs)

        # Test float
        obs = _process_semibin2_arg("some_float", 0.5)
        exp = ["--some-float", "0.5"]
        self.assertListEqual(exp, obs)


if __name__ == "__main__":
    unittest.main()