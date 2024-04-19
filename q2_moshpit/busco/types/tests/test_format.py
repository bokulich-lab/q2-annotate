# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from qiime2.core.exceptions import ValidationError
from qiime2.plugin.testing import TestPluginBase

from q2_moshpit.busco.types import BUSCOResultsFormat


class TestBUSCOFormats(TestPluginBase):
    package = "q2_moshpit.busco.types.tests"

    def test_busco_results_format_ok(self):
        results = BUSCOResultsFormat(
            self.get_data_path('busco_results.tsv'),
            mode='r'
        )
        results.validate(level='min')
        results.validate(level='max')

    def test_busco_results_format_error_header(self):
        results = BUSCOResultsFormat(
            self.get_data_path('busco_results_broken_header.tsv'),
            mode='r'
        )
        with self.assertRaisesRegex(ValidationError, 'Invalid header'):
            results.validate()

    def test_busco_results_format_error_values(self):
        results = BUSCOResultsFormat(
            self.get_data_path('busco_results_broken_values.tsv'),
            mode='r'
        )
        with self.assertRaisesRegex(
                ValidationError,
                'Line 4 has 14 columns'
        ):
            results.validate()
