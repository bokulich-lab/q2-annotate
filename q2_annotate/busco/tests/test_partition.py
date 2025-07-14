# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os

import pandas as pd
from qiime2.plugin.testing import TestPluginBase

from q2_annotate.busco.partition import collate_busco_results
from q2_annotate.busco.types import BUSCOResultsDirectoryFormat


class TestBUSCOPlots(TestPluginBase):
    package = "q2_annotate.busco.tests"

    def test_collate_busco_results(self):
        p1 = self.get_data_path("busco_results/partition1")
        p2 = self.get_data_path("busco_results/partition2")

        busco_results = [
            BUSCOResultsDirectoryFormat(p1, mode="r"),
            BUSCOResultsDirectoryFormat(p2, mode="r"),
        ]

        collated_busco_result = collate_busco_results(busco_results)

        obs = pd.read_csv(os.path.join(str(collated_busco_result), "busco_results.tsv"))
        exp = pd.read_csv(
            self.get_data_path("busco_results/results_all/busco_results.tsv")
        )

        pd.testing.assert_frame_equal(obs, exp)
