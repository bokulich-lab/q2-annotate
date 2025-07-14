# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from qiime2.plugin.testing import TestPluginBase
from q2_types.reference_db import ReferenceDB
from q2_annotate.busco.types._format import (
    BuscoDatabaseDirFmt,
    BUSCOResultsDirectoryFormat,
)
from q2_annotate.busco.types._type import BUSCOResults, BUSCO


class TestBuscoTypes(TestPluginBase):
    package = "q2_annotate.busco.types.tests"

    def test_feature_data_semantic_type_registration(self):
        self.assertRegisteredSemanticType(BUSCOResults)

    def test_sequence_semantic_type_to_format_registration(self):
        self.assertSemanticTypeRegisteredToFormat(
            BUSCOResults, BUSCOResultsDirectoryFormat
        )

    def test_BuscoDatabaseDirFmt_registration(self):
        self.assertRegisteredSemanticType(BUSCO)

    def test_BuscoDatabaseDirFmt_semantic_type_registered_to_DirFmt(self):
        self.assertSemanticTypeRegisteredToFormat(
            ReferenceDB[BUSCO], BuscoDatabaseDirFmt
        )
