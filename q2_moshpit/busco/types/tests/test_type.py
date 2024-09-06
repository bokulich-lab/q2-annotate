# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from qiime2.plugin.testing import TestPluginBase
from q2_types.reference_db import ReferenceDB
from q2_moshpit.busco.types._format import (
    BuscoDatabaseDirFmt, BUSCOResultsDirectoryFormat
)
from q2_moshpit.busco.types._type import (
    BUSCOResults, BuscoDB
)


class TestBuscoTypes(TestPluginBase):
    package = "q2_moshpit.busco.types.tests"

    def test_feature_data_semantic_type_registration(self):
        self.assertRegisteredSemanticType(BUSCOResults)

    def test_sequence_semantic_type_to_format_registration(self):
        self.assertSemanticTypeRegisteredToFormat(
            BUSCOResults, BUSCOResultsDirectoryFormat
        )

    def test_BuscoDatabaseDirFmt_registration(self):
        self.assertRegisteredSemanticType(BuscoDB)

    def test_BuscoDatabaseDirFmt_semantic_type_registered_to_DirFmt(self):
        self.assertSemanticTypeRegisteredToFormat(
                ReferenceDB[BuscoDB],
                BuscoDatabaseDirFmt)
