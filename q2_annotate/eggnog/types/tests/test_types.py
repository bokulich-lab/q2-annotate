# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from qiime2.plugin.testing import TestPluginBase
from q2_annotate.eggnog.types import EggnogHmmerIdmap, EggnogHmmerIdmapDirectoryFmt


class TestEggnogHmmerIdmap(TestPluginBase):
    package = "q2_annotate.eggnog.tests"

    def test_hmmer_registration(self):
        self.assertRegisteredSemanticType(EggnogHmmerIdmap)

    def test_SingleAmino_semantic_type_registered_to_DirFmt(self):
        self.assertSemanticTypeRegisteredToFormat(
            EggnogHmmerIdmap, EggnogHmmerIdmapDirectoryFmt
        )
