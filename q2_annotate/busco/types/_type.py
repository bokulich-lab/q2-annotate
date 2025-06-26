# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from qiime2.core.type import SemanticType
from q2_types.reference_db import ReferenceDB

BUSCOResults = SemanticType("BUSCOResults")
BUSCO = SemanticType("BUSCO", variant_of=ReferenceDB.field["type"])
