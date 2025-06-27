# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import shutil

from q2_types.kraken2 import Kraken2ReportDirectoryFormat, Kraken2OutputDirectoryFormat


def collate_kraken2_reports(
    reports: Kraken2ReportDirectoryFormat,
) -> Kraken2ReportDirectoryFormat:
    collated_kraken2_reports = Kraken2ReportDirectoryFormat()
    _collate_kraken_tsvs(reports, collated_kraken2_reports)
    return collated_kraken2_reports


def collate_kraken2_outputs(
    outputs: Kraken2OutputDirectoryFormat,
) -> Kraken2OutputDirectoryFormat:
    collated_kraken2_outputs = Kraken2OutputDirectoryFormat()
    _collate_kraken_tsvs(outputs, collated_kraken2_outputs)
    return collated_kraken2_outputs


def _collate_kraken_tsvs(results, output):
    for result in results:
        for fp in result.path.iterdir():
            shutil.move(str(fp), output.path / os.path.basename(fp))
