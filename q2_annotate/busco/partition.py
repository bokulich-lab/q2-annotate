# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os

import pandas as pd

from q2_annotate.busco.types import BUSCOResultsDirectoryFormat


def collate_busco_results(
    results: BUSCOResultsDirectoryFormat,
) -> BUSCOResultsDirectoryFormat:
    collated_results = BUSCOResultsDirectoryFormat()

    result_dfs = []
    for result in results:
        df = pd.read_csv(result.path / "busco_results.tsv", sep="\t", index_col=0)
        result_dfs.append(df)

    pd.concat(result_dfs).to_csv(
        os.path.join(collated_results.path, "busco_results.tsv"),
        sep="\t",
        index=True,
        header=True,
    )

    return collated_results
