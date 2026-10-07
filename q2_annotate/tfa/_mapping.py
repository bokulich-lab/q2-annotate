# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import biom
import pandas as pd


def _align_gene_taxonomy(
    feature_load: biom.Table, gene_taxonomy: pd.DataFrame
) -> pd.DataFrame:
    """Join a mapping to table rows, allowing mappings for filtered-out rows."""
    feature_ids = list(feature_load.ids(axis="observation"))
    missing = [
        identifier
        for identifier in feature_ids
        if identifier not in gene_taxonomy.index
    ]
    if missing:
        raise ValueError(
            f"Missing gene taxonomy mappings for feature IDs: {missing!r}."
        )
    return gene_taxonomy.loc[feature_ids]
