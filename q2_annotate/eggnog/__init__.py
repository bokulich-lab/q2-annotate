# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from .orthologs import (
    search_orthologs_diamond,
    _eggnog_diamond_search,
    search_orthologs_hmmer,
    _eggnog_hmmer_search,
    _eggnog_feature_table,
)
from .annotation import map_eggnog, _eggnog_annotate, extract_annotations
from .dbs import (
    fetch_eggnog_db,
    fetch_diamond_db,
    build_custom_diamond_db,
    fetch_eggnog_proteins,
    build_eggnog_diamond_db,
    fetch_ncbi_taxonomy,
    fetch_eggnog_hmmer_db,
)

__all__ = [
    "search_orthologs_diamond",
    "_eggnog_diamond_search",
    "map_eggnog",
    "_eggnog_feature_table",
    "fetch_eggnog_db",
    "fetch_diamond_db",
    "build_custom_diamond_db",
    "fetch_eggnog_proteins",
    "build_eggnog_diamond_db",
    "fetch_ncbi_taxonomy",
    "_eggnog_annotate",
    "fetch_eggnog_hmmer_db",
    "search_orthologs_hmmer",
    "_eggnog_hmmer_search",
    "extract_annotations",
]
