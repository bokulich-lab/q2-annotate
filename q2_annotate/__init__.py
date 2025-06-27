# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from . import abundance
from . import busco
from . import eggnog
from . import prodigal
from .dereplication import dereplicate_mags
from .filtering import filter_derep_mags, filter_mags, filter_reads_pangenome
from .kaiju import classification as kaiju_class, database as kaiju_db
from .kraken2 import (
    classification as kraken_class,
    database as kraken_db,
    bracken,
    helpers as kraken_helpers,
)
from .metabat2 import metabat2
from ._utils import (
    get_feature_lengths,
    multiply_tables,
    _multiply_tables,
    _multiply_tables_pa,
    _multiply_tables_relative,
)

try:
    from ._version import __version__
except ModuleNotFoundError:
    __version__ = "0.0.0+notfound"

__all__ = [
    "metabat2",
    "bracken",
    "kraken_class",
    "kraken_db",
    "kaiju_class",
    "kaiju_db",
    "dereplicate_mags",
    "eggnog",
    "busco",
    "prodigal",
    "kraken_helpers",
    "filter_derep_mags",
    "filter_mags",
    "get_feature_lengths",
    "multiply_tables",
    "_multiply_tables",
    "_multiply_tables_pa",
    "_multiply_tables_relative",
    "abundance",
    "filter_reads_pangenome",
]
