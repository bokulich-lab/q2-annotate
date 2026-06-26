# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from . import eggnog
from . import prodigal
from .kaiju import classification as kaiju_class, database as kaiju_db
from .kraken2 import (
    classification as kraken_class,
    database as kraken_db,
    bracken,
    helpers as kraken_helpers,
)

try:
    from ._version import __version__
except ModuleNotFoundError:
    __version__ = "0.0.0+notfound"

__all__ = [
    "bracken",
    "kraken_class",
    "kraken_db",
    "kaiju_class",
    "kaiju_db",
    "eggnog",
    "prodigal",
    "kraken_helpers",
]
