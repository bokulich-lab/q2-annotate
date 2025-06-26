# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import tempfile
from functools import partial
from typing import Union

import pandas as pd

from q2_annotate.eggnog.orthologs.common import (
    _run_eggnog_search_pipeline,
    _eggnog_search,
    _search_runner,
)
from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_types.genome_data import SeedOrthologDirFmt, LociDirectoryFormat
from q2_types.per_sample_sequences import ContigSequencesDirFmt, MultiMAGSequencesDirFmt
from q2_types.reference_db import DiamondDatabaseDirFmt


def _eggnog_diamond_search(
    seqs: Union[ContigSequencesDirFmt, MultiMAGSequencesDirFmt, MAGSequencesDirFmt],
    db: DiamondDatabaseDirFmt,
    num_cpus: int = 1,
    db_in_memory: bool = False,
) -> (SeedOrthologDirFmt, pd.DataFrame, LociDirectoryFormat):
    with tempfile.TemporaryDirectory() as output_loc:
        db_fp = os.path.join(str(db), "ref_db.dmnd")
        search_runner = partial(
            _search_runner,
            output_loc=str(output_loc),
            num_cpus=num_cpus,
            db_in_memory=db_in_memory,
            runner_args=["diamond", "--dmnd_db", str(db_fp)],
        )
        result, ft, loci = _eggnog_search(seqs, search_runner, str(output_loc))
    return result, ft, loci


def search_orthologs_diamond(
    ctx, seqs, db, num_cpus=1, db_in_memory=False, num_partitions=None
):
    collated_hits, collated_tables, loci = _run_eggnog_search_pipeline(
        ctx,
        seqs,
        [db],
        num_cpus,
        db_in_memory,
        num_partitions,
        "_eggnog_diamond_search",
    )
    return collated_hits, collated_tables, loci
