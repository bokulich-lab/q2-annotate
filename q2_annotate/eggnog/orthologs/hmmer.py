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
    _create_symlinks,
    _eggnog_search,
    _search_runner,
)
from q2_annotate.eggnog.types import EggnogHmmerIdmapDirectoryFmt
from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_types.genome_data import (
    ProteinsDirectoryFormat,
    SeedOrthologDirFmt,
    LociDirectoryFormat,
)
from q2_types.per_sample_sequences import ContigSequencesDirFmt, MultiMAGSequencesDirFmt
from q2_types.profile_hmms import PressedProfileHmmsDirectoryFmt


def _eggnog_hmmer_search(
    seqs: Union[ContigSequencesDirFmt, MultiMAGSequencesDirFmt, MAGSequencesDirFmt],
    idmap: EggnogHmmerIdmapDirectoryFmt,
    pressed_hmm_db: PressedProfileHmmsDirectoryFmt,
    seed_alignments: ProteinsDirectoryFormat,
    num_cpus: int = 1,
    db_in_memory: bool = False,
) -> (SeedOrthologDirFmt, pd.DataFrame, LociDirectoryFormat):
    with tempfile.TemporaryDirectory() as output_loc:
        taxon_id = os.listdir(idmap.path)[0].split(".")[0]
        tmp_subdir = f"{output_loc}/hmmer/{taxon_id}"
        os.makedirs(tmp_subdir)
        _create_symlinks([pressed_hmm_db, idmap, seed_alignments], tmp_subdir)
        search_runner = partial(
            _search_runner,
            output_loc=output_loc,
            num_cpus=num_cpus,
            db_in_memory=db_in_memory,
            runner_args=[
                "hmmer",
                "--data_dir",
                output_loc,
                "-d",
                taxon_id,
                "--genepred",
                "prodigal",  # default incompatible with HMMER
            ],
        )
        result, ft, loci = _eggnog_search(seqs, search_runner, output_loc)
    return result, ft, loci


def search_orthologs_hmmer(
    ctx,
    seqs,
    pressed_hmm_db,
    idmap,
    seed_alignments,
    num_cpus=1,
    db_in_memory=False,
    num_partitions=None,
):
    collated_hits, collated_tables, collated_loci = _run_eggnog_search_pipeline(
        ctx,
        seqs,
        [idmap, pressed_hmm_db, seed_alignments],
        num_cpus,
        db_in_memory,
        num_partitions,
        "_eggnog_hmmer_search",
    )
    return collated_hits, collated_tables, collated_loci
