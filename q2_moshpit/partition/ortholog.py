# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import os
import shutil
import warnings

import numpy as np
from q2_types.genome_data import SeedOrthologDirFmt
from qiime2.util import duplicate


def collate_orthologs(orthologs: SeedOrthologDirFmt) -> SeedOrthologDirFmt:
    result = SeedOrthologDirFmt()

    for ortholog in orthologs:
        for fp in ortholog.path.iterdir():
            shutil.move(fp, result.path / os.path.basename(fp))

    return result


def partition_orthologs(
        orthologs: SeedOrthologDirFmt, num_partitions: int = None
) -> SeedOrthologDirFmt:
    """
    Returns a dictionary where each key is either the sample_id and
    values are the new objects with the orthologs.
    """
    partitioned_orthologs = {}

    # TODO: this logic should move to the format itself
    orthologs = glob.glob(os.path.join(str(orthologs), "*.seed_orthologs"))
    names = [
        os.path.basename(x).replace(".emapper.seed_orthologs", "")
        for x in orthologs
    ]
    orthologs = list(zip(names, orthologs))

    num_samples = len(orthologs)
    if num_partitions is None:
        num_partitions = num_samples
    elif num_partitions > num_samples:
        warnings.warn(
            "You have requested a number of partitions"
            f" '{num_partitions}' that is greater than your number"
            f" of samples '{num_samples}.' Your data will be"
            f" partitioned by sample into '{num_samples}'"
            " partitions."
        )
        num_partitions = num_samples

    orthologs = np.array_split(orthologs, num_partitions)
    for i, samples in enumerate(orthologs, 1):
        result = SeedOrthologDirFmt()

        for sample_id, sample_fp in samples:
            duplicate(sample_fp, result.path / os.path.basename(sample_fp))

        # If num_partitions == num_samples we will only have gone through one
        # sample in the above loop and will use its id as a key. Otherwise we
        # may have gone through multiple samples in the above loop and will be
        # using indices for keys
        if num_partitions == num_samples:
            partitioned_orthologs[sample_id] = result
        else:
            partitioned_orthologs[i] = result

    return partitioned_orthologs
