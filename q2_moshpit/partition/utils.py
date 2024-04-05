# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from typing import List
import warnings


def _validate_num_partitions(num_mags: int, num_partitions: int) -> int:

    if num_partitions is None:
        return num_mags
    elif num_partitions > num_mags:
        warnings.warn(
            "You have requested a number of partitions"
            f" '{num_partitions}' that is greater than your number"
            f" of MAGs '{num_mags}.' Your data will be"
            f" partitioned by MAG into '{num_mags}'"
            " partitions."
        )
        return num_mags
    else:
        return num_partitions


def _validate_mag_ids(
    num_partitions: int, num_mags: int, mags_all: List[tuple]
):
    # If num_partitions == num_mags and MAG ids are not unique
    # the output will be missing these duplicated-id MAGs.
    # While this is technically impossible since
    # MAGs should have unique IDs by construction, it could still happen that a
    # used imports MAGs with non-unique IDs. In such case this test would be
    # useful.

    if num_partitions == num_mags:
        mag_ids = [mag_id[1] for mag_id in mags_all]
        duplicates = [
            mag_id for mag_id in mag_ids if mag_ids.count(mag_id) > 1
        ]
        if len(duplicates) > 0:
            raise ValueError(
                "MAG IDs are not unique. "
                "They must be unique in order to output all partitions "
                "correctly. Printing duplicate MAG IDs: "
                f"{set(duplicates)}"
            )
