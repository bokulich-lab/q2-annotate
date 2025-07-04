# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from qiime2 import Metadata


def _filter_ids(
    ids: set, metadata: Metadata = None, where: str = None, exclude_ids: bool = False
) -> set:
    """
    Filters IDs based on the provided metadata.

    Parameters:
        ids (set): The set of IDs to filter.
        metadata (Metadata, optional): The metadata to use for filtering.
            Defaults to None.
        where (str, optional): The condition to use for filtering.
            Defaults to None.
        exclude_ids (bool, optional): Whether to exclude the IDs that
            match the condition. Defaults to False.

    Returns:
        set: The filtered set of IDs.
    """
    selected_ids = metadata.get_ids(where=where)
    if not selected_ids:
        print("The filter query returned no IDs to filter out.")
    else:
        if exclude_ids:
            ids -= set(selected_ids)
        else:
            ids &= set(selected_ids)
    print(f"Found {len(ids)} IDs to keep.")

    return ids


def _validate_parameters(metadata, where, remove_empty):
    if not any([metadata, remove_empty]):
        raise ValueError(
            "At least one of the following parameters must be provided: "
            "metadata, remove_empty."
        )

    if metadata and not where:
        raise ValueError(
            "A filter query must be provided through the 'where' parameter "
            "when filtering by metadata."
        )
