# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os

import pandas as pd
from qiime2 import Metadata
from qiime2.util import duplicate

from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt


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


def _filter_manifest(
    manifest: pd.DataFrame, ids_to_keep: set, on: str = "mag"
) -> pd.DataFrame:
    """
    Filters a manifest DataFrame based on a set of IDs.

    Parameters:
        manifest (pd.DataFrame): The manifest DataFrame to filter.
        ids_to_keep (set): The set of IDs to keep.
        on (str): The level on which to filter ('mag' or 'sample').
            Defaults to 'mag'.

    Returns:
        pd.DataFrame: The filtered manifest DataFrame.
    """
    if on == "mag":
        lvl = "mag-id"
    elif on == "sample":
        lvl = "sample-id"
    else:
        raise ValueError(f"Invalid value for 'on' parameter: {on}")

    manifest["filename"] = (
        manifest.index.get_level_values("sample-id")
        + "/"
        + manifest.index.get_level_values("mag-id")
        + ".fasta"
    )

    return manifest[manifest.index.get_level_values(lvl).isin(ids_to_keep)]


def _mags_to_df(mags: MultiMAGSequencesDirFmt, on: str):
    """
    Converts a MultiMAGSequencesDirFmt object to a DataFrame.

    Parameters:
        mags (MultiMAGSequencesDirFmt): The MultiMAGSequencesDirFmt
            object to convert.
        on (str): The level on which to index the DataFrame
            ('sample' or 'mag').

    Returns:
        pd.DataFrame: The converted DataFrame.
    """
    mags_df = pd.DataFrame.from_dict(mags.sample_dict(), orient="index")
    mags_df = mags_df.stack().reset_index()
    mags_df.columns = ["sample_id", "mag_id", "mag_fp"]
    if on == "sample":
        mags_df.set_index("sample_id", inplace=True)
    elif on == "mag":
        mags_df.set_index("mag_id", inplace=True)
    return mags_df


def filter_derep_mags(
    mags: MAGSequencesDirFmt,
    metadata: Metadata,
    where: str = None,
    exclude_ids: bool = False,
) -> MAGSequencesDirFmt:
    results = MAGSequencesDirFmt()
    features = mags.feature_dict()
    ids_to_keep = _filter_ids(set(features.keys()), metadata, where, exclude_ids)
    try:
        for _id in ids_to_keep:
            duplicate(features[_id], os.path.join(str(results), f"{_id}.fasta"))
    except KeyError:
        raise ValueError(f"{_id!r} is not a MAG present in the input data.")

    return results


def filter_mags(
    mags: MultiMAGSequencesDirFmt,
    metadata: Metadata,
    where: str = None,
    exclude_ids: bool = False,
    on: str = "mag",
) -> MultiMAGSequencesDirFmt:
    results = MultiMAGSequencesDirFmt()
    mags_df = _mags_to_df(mags, on)

    ids_to_keep = _filter_ids(set(mags_df.index), metadata, where, exclude_ids)

    filtered_mags = mags_df[mags_df.index.isin(ids_to_keep)]
    filtered_manifest = _filter_manifest(
        mags.manifest.view(pd.DataFrame), ids_to_keep, on=on
    )
    filtered_manifest.to_csv(os.path.join(str(results), "MANIFEST"), sep=",")
    try:
        for _id, row in filtered_mags.iterrows():
            if on == "mag":
                sample_dir = os.path.join(str(results), row["sample_id"])
                mag_dest = os.path.join(sample_dir, f"{_id}.fasta")
            else:
                sample_dir = os.path.join(str(results), _id)
                mag_dest = os.path.join(sample_dir, f"{row['mag_id']}.fasta")
            os.makedirs(sample_dir, exist_ok=True)
            duplicate(row["mag_fp"], mag_dest)
    except KeyError:
        raise ValueError(f"{_id!r} is not a MAG present in the input data.")

    return results
