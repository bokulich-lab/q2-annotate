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

from q2_annotate.filtering.utils import _filter_ids


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


def _find_empty_mags(file_dict: dict) -> set:
    """
    Identifies empty MAG files (0-byte) in a file dict.

    Parameters:
        file_dict (dict): A nested dictionary with keys and full paths created by
                          sample_dict or feature_dict functions.

    Returns:
        set: A set of MAG IDs corresponding to empty files.
    """
    empty_mags = set()
    for sample_id, mag_dict in file_dict.items():
        for mag_id, path in mag_dict.items():
            if os.path.getsize(path) == 0:
                empty_mags.add(mag_id)
    return empty_mags


def filter_derep_mags(
    mags: MAGSequencesDirFmt,
    metadata: Metadata = None,
    where: str = None,
    exclude_ids: bool = False,
    remove_empty: bool = False,
) -> MAGSequencesDirFmt:
    results = MAGSequencesDirFmt()
    features = mags.feature_dict()
    ids_to_keep = set(features.keys())

    if metadata is not None:
        ids_to_keep = _filter_ids(set(features.keys()), metadata, where, exclude_ids)

    if remove_empty:
        empty_mags = _find_empty_mags({"": features})
        ids_to_keep -= empty_mags

    try:
        for _id in ids_to_keep:
            duplicate(features[_id], os.path.join(str(results), f"{_id}.fasta"))
    except KeyError:
        raise ValueError(f"{_id!r} is not a MAG present in the input data.")

    return results


def filter_mags(
    mags: MultiMAGSequencesDirFmt,
    metadata: Metadata = None,
    where: str = None,
    exclude_ids: bool = False,
    on: str = "mag",
    remove_empty: bool = False,
) -> MultiMAGSequencesDirFmt:
    results = MultiMAGSequencesDirFmt()
    mags_df = _mags_to_df(mags, on)
    sample_dict = mags.sample_dict()
    ids_to_keep = set(mags_df.index)

    if metadata is not None:
        ids_to_keep = _filter_ids(set(mags_df.index), metadata, where, exclude_ids)

    if remove_empty:
        empty_mags = _find_empty_mags(sample_dict)
        ids_to_keep -= empty_mags

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
