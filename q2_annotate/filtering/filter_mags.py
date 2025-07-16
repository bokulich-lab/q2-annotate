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

    if exclude_ids:
        ids -= set(selected_ids)
    else:
        ids &= set(selected_ids)
    print(f"Found {len(ids)} IDs to keep.")

    return ids


def _filter_manifest(manifest: pd.DataFrame, ids_to_keep: set) -> pd.DataFrame:
    """
    Filters a manifest DataFrame based on a set of IDs.

    Parameters:
        manifest (pd.DataFrame): The manifest DataFrame to filter.
        ids_to_keep (set): The set of IDs to keep.

    Returns:
        pd.DataFrame: The filtered manifest DataFrame.
    """

    manifest["filename"] = (
        manifest.index.get_level_values("sample-id")
        + "/"
        + manifest.index.get_level_values("mag-id")
        + ".fasta"
    )

    return manifest[manifest.index.get_level_values("mag-id").isin(ids_to_keep)]


def _mags_to_df(sample_dict: dict):
    """
    Converts a sample dict to a DataFrame. The sample dict can be created with the
    function sample_dict() of the class MultiMAGSequencesDirFmt.

    Parameters:
        sample_dict (dict): The sample dict.

    Returns:
        pd.DataFrame: The converted DataFrame.
    """
    mags_df = pd.DataFrame.from_dict(sample_dict, orient="index")
    mags_df = mags_df.stack().reset_index()
    mags_df.columns = ["sample_id", "mag_id", "mag_fp"]
    mags_df.set_index("mag_id", inplace=True)
    return mags_df


def _find_empty_mags(mag_df) -> set:
    """
    Identifies empty FASTA files (0-byte) from a list of paths provided in a DataFrame.

    Parameters:
        mag_df (DataFrame): A DataFrame with MAG IDs as the index and a column with
        full paths to the FASTA files.

    Returns:
        set: A set of MAG IDs corresponding to empty files.
    """
    empty_mags = set()
    for mag_id, row in mag_df.iterrows():
        if os.path.getsize(row["mag_fp"]) == 0:
            empty_mags.add(mag_id)
    return empty_mags


def _validate_parameters(metadata, remove_empty):
    if not any([metadata, remove_empty]):
        raise ValueError(
            "At least one of the following parameters must be provided: "
            "metadata, remove_empty."
        )


def filter_derep_mags(
    mags: MAGSequencesDirFmt,
    metadata: Metadata = None,
    where: str = None,
    exclude_ids: bool = False,
    remove_empty: bool = False,
) -> MAGSequencesDirFmt:
    _validate_parameters(metadata, remove_empty)

    results = MAGSequencesDirFmt()
    mags_df = _mags_to_df({"": mags.feature_dict()})
    ids_to_keep = set(mags_df.index)

    if metadata is not None:
        ids_to_keep = _filter_ids(ids_to_keep, metadata, where, exclude_ids)

    if remove_empty:
        empty_mags = _find_empty_mags(mags_df)
        ids_to_keep -= empty_mags

    if len(ids_to_keep) == 0:
        raise ValueError("No MAGs remain after filtering.")

    try:
        for _id in ids_to_keep:
            duplicate(
                mags_df.loc[_id, "mag_fp"], os.path.join(str(results), f"{_id}.fasta")
            )
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
    _validate_parameters(metadata, remove_empty)

    results = MultiMAGSequencesDirFmt()
    mags_df = _mags_to_df(mags.sample_dict())
    ids_to_keep = set(mags_df.index)

    if metadata is not None:
        if on == "mag":
            ids_to_keep = _filter_ids(set(mags_df.index), metadata, where, exclude_ids)
        else:
            samples_to_keep = _filter_ids(
                set(mags_df["sample_id"]), metadata, where, exclude_ids
            )
            ids_to_keep = set(mags_df[mags_df["sample_id"].isin(samples_to_keep)].index)

    if remove_empty:
        empty_mags = _find_empty_mags(mags_df)
        ids_to_keep -= empty_mags

    if len(ids_to_keep) == 0:
        raise ValueError("No MAGs remain after filtering.")

    filtered_mags = mags_df[mags_df.index.isin(ids_to_keep)]
    filtered_manifest = _filter_manifest(mags.manifest.view(pd.DataFrame), ids_to_keep)

    filtered_manifest.to_csv(os.path.join(str(results), "MANIFEST"), sep=",")
    try:
        for _id, row in filtered_mags.iterrows():
            sample_dir = os.path.join(str(results), row["sample_id"])
            os.makedirs(sample_dir, exist_ok=True)
            duplicate(row["mag_fp"], os.path.join(sample_dir, f"{_id}.fasta"))
    except KeyError:
        raise ValueError(f"{_id!r} is not a MAG present in the input data.")

    return results
