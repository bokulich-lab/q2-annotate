# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import gzip
import os

from qiime2 import Metadata
from qiime2.util import duplicate

from q2_types.per_sample_sequences import CasavaOneEightSingleLanePerSampleDirFmt

from q2_annotate._utils import colorify
from q2_annotate.filtering.utils import _filter_ids


def _filter_empty(manifest):
    """
    Identify and return sample names from a manifest DataFrame where at least one FASTQ
    file is empty.

    Parameters:
        manifest (pandas.DataFrame): A DataFrame where the index contains sample names
            and the values are file paths pointing to FASTQ files.

    Returns:
        list: A list of sample names for which at least one associated file is empty
            (no content on the first line).
    """
    empty_samples = []
    for sample, row in manifest.iterrows():
        for path in row:
            with gzip.open(path, "rt") as f:
                if not f.readline():
                    empty_samples.append(sample)
                    break
    return empty_samples


def filter_reads(
    reads: CasavaOneEightSingleLanePerSampleDirFmt,
    metadata: Metadata = None,
    where: str = None,
    exclude_ids: bool = False,
    remove_empty: bool = False,
) -> CasavaOneEightSingleLanePerSampleDirFmt:

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

    results = CasavaOneEightSingleLanePerSampleDirFmt()
    manifest = reads.manifest
    samples_to_keep = set(manifest.index)

    # If metadata is provided, filter the IDs based on the metadata
    if metadata is not None:
        samples_to_keep = _filter_ids(set(manifest.index), metadata, where, exclude_ids)

    # Remove empty samples if requested
    if remove_empty:
        empty_samples = _filter_empty(manifest)
        samples_to_keep = [i for i in samples_to_keep if i not in empty_samples]
        print(
            colorify(
                f"Removed {len(empty_samples)} empty samples: "
                f"{', '.join(empty_samples)}"
            )
        )

    if not samples_to_keep:
        raise ValueError(
            "There are no samples left after filtering. Please check your metadata and "
            "filtering criteria."
        )

    # Copy the files that are kept to the results directory
    for sample, row in manifest.iterrows():
        for path in row:
            if sample in samples_to_keep:
                try:
                    duplicate(path, os.path.join(results.path, os.path.basename(path)))
                except KeyError:
                    raise ValueError(f"{sample!r} is not a sample present in the data.")

    return results
