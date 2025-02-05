# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os

from q2_types.kraken2 import Kraken2ReportDirectoryFormat
from qiime2 import Metadata
from qiime2.util import duplicate


def _validate_parameters(metadata, remove_empty, where, exclude_ids):
    if not metadata and not remove_empty:
        raise ValueError('Please specify parameters "--m-metadata-file" or '
                         '"--p-remove-empty"  to filter accordingly.')

    if where and not metadata:
        raise ValueError('The parameter "--p-where" can only be specified in '
                         'combination with the parameter '
                         '"--m-metadata-file".')

    if exclude_ids and not metadata:
        raise ValueError('The parameter "--p-exclude-ids" can only be '
                         'specified in combination with the parameter '
                         '"--m-metadata-file".')


def _find_empty_reports(file_dict: dict) -> set:
    empty_ids = set()
    for inner_dict in file_dict.values():
        for inner_id, file_fp in inner_dict.items():
            with open(file_fp, 'r') as file:
                # Read the first line and check if there's a second line
                first_line = file.readline().strip()
                second_line = file.readline()

                # Only process if the file has exactly one line
                if not second_line:
                    columns = first_line.split('\t')

                    # Check if the 6th column contains "unclassified" or
                    # "root"
                    if len(columns) > 5 and columns[5] in ["unclassified",
                                                           "root"]:
                        empty_ids.add(inner_id)

    return empty_ids


def _create_filtered_results(file_dict, ids_to_keep):
    # Specify output format and file name suffix
    results = Kraken2ReportDirectoryFormat()
    suffix = ".report"

    for outer_id, inner_dict in file_dict.items():
        for inner_id, file_fp in inner_dict.items():
            if inner_id in ids_to_keep:
                if outer_id:
                    os.makedirs(os.path.join(str(results), outer_id),
                                exist_ok=True)
                duplicate(
                    src=file_dict[outer_id][inner_id],
                    dst=os.path.join(str(results), outer_id,
                                     f"{inner_id}{suffix}.txt")
                )

    return results


def filter_kraken_reports(
        reports: Kraken2ReportDirectoryFormat,
        metadata: Metadata = None,
        where: str = None,
        exclude_ids: bool = False,
        remove_empty: bool = False,
) -> Kraken2ReportDirectoryFormat:
    # Validate parameters
    _validate_parameters(metadata, remove_empty, where, exclude_ids)

    # Create file_dict
    file_dict = reports.file_dict()

    # Create fake outer ID if there is none, to make it easier to iterate
    if not any(isinstance(value, dict) for value in file_dict.values()):
        file_dict = {"": file_dict}

    # Extract all inner IDs
    ids_to_keep = {key for inner in file_dict.values() for key in inner}

    # Remove IDs that are linked to an empty report
    if remove_empty:
        ids_to_remove = _find_empty_reports(file_dict)
        ids_to_keep -= ids_to_remove
        if ids_to_remove:
            print(f"Removing empty IDs: {', '.join(sorted(ids_to_remove))}")

    if metadata:
        selected_ids = metadata.get_ids(where=where)
        if not selected_ids:
            print("The filter query returned no IDs to filter out.")

        if not (set(selected_ids) - ids_to_keep):
            print(f"IDs {', '.join(sorted(set(selected_ids) - ids_to_keep))} "
                  f"are not present in the data.")

        if exclude_ids:
            ids_to_keep -= set(selected_ids)
        else:
            ids_to_keep &= set(selected_ids)

    if len(ids_to_keep) == 0:
        raise ValueError("No IDs remain after filtering.")

    return _create_filtered_results(file_dict, ids_to_keep)
