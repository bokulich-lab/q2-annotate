# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os

from q2_types.kraken2 import Kraken2ReportDirectoryFormat, Kraken2OutputDirectoryFormat
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


def _create_filtered_results(suffix, file_dict, ids_to_keep):
    if suffix == "report":
        fmt = Kraken2ReportDirectoryFormat()
    else:
        fmt = Kraken2OutputDirectoryFormat()

    # Recreate the directory structure with only the specified ids
    for outer_id, inner_dict in file_dict.items():
        for inner_id, file_fp in inner_dict.items():
            if inner_id in ids_to_keep:
                if outer_id:
                    os.makedirs(os.path.join(str(fmt), outer_id),
                                exist_ok=True)
                duplicate(
                    src=file_dict[outer_id][inner_id],
                    dst=os.path.join(str(fmt), outer_id, f"{inner_id}.{suffix}.txt")
                )
    return fmt


def _validate_ids(file_dict_reports, file_dict_outputs):
    # Extract all inner IDs of file dicts
    inner_ids_reports = {key for inner in file_dict_reports.values() for key in inner}
    inner_ids_outputs = {key for inner in file_dict_outputs.values() for key in inner}

    # Check for ID mismatches between reports and outputs
    missing_in_reports = inner_ids_outputs - inner_ids_reports
    missing_in_outputs = inner_ids_reports - inner_ids_outputs

    if missing_in_reports or missing_in_outputs:
        error_message = (
            "There is a mismatch of IDs in the provided Kraken2 outputs and reports:\n"
        )
        if missing_in_reports:
            error_message += (
                f"IDs in outputs but missing in reports: {missing_in_reports}\n"
            )
        if missing_in_outputs:
            error_message += (
                f"IDs in reports but missing in outputs: {missing_in_outputs}\n"
            )

        raise ValueError(error_message.strip())

    return inner_ids_reports


def filter_kraken2_results(
        reports: Kraken2ReportDirectoryFormat,
        outputs: Kraken2OutputDirectoryFormat,
        metadata: Metadata = None,
        where: str = None,
        exclude_ids: bool = False,
        remove_empty: bool = False,
) -> (Kraken2ReportDirectoryFormat, Kraken2OutputDirectoryFormat):
    # Validate parameters
    _validate_parameters(metadata, remove_empty, where, exclude_ids)

    # Create file_dict for reports and outputs
    file_dict_reports = reports.file_dict()
    file_dict_outputs = outputs.file_dict()

    # Create fake outer ID if there is none
    if not any(isinstance(value, dict) for value in file_dict_reports.values()):
        file_dict_reports = {"": file_dict_reports}
        file_dict_outputs = {"": file_dict_outputs}

    # Get and validate IDs
    ids_to_keep = _validate_ids(file_dict_reports, file_dict_outputs)

    # Remove IDs that are linked to an empty report
    if remove_empty:
        ids_to_remove = _find_empty_reports(file_dict_reports)
        ids_to_keep -= ids_to_remove
        if ids_to_remove:
            print(f"Removing empty IDs: {', '.join(sorted(ids_to_remove))}")

    # Filter by metadata
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

    # Error if no IDs remain after filtering
    if len(ids_to_keep) == 0:
        raise ValueError("No IDs remain after filtering.")

    # Create filtered reports and outputs
    filtered_reports = _create_filtered_results(
        "report", file_dict_reports, ids_to_keep
    )
    filtered_outputs = _create_filtered_results(
        "output", file_dict_outputs, ids_to_keep
    )

    return filtered_reports, filtered_outputs
