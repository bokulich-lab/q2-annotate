# ----------------------------------------------------------------------------
# Copyright (c) 2022-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from copy import copy
from pathlib import Path
import os

from skbio import TreeNode
import pandas as pd

from qiime2 import Metadata
from qiime2.util import duplicate
from q2_types.kraken2 import (
    Kraken2ReportDirectoryFormat,
    Kraken2OutputDirectoryFormat,
    Kraken2ReportFormat,
    Kraken2OutputFormat,
)

from q2_annotate.kraken2.select import _get_indentation


class ReportFilterError(ValueError):
    pass


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
                third_line = file.readline()

                # Only process if the file has exactly one line
                if not second_line or not third_line:
                    if not second_line:
                        next_line = first_line
                    else:
                        next_line = second_line
                    columns = next_line.split('\t')

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


def _filter_kraken2_results_by_metadata(
    reports: Kraken2ReportDirectoryFormat,
    outputs: Kraken2OutputDirectoryFormat,
    metadata: Metadata | None = None,
    where: str | None = None,
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


def filter_kraken2_results(
    ctx,
    reports,
    outputs,
    metadata=None,
    where=None,
    exclude_ids=False,
    remove_empty=False,
    abundance_threshold=None,
):
    '''
    Filters kraken2 reports and outputs by metadata and filters taxa within
    reports by relative abundance.

    Parameters
    ----------
    ctx : qiime2.sdk.Context
        The pipeline context object.
    reports : Kraken2ReportDirectoryFormat
        The kraken2 reports to be filtered.
    outputs : Kraken2OutputDirectoryFormat
        The kraken2 outputs to be filtered.
    metadata: qiime2.Metadata | None
        The per-sample metadata.
    where : str | None
        A SQLite where clause specifying which samples to retain.
    exlude_ids : bool
        Whether the sample IDs selected by the where clause should be exluded,
        instead of retained.
    remove_empty : bool
        Whether to remove reports that have 100% unclassified reads.
    abundance_threshold : float | None
        The relative abundance threshold beneath which taxa in each kraken2
        report will be filtered.

    Returns
    -------
    tuple[Kraken2ReportsDirectoryFormat, Kraken2OutputsDirectoryFormat]
        The filtered sets of kraken2 reports and outputs.

    Raises
    ------
    ValueError
        If all filtering parameters were left as default and no filtering can
        be perfomed.

    ValueError
        If the set of sample IDs in the kraken2 reports is not equal to the set
        of sample IDs in the kraken2 outputs.
    '''
    if (
        metadata is None and
        remove_empty is False and
        abundance_threshold is None
    ):
        msg = (
            'None of metadata-based filtering, empty report filtering, or '
            'abundance filtering were specified. Please specify at least one '
            'of these.'
        )
        raise ValueError(msg)

    if (
        set(reports.view(Kraken2ReportDirectoryFormat).file_dict().keys()) !=
        set(outputs.view(Kraken2OutputDirectoryFormat).file_dict().keys())
    ):
        msg = (
            'The sample IDs present in the input kraken2 reports did not '
            'match the sample IDs present in the input kraken2 outputs.'
        )
        raise ValueError(msg)

    _filter_kraken2_results_by_metadata = ctx.get_action(
        'annotate', '_filter_kraken2_results_by_metadata'
    )
    _filter_kraken2_reports_by_abundance = ctx.get_action(
        'annotate', '_filter_kraken2_reports_by_abundance'
    )
    _align_outputs_with_reports = ctx.get_action(
        'annotate', '_align_outputs_with_reports'
    )

    if metadata is not None or remove_empty:
        # not parallelizable because errors on all reports being filtered
        post_md_reports, post_md_outputs = _filter_kraken2_results_by_metadata(
            reports,
            outputs,
            metadata,
            where,
            exclude_ids,
            remove_empty
        )
    else:
        post_md_reports = reports
        post_md_outputs = outputs

    if abundance_threshold is None:
        return post_md_reports, post_md_outputs

    abun_f_reports, = _filter_kraken2_reports_by_abundance(
        post_md_reports, abundance_threshold, remove_empty
    )
    abun_f_outputs, = _align_outputs_with_reports(
        post_md_outputs, abun_f_reports
    )

    return abun_f_reports, abun_f_outputs


def _filter_kraken2_reports_by_abundance(
    reports: Kraken2ReportDirectoryFormat, abundance_threshold: float,
    remove_empty: bool = False
) -> Kraken2ReportDirectoryFormat:
    '''
    Filters all nodes in a kraken2 report with a relative abundance that is
    less than `abundance_threshold`. Relative abundance is calculated as the
    number of reads mapped to a tnode divided by the total number of classified
    reads in the report.

    Parameters
    ----------
    reports : Kraken2ReportDirectoryFormat
        The reports to be filtered.
    abundance_threshold : float
        The relative abundance threshold beneath which a node will be filtered.

    Returns
    -------
    Kraken2ReportDirectoryFormat
        The filtered reports.
    '''
    filtered_reports = Kraken2ReportDirectoryFormat()

    has_samples = False

    for report_filename, report in reports.reports.iter_views(
        Kraken2ReportFormat
    ):
        report = report.view(pd.DataFrame)

        # parse report into tree
        root, unclassified_node = _report_df_to_tree(report)

        if root is not None:
            total_reads = root._kraken_data['n_frags_covered']

            # trim tree
            try:
                root = _trim_tree_dfs(
                    root,
                    abundance_threshold=abundance_threshold,
                    total_reads=total_reads
                )
            except ReportFilterError:
                if remove_empty:
                    continue
                else:
                    raise

        # dump to report
        report = _dump_tree_to_report(root, unclassified_node)

        # add report to output format
        output_fp = Path(filtered_reports.path) / report_filename
        report.to_csv(output_fp, sep='\t', header=None, index=None)
        has_samples = True

    if not has_samples:
        raise ReportFilterError('All samples were filtered out.')

    # return directory format
    return filtered_reports


def _report_df_to_tree(
    report: pd.DataFrame
) -> tuple[TreeNode | None, TreeNode | None]:
    '''
    Parses a kraken2 report dataframe into a tree. Returns the root node of the
    tree and an unclassified node if unclassified reads are present.

    Parameters
    ----------
    report : pd.DataFrame
        The kraken2 report.

    Returns
    -------
    tuple[TreeNode | None, TreeNode | None]
        The root of the parsed tree, or None if the report contains only
        unclassified reads, and a node representing the unclassified
        rows, or None if no unclassified reads exist.
    '''
    most_recent_parents = []
    unclassified_node = None
    for _, row in report.iterrows():
        node = TreeNode()

        # add row data to node object
        data = {col: row[col] for col in report.columns}
        node._kraken_data = data

        # don't incorporate unclassified node into tree
        if data['name'] == 'unclassified':
            unclassified_node = node
            continue

        # find taxonomy depth
        taxonomy_depth = _get_indentation(node._kraken_data['name'], indent=2)

        # update most_recent_parents
        if len(most_recent_parents) <= taxonomy_depth:
            most_recent_parents.append(node)
        else:
            most_recent_parents[taxonomy_depth] = node

        if taxonomy_depth == 0:
            continue

        parent = most_recent_parents[taxonomy_depth - 1]
        parent.append(node)

    if len(most_recent_parents) == 0:
        return None, unclassified_node

    return most_recent_parents[0], unclassified_node


def _trim_tree_dfs(
    node: TreeNode, abundance_threshold: float, total_reads: int
) -> TreeNode | int:
    '''
    Filters nodes beneath an abundance threshold from a tree. Filtered nodes
    are removed from the tree.

    Parameters
    ----------
    node : TreeNode
        The root of the report tree.
    abundance_threshold : float
        The minimum abundance threshold required for a node to be retained.
    total_reads : int
        The total number of reads classified to the tree that `node` is in.

    Returns
    -------
    TreeNode
        The root of the filtered tree.
    '''
    relative_abundance = node._kraken_data['n_frags_covered'] / total_reads

    if relative_abundance < abundance_threshold:
        ancestors = node.ancestors()

        # unlink node and its descendants from tree
        node.parent.remove(node)

        # find the number of reads that were lost in the removed subtree
        reads_removed = node._kraken_data['n_frags_covered']

        ancestors_removed = 0
        for ancestor in ancestors:
            # adjust the ancestor's subtree read count
            ancestor._kraken_data['n_frags_covered'] -= reads_removed

            # check if ancestor has fallen beneath abundance threshold
            n_frags_covered = ancestor._kraken_data['n_frags_covered']
            ancestor_relative_abundance = n_frags_covered / total_reads

            if ancestor_relative_abundance < abundance_threshold:
                if ancestor.is_root():
                    raise ReportFilterError(
                        'Root taxon had all reads removed in at least one '
                        'of the samples. To exclude these sample(s) '
                        'from the output, re-run with the '
                        '`remove_empty` parameter.')

                # unlink the ancestor and increase the removed subtree
                # read count
                ancestor.parent.remove(ancestor)
                ancestors_removed += 1
                reads_removed += n_frags_covered

        # short circuit recursion because subtree has already been trimmed
        return ancestors_removed

    for child in copy(node.children):
        status = _trim_tree_dfs(child, abundance_threshold, total_reads)

        # short circuit to parent level because we are in a trimmed subtree
        if isinstance(status, int) and status > 0:
            return status - 1

    return node


def _dump_tree_to_report(
    root: TreeNode | None, unclassified_node: TreeNode | None
) -> pd.DataFrame:
    '''
    Recreates the kraken2 report from the filtered tree and optional
    unlcassifed node.

    Parameters
    ----------
    root : TreeNode | None
        The root node of the report tree, or None if the report is empty or
        there is only an unclassified node.
    unclassified_node : TreeNode | None
        The node represeting the unclassified row of the report, if one was
        present.

    Returns
    -------
    pd.DataFrame
        The recreated kraken2 report.
    '''
    if root is None and unclassified_node is None:
        return pd.DataFrame()

    total_reads = 0
    report = None
    if root is not None:
        report = pd.DataFrame(columns=root._kraken_data.keys())
        total_reads += root._kraken_data['n_frags_covered']
    if unclassified_node is not None:
        if report is None:
            report = pd.DataFrame(columns=unclassified_node._kraken_data.keys())

        total_reads += unclassified_node._kraken_data['n_frags_covered']

    report._kraken_total_reads = total_reads

    if unclassified_node is not None:
        _write_node_to_report(unclassified_node, report)

    if root is not None:
        _write_report_dfs(root, report)

    return report


def _write_report_dfs(node: TreeNode, report: pd.DataFrame) -> None:
    '''
    Writes nodes to the report using a depth-first search of the report tree.

    NOTE: Children are sorted in decreasing order of n_frags_covered to be
    consistent with the manner in which kraken2 structures its report files.

    Parameters
    ----------
    node : TreeNode
        The current node being explored.
    report : pd.DataFrame
        The kraken2 report dataframe.

    Returns
    -------
    None
        Writes to the `report` dataframe.
    '''
    _write_node_to_report(node, report)

    children = node.children
    children.sort(
        key=lambda n: n._kraken_data['n_frags_covered'], reverse=True
    )

    for child in children:
        _write_report_dfs(child, report)


def _write_node_to_report(node: TreeNode, report: pd.DataFrame):
    '''
    Writes a single node to a kraken2 report represented as a dataframe.

    Parameters
    ----------
    node : TreeNode
        The node for which a record will be added to the report dataframe.
    report : pd.DataFrame
        The kraken2 report represented as a dataframe.
    '''
    row = pd.Series(node._kraken_data)

    # calculate perc_frags_covered for node
    row['perc_frags_covered'] = round(
        (row['n_frags_covered'] / report._kraken_total_reads) * 100, 2
    )

    report.loc[len(report)] = row


def _align_outputs_with_reports(
    outputs: Kraken2OutputDirectoryFormat,
    reports: Kraken2ReportDirectoryFormat,
) -> Kraken2OutputDirectoryFormat:
    '''
    Filters kraken2 outputs to align them with filtered kraken2 reports.
    Assumes that the set of report sample IDs and output sample IDs are
    exactly the same.

    Parameters
    ----------
    outputs : Kraken2OutputDirectoryFormat
        The kraken2 outputs to filter.
    reports : Kraken2ReportDirectoryFormat
        The filtered kraken2 reports.

    Returns
    -------
    Kraken2OutputDirectoryFormat
        The aligned kraken2 outputs.
    '''
    aligned_outputs = Kraken2OutputDirectoryFormat()

    output_sample_map = outputs.file_dict()
    report_sample_map = reports.file_dict()
    for sample_id, report_path in report_sample_map.items():
        output_path = output_sample_map[sample_id]
        output = Kraken2OutputFormat(output_path, mode='r')
        report = Kraken2ReportFormat(report_path, mode='r')

        _align_single_output_with_report(
            output, report, aligned_outputs, sample_id
        )

    return aligned_outputs


def _align_single_output_with_report(
    output: Kraken2OutputFormat,
    report: Kraken2ReportFormat,
    output_dir_fmt: Kraken2OutputDirectoryFormat,
    sample_id: str,
) -> None:
    '''
    Filters a kraken2 output to align it with a filtered kraken2 report.
    Rows in the output that have a taxon id that is not present in the report
    are removed.

    Parameters
    ----------
    output : Kraken2OutputFormat
        The kraken2 output to filter.
    report : Kraken2ReportFormat
        The filtered kraken2 report.
    output_dir_fmt : Kraken2OutputDirectoryFormat
        The directory format to which to write the filtered output.
    sample_id : str
        The sample ID of the the report and output.

    Returns
    -------
    None
        Writes the filtered output to `output_dir_fmt`.
    '''
    report_df = report.view(pd.DataFrame)
    retained_ids = set(report_df['taxon_id'])
    output_df = output.view(pd.DataFrame)
    output_df = output_df[output_df['taxon_id'].isin(retained_ids)]

    output_fp = os.path.join(output_dir_fmt.path, f"{sample_id}.output.txt")
    output_df.to_csv(output_fp, sep='\t', header=None, index=None)
