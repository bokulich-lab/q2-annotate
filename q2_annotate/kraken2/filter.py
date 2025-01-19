# ----------------------------------------------------------------------------
# Copyright (c) 2022-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from copy import copy
from pathlib import Path

from skbio import TreeNode
import pandas as pd

from q2_types.kraken2 import (
    Kraken2ReportDirectoryFormat,
    Kraken2OutputDirectoryFormat,
    Kraken2ReportFormat,
)

from q2_annotate.kraken2.select import _get_indentation


def filter_kraken2_reports(
    reports: Kraken2ReportDirectoryFormat,
    abundance_threshold: float,
) -> Kraken2ReportDirectoryFormat:
    '''

    Parameters
    ----------

    Returns
    -------
    '''
    # make output directory format
    output_dir_format = Kraken2ReportDirectoryFormat()

    # iterate over input formats
    for report_filename, report in reports.reports.iter_views(
        Kraken2ReportFormat
    ):
        report = report.view(pd.DataFrame)

        # parse report into tree
        root, unclassified_node = _report_df_to_tree(report)

        # calculate total reads
        total_reads = root._kraken_data['n_frags_covered']

        # trim tree
        root_trimmed = _trim_tree_dfs(
            root,
            abundance_threshold=abundance_threshold,
            total_reads=total_reads
        )

        # dump to report
        trimmed_report = _dump_tree_to_report(
            root_trimmed, unclassified_node
        )

        # add report to output format
        output_fp = Path(output_dir_format.path) / report_filename
        trimmed_report.to_csv(
            output_fp, sep='\t', header=None, index=None
        )

    # return directory format
    return output_dir_format


def _report_df_to_tree(
    report: pd.DataFrame
) -> tuple[TreeNode, TreeNode | None]:
    '''
    Parses a kraken2 report dataframe into a tree. Returns the root node of the
    tree and an unclassified node if unclassified reads are present.

    Parameters
    ----------
    report : pd.DataFrame
        The kraken2 report.

    Returns
    -------
    tuple[TreeNode, TreeNode | None]
        The root of the parsed tree and a node representing the unclassified
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

    return most_recent_parents[0], unclassified_node


def _trim_tree_dfs(
    node: TreeNode, abundance_threshold: float, total_reads: int
) -> TreeNode | None:
    '''
    Filters nodes beneath an abundance threshold from a tree. Filtered nodes
    are removed from the tree.

    Parameters
    ----------
    root : TreeNode
        The root of the report tree.
    abundance_threshold : float
        The minimum abundance threshold required for a node to be retained.

    Returns
    -------
    TreeNode
        The root of the filtered tree.
    '''
    relative_abundance = node._kraken_data['n_frags_covered'] / total_reads

    if relative_abundance < abundance_threshold:
        ancestors = node.ancestors()

        # unlink node and its descendants from tree
        root = list(node.ancestors())[-1]

        node.parent.remove(node)

        # find the number of reads that were lost in the removed subtree
        reads_removed = node._kraken_data['n_frags_covered']

        for ancestor in ancestors:
            # adjust the ancestor's subtree read count
            ancestor._kraken_data['n_frags_covered'] -= reads_removed

            # check if ancestor has fallen beneath abundance threshold
            n_frags_covered = ancestor._kraken_data['n_frags_covered']
            ancestor_relative_abundance = n_frags_covered / total_reads

            if ancestor_relative_abundance < abundance_threshold:
                if ancestor.is_root():
                    raise ValueError('Root taxon had all reads removed.')

                # unlink the ancestor and increase the removed subtree
                # read count
                ancestor.parent.remove(ancestor)
                reads_removed += n_frags_covered

        # short circuit recursion because subtree has already been trimmed
        return


    for child in copy(node.children):
        _trim_tree_dfs(child, abundance_threshold, total_reads)

    return node

def _dump_tree_to_report(
    root: TreeNode, unclassified_node: TreeNode | None
) ->  pd.DataFrame:
    '''
    Recreates the kraken2 report from the filtered tree and optional
    unlcassifed node.

    Parameters
    ----------
    root : TreeNode
        The root node of the report tree.
    unclassified_node : TreeNode | None
        The node represeting the unclassified row of the report, if one was
        present.

    Returns
    -------
    pd.DataFrame
        The recreated kraken2 report.
    '''
    # create dataframe
    report = pd.DataFrame(columns=root._kraken_data.keys())

    # if unclassified exists write to dataframe
    if unclassified_node:
        _write_node_to_report(unclassified_node, report)

    # calculate denominator for perc_frags_covered_column
    total_reads = root._kraken_data['n_frags_covered']
    if unclassified_node:
        total_reads += unclassified_node._kraken_data['n_frags_covered']

    report._kraken_total_reads = total_reads

    # write tree to report
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

    Parameters
    ----------

    Returns
    -------
    '''
    row = pd.Series(node._kraken_data)

    # calculate perc_frags_covered for node
    row['perc_frags_covered'] = round(
        (row['n_frags_covered'] / report._kraken_total_reads) * 100, 2
    )

    report.loc[len(report)] = row
