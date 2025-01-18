# ----------------------------------------------------------------------------
# Copyright (c) 2022-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from skbio import TreeNode
import pandas as pd

from q2_types.kraken2 import (
    Kraken2ReportDirectoryFormat,
    Kraken2OutputDirectoryFormat,
    Kraken2ReportFormat,
)

from q2_annotate.kraken2.select import _get_indentation


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


def _filter_tree(root: TreeNode, abundance_threshold: float) -> TreeNode:
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
    # calculate denominator for relative abundance
    total_reads = sum(
        [node._kraken_data['n_frags_assigned'] for node in root.traverse()]
    )

    for node in root.traverse():
        relative_abundance = node._kraken_data['n_frags_covered'] / total_reads

        if node.is_root():
            continue

        if relative_abundance < abundance_threshold:
            ancestors = node.ancestors()

            # unlink node and its descendants from tree
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

    return root


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
    columns = root._kraken_data.keys()
    report = pd.DataFrame(data=0, columns=columns)

    # if unclassified exists write to dataframe
    if unclassified_node:
        _write_node_to_report(unclassified_node, report)

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
    # write node to report
    row = pd.Series(node._kraken_data)
    report = pd.concat(report, row, axis='index')

    children = node.children
    children.sort(
        key=lambda n: n._kraken_data['n_frags_covered'], reverse=True
    )

    for child in children:
        _write_report_dfs(child, report)
