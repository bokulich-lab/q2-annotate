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
    Kraken2ReportDirectoryFormat, Kraken2OutputDirectoryFormat,
    Kraken2ReportFormat
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
        data['name'] = data['name'].strip()
        node.kraken_data = data

        # don't incorporate unclassified node into tree
        if data['name'] == 'unclassified':
            unclassified_node = node
            continue

        # find taxonomy depth
        indented_name = row['name']
        indent = 2
        indentation = _get_indentation(indented_name, indent)
        taxonomy_depth = indentation / indent

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
        [node.kraken_data['n_frags_assigned'] for node in root.traverse()]
    )

    for node in root.traverse():
        relative_abundance = node.kraken_data['n_frags_covered'] / total_reads

        if node.is_root():
            continue

        if relative_abundance < abundance_threshold:
            ancestors = node.ancestors()

            # unlink node and its descendants from tree
            node.parent.remove(node)

            # find the number of reads that were lost in the removed subtree
            reads_removed = node.kraken_data['n_frags_covered']

            for ancestor in ancestors:
                # adjust the ancestor's subtree read count
                ancestor.kraken_data['n_frags_covered'] -= reads_removed

                # check if ancestor has fallen beneath abundance threshold
                n_frags_covered = ancestor.kraken_data['n_frags_covered']
                ancestor_relative_abundance = n_frags_covered / total_reads
                if ancestor_relative_abundance < abundance_threshold:
                    if ancestor.is_root():
                        raise ValueError('Root taxon had all reads removed.')

                    # unlink the ancestor and increase the removed subtree
                    # read count
                    ancestor.parent.remove(ancestor)
                    reads_removed += n_frags_covered

    return root

def _dump_tree_to_dataframe(
    root: TreeNode, unclassified_node: TreeNode | None
) ->  pd.DataFrame:
    '''

    Parameters
    ----------

    Returns
    -------
    '''
    # create dataframe

    # if unclassified exists write to dataframe


    # perform DFS of tree
        # use n_frags_covered to determine which node to recurse to
        # (use highest first)

        # write dataframe row

def _taxonomy_dfs(node: TreeNode):
    # at given level find siblings and sort in order of decreasing
    # n_frags_covered

    # iterate over siblings and call dfs recursively

    return node
