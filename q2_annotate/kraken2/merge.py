# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import functools

import pandas as pd
from skbio import TreeNode

from q2_types.kraken2 import (
    Kraken2ReportDirectoryFormat,
    Kraken2OutputDirectoryFormat,
    Kraken2ReportFormat,
    Kraken2OutputFormat,
)

from q2_annotate.kraken2.filter import (
    _report_df_to_tree,
    _dump_tree_to_report,
)


def _merge_kraken2_results(
    reports: Kraken2ReportDirectoryFormat,
    outputs: Kraken2OutputDirectoryFormat
) -> (Kraken2ReportDirectoryFormat, Kraken2OutputDirectoryFormat):
    '''
    Merges kraken2 reports and outputs on a per-sample-id basis.

    Parameters
    ----------
    reports : list[Kraken2ReportDirectoryFormat]
        The kraken2 reports.
    outputs : list[Kraken2OutputDirectoryFromat]
        The kraken2 outputs.

    Returns
    -------
    tuple[Kraken2ReportDirectoryFormat, Kraken2OutputDirectoryFormat]
        The merged reports and formats.
    '''
    merged_reports = Kraken2ReportDirectoryFormat()
    merged_outputs = Kraken2OutputDirectoryFormat()

    mags = isinstance(next(iter(reports[0].file_dict().values())), dict)

    report_mapping, output_mapping = _condense_formats(reports, outputs, mags)

    for sample_id in report_mapping:
        if mags:
            for filename, report in report_mapping[sample_id].items():
                merged_reports.reports.write_data(
                    view=report,
                    view_type=Kraken2ReportFormat,
                    sample_id=sample_id,
                    mag_id=filename
                )

            for filename, output in output_mapping[sample_id].items():
                merged_outputs.outputs.write_data(
                    view=output,
                    view_type=Kraken2OutputFormat,
                    sample_id=sample_id,
                    mag_id=filename
                )
        else:
            reports = report_mapping[sample_id]
            merged_report = _merge_reports(reports)
            merged_reports.reports.write_data(
                view=merged_report,
                view_type=Kraken2ReportFormat,
                sample_id=sample_id
            )

            outputs = output_mapping[sample_id]
            merged_output = _merge_outputs(outputs)
            merged_outputs.outputs.write_data(
                view=merged_output,
                view_type=Kraken2OutputFormat,
                sample_id=sample_id
            )

    return merged_reports, merged_outputs


def _condense_formats(
    reports: list[Kraken2ReportDirectoryFormat],
    outputs: list[Kraken2OutputDirectoryFormat],
    mags: bool
) -> tuple[dict, dict]:
    '''
    Condenses multiple report and output directory formats into a single
    output mapping and a single report mapping. The structure is
    sample_id -> list[format] for reads/contigs and
    sample_id -> {filename -> format} for mags.

    Note that MAGs will never be merged on a per-MAG basis, only on a
    per-sample-id basis.

    Parameters
    ----------
    reports : list[Kraken2ReportDirectoryFormat]
        The kraken2 reports.
    outputs : list[Kraken2OutputDirectoryFormat]
        The kraken2 outputs.
    mags : bool
        Whether the directory formats represent MAG results.

    Returns
    -------
    tuple[dict, dict]
        A tuple of mappings as described above. The first contains the kraken2
        report mapping and the second the kraken2 output mapping.

    Raises
    ------
    ValueError
        If two MAGs with the same uuid are detected.
    ValueError
        If two reports with the same sample ID are to be merged but the
        minimizers columns are present in the reports.
    '''
    minimizers_present = _check_for_minimizers(reports)

    chained_reports = []
    for report in reports:
        for sample_id, value in report.file_dict().items():
            chained_reports.append((sample_id, value))

    chained_outputs = []
    for output in outputs:
        for sample_id, value in output.file_dict().items():
            chained_outputs.append((sample_id, value))

    def _update_mapping(sample_id, value, mapping, Format):
        if mags:
            filename_to_filepath = value
            for filename, filepath in filename_to_filepath.items():
                format = Format(filepath, mode='r')
                if sample_id not in mapping:
                    mapping[sample_id] = {filename: format}
                else:
                    if filename in mapping[sample_id]:
                        msg = (
                            'Two MAGs with the same uuid were detected. '
                            f'Duplicated uuid: {filename}.'
                        )
                        raise ValueError(msg)

                    mapping[sample_id][filename] = format
        else:
            filepath = value
            format = Format(filepath, mode='r')
            if sample_id not in mapping:
                mapping[sample_id] = [format]
            elif minimizers_present:
                msg = (
                    'Two or more reports with the same sample ID were '
                    'attempted to be merged but the option to capture minimzer '
                    'data was enabled. It is not possible to merge kraken2 '
                    'reports that contain minimizer information.'
                )
                raise ValueError(msg)
            else:
                mapping[sample_id].append(format)

    report_mapping = {}
    for sample_id, value in chained_reports:
        _update_mapping(
            sample_id, value, report_mapping, Kraken2ReportFormat
        )

    output_mapping = {}
    for sample_id, value in chained_outputs:
        _update_mapping(
            sample_id, value, output_mapping, Kraken2OutputFormat
        )

    return report_mapping, output_mapping


def _check_for_minimizers(reports: list[Kraken2ReportDirectoryFormat]) -> bool:
    '''
    Checks whether the reports being processed include the optional minimizer
    columns. This determines whether reports with a shared sample ID can be
    merged.

    Parameters
    ----------
    reports : list[Kraken2ReportDirectoryFormat]
        The kraken2 report directory formats to examine for minimzer columns.

    Returns
    -------
    bool
        Wether minimizer columns are present in the reports.
    '''
    for report_dir_format in reports:
        _, report = list(
            report_dir_format.reports.iter_views(Kraken2ReportFormat)
        )[0]

        if 'n_read_minimizers' in report.view(pd.DataFrame):
            return True

    return False


def _merge_reports(
    reports: list[Kraken2ReportFormat]
) -> Kraken2ReportFormat:
    '''
    Merges two or more kraken2 reports into a single report. See
    `_merge_trees` for the actual merging algorithm.

    Parameters
    ----------
    reports : list[Kraken2ReportFormat]
        Kraken2 reports belonging to the same sample ID.

    Returns
    -------
    Kraken2ReportFormat
        The merged kraken2 report format.
    '''
    if len(reports) == 1:
        return reports[0]

    merged = functools.reduce(
        _merge_trees,
        [_report_df_to_tree(report.view(pd.DataFrame)) for report in reports]
    )

    merged_report_df = _dump_tree_to_report(*merged)
    merged_report = Kraken2ReportFormat()
    merged_report_df.to_csv(
        str(merged_report), sep='\t', header=False, index=False
    )

    return merged_report


def _merge_outputs(
    outputs: list[Kraken2OutputFormat]
) -> Kraken2OutputFormat:
    '''
    Merges two or more kraken2 outputs into a single output. Output files
    are merged by concatenation.

    Parameters
    ----------
    outputs : list[Kraken2OutputFormat]
        Kraken2 outputs belonging to the same sample ID.

    Returns
    -------
        The merged kraken2 output format.
    '''
    if len(outputs) == 1:
        return outputs[0]

    merged_output = Kraken2OutputFormat()

    with open(str(merged_output), 'w') as merged_fh:
        while outputs:
            with open(str(outputs.pop()), 'r') as fh:
                while buffer := fh.read(4096):
                    merged_fh.write(buffer)

    return merged_output


def _merge_trees(
    first: tuple[TreeNode | None, TreeNode | None],
    second: tuple[TreeNode | None, TreeNode | None]
) -> tuple[TreeNode | None, TreeNode | None]:
    '''
    Merges two trees each representing a kraken2 report into a single tree.
    The number of reads assigned to each node are summed where nodes overlap,
    and new nodes are inserted into the tree where they don't. The proportions
    of assigned and covered reads are then updated in a final passover.

    Parameters
    ----------
    first : tuple[TreeNode | None, TreeNode | None]
        The first report tree, where the first node in the tuple represents the
        tree and the second node represents an optional unclassified node.
    second : tuple[TreeNode | None, TreeNode | None]
        The second report tree, where the first node in the tuple represents the
        tree and the second node represents an optional unclassified node.

    Returns
    -------
    tuple[TreeNode | None, TreeNode | None]
        The merged tree.
    '''
    first_tree, first_unclassified_node = first
    second_tree, second_unclassified_node = second

    # merge trees (with respect to `n_frags_assigned`, `n_frags_covered`)
    if first_tree is None and second_tree is None:
        merged_tree = None
    elif first_tree is None:
        merged_tree = second_tree
    elif second_tree is None:
        merged_tree = first_tree
    else:
        _merge_trees_recursively(from_node=first_tree, into_node=second_tree)
        merged_tree = second_tree

    unclassified_node = _merge_unclassified_nodes(
        first_unclassified_node, second_unclassified_node
    )

    # final passover to update `perc_frags_covered`
    if merged_tree is not None:
        total_reads = merged_tree._kraken_data['n_frags_covered']
        if unclassified_node is not None:
            total_reads += unclassified_node._kraken_data['n_frags_covered']

        for node in merged_tree.traverse():
            node._kraken_data['perc_frags_covered'] = round(
                (node._kraken_data['n_frags_covered'] / total_reads) * 100, 2
            )

    return merged_tree, unclassified_node


def _merge_trees_recursively(
    from_node: TreeNode, into_node: TreeNode
) -> None:
    '''
    Merges the tree rooted at `from_node` into the tree rooted at `into_node`.
    After this function completes, the `into_node` node represents the merged
    tree.

    Parameters
    ----------
    from_node : TreeNode
        The root of a tree to be merged into the tree rooted by `into_node`.
    into_node : TreeNode
        The root of a tree into which the tree represented by `from_node`
        will be merged.
    '''
    into_node._kraken_data['n_frags_covered'] += \
        from_node._kraken_data['n_frags_covered']
    into_node._kraken_data['n_frags_assigned'] += \
        from_node._kraken_data['n_frags_assigned']

    for child in from_node.children:
        match = _find_node_match(child, into_node.children)
        if match is None:
            into_node.append(child)
            child.parent = into_node
        else:
            _merge_trees_recursively(child, match)


def _find_node_match(
    search_node: TreeNode, children: list[TreeNode]
) -> TreeNode | None:
    '''
    Searches for a match to `search_node` in `children`. The nodes' taxon ids
    are used to define matching.

    Parameters
    ----------
    search_node: TreeNode
        The node for which to search in `children`.
    children: list[TreeNode]
        The nodes in which to search for `search_node`.

    Returns
    -------
    TreeNode | None
        The matching node or None if no match is found.

    Raises
    ------
    ValueError
        If more than one matching node is found.
    '''
    search_taxon_id = search_node._kraken_data['taxon_id']

    matches = []
    for child in children:
        child_taxon_id = child._kraken_data['taxon_id']
        if search_taxon_id == child_taxon_id:
            matches.append(child)

    if len(matches) > 1:
        raise ValueError('Did not expect more than one taxon id match.')
    elif len(matches) == 1:
        return matches[0]
    else:
        return None


def _merge_unclassified_nodes(
    first: TreeNode | None, second: TreeNode | None
) -> TreeNode | None:
    '''
    Merges two TreeNodes representing unclassified nodes from kraken2 reports.

    Parameters
    ----------
    first : TreeNode | None
        The first unclassified node, or None if no such node exists.
    second : TreeNode | None
        The second unclassified node, or None if no such node exists.

    Returns
    -------
    TreeNode | None
        The merged unclassified node, or None if both inputs are None.
    '''
    if first is None and second is None:
        unclassified_node = None
    elif first is None:
        unclassified_node = second
    elif second is None:
        unclassified_node = first
    else:
        unclassified_node = second
        unclassified_node._kraken_data['n_frags_assigned'] += \
            first._kraken_data['n_frags_assigned']
        unclassified_node._kraken_data['n_frags_covered'] += \
            first._kraken_data['n_frags_covered']

    return unclassified_node
