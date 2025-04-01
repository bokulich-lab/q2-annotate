# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import itertools

from skbio import TreeNode

from q2_types.kraken2 import (
    Kraken2ReportDirectoryFormat,
    Kraken2OutputDirectoryFormat,
    Kraken2ReportFormat,
    Kraken2OutputFormat,
)


def _merge_kraken2_results(
    reports: list[Kraken2ReportDirectoryFormat],
    outputs: list[Kraken2OutputDirectoryFormat]
) -> tuple[Kraken2ReportDirectoryFormat, Kraken2OutputDirectoryFormat]:
    '''
    Merges kraken2 reports and outputs into a single report and format per
    unique sample id.

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

    mags = False
    if isinstance(next(iter(reports[0].file_dict().values())), dict):
        mags = True

    report_mapping, output_mapping = _condense_formats(reports, outputs, mags)

    def _merge_formats(sample_id, merged_formats, mapping, merger, Format):
        for sample_id in report_mapping:
            if mags:
                for filename, formats in mapping[sample_id].items():
                    merged_format = merger(formats)
                    merged_reports.write_data(
                        merged_format,
                        Format,
                        sample_id=sample_id,
                        mag_id=filename
                    )
            else:
                formats = format[sample_id]
                merged_format = merger(formats)
                merged_formats.write_data(
                    merged_format, Format, sample_id=sample_id
                )

    for sample_id in report_mapping:
        _merge_formats(
            sample_id,
            merged_reports,
            report_mapping,
            _merge_reports,
            Kraken2ReportFormat
        )
        _merge_formats(
            sample_id,
            merged_outputs,
            output_mapping,
            _merge_outputs,
            Kraken2OutputFormat
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
    sample_id -> {filename -> list[format]} for mags.

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
    '''
    chained_reports = itertools.chain(
        *[report.file_dict() for report in reports]
    )
    chained_outputs = itertools.chain(
        *[output.file_dict() for output in outputs]
    )

    def _update_mapping(sample_id, chain, mapping, Format):
        if sample_id not in mapping:
            if mags:
                mapping[sample_id] = {}
                for filename, filepath in chain[sample_id].items():
                    format = Format(filepath, mode='r')
                    mapping[sample_id][filename] = [format]
            else:
                format = Format(chain[sample_id], mode='r')
                mapping[sample_id] = [format]
        else:
            if mags:
                for filename, filepath in chain[sample_id].items():
                    format = Format(filepath, mode='r')
                    mapping[sample_id][filename].append(format)
            else:
                format = Format(chain[sample_id], mode='r')
                mapping[sample_id].append(format)

    report_mapping = {}
    output_mapping = {}

    for sample_id in chained_reports:
        _update_mapping(
            sample_id, chained_reports, report_mapping, Kraken2ReportFormat
        )
        _update_mapping(
            sample_id, chained_outputs, output_mapping, Kraken2OutputFormat
        )

    return report_mapping, output_mapping


def _merge_reports(
    report_mapping: dict[str, Kraken2ReportFormat]
) -> Kraken2ReportDirectoryFormat:
    '''
    '''
    pass


def _merge_outputs(
    output_mapping: dict[str, Kraken2OutputFormat]
) -> Kraken2OutputDirectoryFormat:
    '''
    '''
    pass


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
        for node in first_tree.levelorder():
            match = _find_node(second_tree, node)
            if match is not None:
                match._kraken_data['n_frags_assigned'] += \
                    node._kraken_data['n_frags_assigned']
                match._kraken_data['n_frags_covered'] += \
                    node._kraken_data['n_frags_assigned']
            else:
                parent = _find_node(second_tree, node.parent)
                new_node = TreeNode()
                new_node._kraken_data = node._kraken_data
                new_node._kraken_data['n_frags_covered'] = \
                    new_node._kraken_data['n_frags_assigned']
                parent.append(new_node)
                match = new_node

            # keep `n_frags_covered` up to date
            for ancestor in match.ancestors():
                ancestor._kraken_data['n_frags_covered'] += \
                    node._kraken_data['n_frags_assigned']

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


def _find_node(tree: TreeNode, node: TreeNode) -> TreeNode | None:
    '''
    Searches for `node` in `tree`, returns a match if there is one, otherwise
    None. Matches are determined based on the kraken2-assigned taxon id.

    Parameters
    ----------
    tree : skbio.TreeNode
        The tree in which to search.
    node : skbio.TreeNode
        The node for which to search.

    Returns
    -------
    skbio.TreeNode | None
        The matching node if one was found, otherwise None.

    Raises
    ------
    ValueError
        If more than one matching node is found. Taxon ids should be unique
        within a tree.
    '''
    def _match_by_taxon_id(n):
        return n._kraken_data['taxon_id'] == node._kraken_data['taxon_id']

    matches = list(tree.find_by_func(_match_by_taxon_id))

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
