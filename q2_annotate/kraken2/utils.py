# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from itertools import takewhile, dropwhile
from re import sub
from typing import List

from q2_annotate._utils import _construct_param


def _process_kraken2_arg(arg_key, arg_val):
    """Creates a list with argument and its value to be consumed by Kraken 2.

    Argument names will be converted to command line parameters by
    appending a '--' prefix and replacing all '_' with '-',
    e.g.: 'some_parameter_x' -> '--some-parameter-x'.

    Args:
        arg_key (str): Argument name.
        arg_val: Argument value.

    Returns:
        [converted_arg, arg_value]: List containing a prepared command line
            parameter and, optionally, its value.
    """
    if isinstance(arg_val, bool) and arg_val:
        return [_construct_param(arg_key)]
    elif not isinstance(arg_val, list):
        return [_construct_param(arg_key), str(arg_val)]
    else:
        raise NotImplementedError(
            f'Parsing arguments of type "{type(arg_val)}" is not supported.'
        )


def _find_lca(taxa):
    """Find least common ancestor between two semicolon-delimited strings.

    This code was copied over from https://github.com/bokulich-lab/RESCRIPt.
    """
    # LCA ends where zipped taxonomy strings no longer converge to len == 1
    taxa_comparison = [set(rank) - {None} for rank in zip(*taxa)]
    return (rank.pop() for rank in takewhile(lambda x: len(x) == 1, taxa_comparison))


# def _find_lca_majority(taxa):
#     """Find least common ancestor by majority.
#
#     This code was copied over and adapted from
#      https://github.com/bokulich-lab/RESCRIPt.
#     """
#     # collapse and count unique labels at each rank
#     # yields list of ('labels', counts) sorted by most to least abundant
#     taxa = [[t for t in r if t not in ['', ]] for r in zip_longest(*taxa)]
#     taxa_comparison = [Counter(t).most_common() for t in taxa]
#     # return majority wherever a clear majority is found
#     # terminate when no majority is found, that's your LCA
#     # propagate empty ranks that maintain majority/consensus by inserting
#     # '' to preserve empty ranks in the original taxonomies (e.g., when
#     # using a rank handle unannotated levels will be removed, but as long
#     # as these pass muster we want to preserve those labels)
#     return [rank[0][0] if rank else '' for rank in takewhile(
#         lambda x: len(x) < 2 or x[0][1] > x[1][1], taxa_comparison)]


def _taxon_to_list(taxon, rank_handle):
    """Split taxonomy string into list of taxonomic labels"""
    return [sub(rank_handle, "", t.strip()) for t in taxon.split(";")]


def _join_ranks(taxonomy: List[str], ranks: List[str]) -> str:
    """Join ranks into a single string annotated with provided rank labels.

    All the terminal None values will be removed.
    """
    taxonomy = list(dropwhile(lambda x: x is None, taxonomy[::-1]))
    taxonomy = ["" if x is None else x for x in taxonomy[::-1]]
    taxonomy = ";".join(["".join(t) for t in zip(ranks, taxonomy)])
    return taxonomy
