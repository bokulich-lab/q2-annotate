# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import base64
import hashlib
import re
from collections import defaultdict
from urllib.parse import quote

import biom
import numpy as np
import pandas as pd
import scipy.sparse as sp

from .types import TFA_TAXONOMY_COLUMNS
from ._mapping import _align_gene_taxonomy


def _gene_taxonomy_id(taxon_id: str, gene_id: str) -> str:
    """Join readable IDs, escaping the separator and literal percent signs."""
    return f"{quote(str(taxon_id), safe='')}|{quote(str(gene_id), safe='')}"


def _map_contigs_to_taxa(
    taxonomy: pd.DataFrame, taxon_to_contig_map: dict
) -> dict[str, str]:
    """Invert assignments for taxa with taxonomy, rejecting ambiguous contigs."""
    contig_to_taxon = {}
    for taxon_id, contigs in taxon_to_contig_map.items():
        if taxon_id not in taxonomy.index:
            continue
        for contig_id in contigs:
            if contig_id in contig_to_taxon and contig_to_taxon[contig_id] != taxon_id:
                raise ValueError(f"Contig {contig_id!r} is assigned to multiple taxa.")
            contig_to_taxon[contig_id] = taxon_id
    return contig_to_taxon


def _align_contig_matrices(
    contig_abundance: biom.Table,
    feature_inventory: biom.Table,
    contig_to_taxon: dict[str, str],
) -> tuple[list[str], sp.csr_matrix, sp.csr_matrix]:
    """Align shared contigs as rows of abundance and gene inventory matrices."""
    abundance_ids = list(contig_abundance.ids(axis="observation"))
    inventory_ids = list(feature_inventory.ids(axis="sample"))
    contig_ids = sorted(set(abundance_ids) & set(inventory_ids) & set(contig_to_taxon))
    abundance_index = {identifier: row for row, identifier in enumerate(abundance_ids)}
    inventory_index = {
        identifier: column for column, identifier in enumerate(inventory_ids)
    }

    # Both matrices have the same contig rows; columns are samples and genes.
    abundances = contig_abundance.matrix_data.tocsr()[
        [abundance_index[identifier] for identifier in contig_ids], :
    ]
    inventory = feature_inventory.matrix_data.tocsc()[
        :, [inventory_index[identifier] for identifier in contig_ids]
    ].T.tocsr()
    inventory.eliminate_zeros()
    return contig_ids, abundances, inventory


def _calculate_taxon_gene_loads(
    abundances: sp.csr_matrix,
    inventory: sp.csr_matrix,
    contig_ids: list[str],
    contig_to_taxon: dict[str, str],
    gene_ids: list[str],
) -> tuple[sp.csr_matrix, list[tuple[str, str]]]:
    """Calculate pair-by-sample loads, omitting pairs with no sample load."""
    contig_rows_by_taxon = defaultdict(list)
    for row, contig_id in enumerate(contig_ids):
        contig_rows_by_taxon[contig_to_taxon[contig_id]].append(row)

    blocks, pairs = [], []
    for taxon_id, contig_rows in sorted(contig_rows_by_taxon.items()):
        taxon_inventory = inventory[contig_rows, :]
        gene_columns = np.unique(taxon_inventory.indices)
        if not len(gene_columns):
            continue
        # Sparse multiplication avoids a dense taxon x gene x sample cube.
        loads = (
            taxon_inventory[:, gene_columns].T @ abundances[contig_rows, :]
        ).tocsr()
        loads.eliminate_zeros()
        nonzero_rows = np.flatnonzero(np.diff(loads.indptr))
        if len(nonzero_rows):
            blocks.append(loads[nonzero_rows, :])
            pairs.extend(
                (taxon_id, gene_ids[gene_columns[row]]) for row in nonzero_rows
            )

    matrix = (
        sp.vstack(blocks, format="csr")
        if blocks
        else sp.csr_matrix((0, abundances.shape[1]))
    )
    if not np.all(np.isfinite(matrix.data)):
        raise ValueError("Estimate loads must contain non-negative floats.")
    return matrix, pairs


def _estimate_tfa(
    contig_abundance: biom.Table,
    feature_inventory: biom.Table,
    taxonomy: pd.DataFrame,
    taxon_to_contig_map: dict,
) -> (biom.Table, pd.DataFrame):
    """Estimate sparse, sample-resolved taxon/function loads."""
    contig_to_taxon = _map_contigs_to_taxa(taxonomy, taxon_to_contig_map)
    contig_ids, abundances, inventory = _align_contig_matrices(
        contig_abundance, feature_inventory, contig_to_taxon
    )
    loads, pairs = _calculate_taxon_gene_loads(
        abundances,
        inventory,
        contig_ids,
        contig_to_taxon,
        list(feature_inventory.ids(axis="observation")),
    )

    if not pairs:
        raise ValueError("No taxon/gene pairs have nonzero loads in the input samples.")

    pair_ids = [_gene_taxonomy_id(taxon_id, gene_id) for taxon_id, gene_id in pairs]
    gene_taxonomy = pd.DataFrame(
        [
            [taxonomy.at[taxon_id, "Taxon"], taxon_id, gene_id]
            for taxon_id, gene_id in pairs
        ],
        index=pd.Index(pair_ids, name="Feature ID"),
        columns=TFA_TAXONOMY_COLUMNS,
    )
    return (
        biom.Table(
            loads,
            observation_ids=pair_ids,
            sample_ids=contig_abundance.ids(axis="sample"),
        ),
        gene_taxonomy,
    )


def _collection_key(identifier: str) -> str:
    """Make a deterministic, filesystem-safe QIIME 2 collection key."""
    if (
        len(identifier) <= 100
        and re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9+_.~-]*", identifier)
        and not identifier.startswith(("id_", "hash_"))
    ):
        return identifier
    encoded = base64.urlsafe_b64encode(identifier.encode("utf-8"))
    key = "id_" + encoded.decode("ascii").rstrip("=")
    if len(key) <= 200:
        return key
    return "hash_" + hashlib.sha256(identifier.encode("utf-8")).hexdigest()


def _select_tfa(
    feature_load: biom.Table,
    gene_taxonomy: pd.DataFrame,
    ids: list[str],
    dimension: int,
) -> dict[str, biom.Table]:
    label = "taxon" if dimension == 0 else "function"
    selected = {_id: ([], []) for _id in ids}
    if not selected:
        raise ValueError(f"Select at least one {label} ID.")

    mapping = _align_gene_taxonomy(feature_load, gene_taxonomy)
    pairs = mapping[["Taxon ID", "Gene ID"]].itertuples(index=False, name=None)
    for index, pair in enumerate(pairs):
        if pair[dimension] in selected:
            rows, other_ids = selected[pair[dimension]]
            rows.append(index)
            other_ids.append(pair[1 - dimension])

    missing = [_id for _id, (rows, _) in selected.items() if not rows]
    if missing:
        raise ValueError(f"No TFA rows found for {label} IDs: {missing!r}.")

    matrix = feature_load.matrix_data.tocsr()
    sample_ids = feature_load.ids(axis="sample")
    return {
        _collection_key(_id): biom.Table(
            matrix[rows, :], observation_ids=other_ids, sample_ids=sample_ids
        )
        for _id, (rows, other_ids) in selected.items()
    }


def tfa_by_taxon(
    feature_load: biom.Table, gene_taxonomy: pd.DataFrame, taxon_ids: list[str]
) -> biom.Table:
    """Return function-by-sample tables for the requested taxon IDs."""
    # QIIME 2 uses the collection member's view type as the return annotation.
    return _select_tfa(feature_load, gene_taxonomy, taxon_ids, 0)


def tfa_by_function(
    feature_load: biom.Table, gene_taxonomy: pd.DataFrame, function_ids: list[str]
) -> biom.Table:
    """Return taxon-by-sample tables for the requested function IDs."""
    return _select_tfa(feature_load, gene_taxonomy, function_ids, 1)


def tfa_collapse_samples(
    feature_load: biom.Table, gene_taxonomy: pd.DataFrame, aggregation: str = "sum"
) -> biom.Table:
    """Collapse samples into one sparse taxon-by-function table."""
    if aggregation not in {"sum", "mean", "median"}:
        raise ValueError(f"Unsupported sample aggregation: {aggregation!r}.")

    matrix = feature_load.matrix_data.tocsr()
    sample_count = matrix.shape[1]
    if sample_count == 0 and aggregation != "sum":
        raise ValueError(f"Cannot calculate {aggregation} without samples.")

    mapping = _align_gene_taxonomy(feature_load, gene_taxonomy)
    pairs = list(mapping[["Taxon ID", "Gene ID"]].itertuples(index=False, name=None))
    taxon_ids = list(dict.fromkeys(taxon_id for taxon_id, _ in pairs))
    function_ids = list(dict.fromkeys(function_id for _, function_id in pairs))
    taxon_indices = {taxon_id: index for index, taxon_id in enumerate(taxon_ids)}
    function_indices = {
        function_id: index for index, function_id in enumerate(function_ids)
    }

    if aggregation == "median":
        collapsed = np.zeros(len(pairs), dtype=float)
        for row in range(len(pairs)):
            data = matrix.data[matrix.indptr[row] : matrix.indptr[row + 1]]
            nonzero = np.sort(data[data != 0])
            zero_count = sample_count - len(nonzero)

            def at_rank(rank):
                return 0.0 if rank < zero_count else nonzero[rank - zero_count]

            collapsed[row] = (
                at_rank((sample_count - 1) // 2) + at_rank(sample_count // 2)
            ) / 2
    else:
        collapsed = np.asarray(matrix.sum(axis=1)).ravel()
        if aggregation == "mean":
            collapsed = collapsed / sample_count

    present = np.flatnonzero(collapsed)
    result = sp.coo_matrix(
        (
            collapsed[present],
            (
                [taxon_indices[pairs[row][0]] for row in present],
                [function_indices[pairs[row][1]] for row in present],
            ),
        ),
        shape=(len(taxon_ids), len(function_ids)),
    ).tocsr()
    return biom.Table(
        result,
        observation_ids=taxon_ids,
        sample_ids=function_ids,
    )


def estimate_tfa(
    ctx,
    contig_abundance,
    feature_inventory,
    taxonomy,
    taxon_to_contig_map,
    metadata=None,
):
    """Estimate sample-resolved TFA and visualize it with optional metadata."""
    tfa, gene_taxonomy = _estimate_tfa(
        contig_abundance=contig_abundance.view(biom.Table),
        feature_inventory=feature_inventory.view(biom.Table),
        taxonomy=taxonomy.view(pd.DataFrame),
        taxon_to_contig_map=taxon_to_contig_map.view(dict),
    )

    tfa = ctx.make_artifact("FeatureTable[Frequency % Properties('tfa')]", tfa)
    gene_taxonomy = ctx.make_artifact(
        "FeatureData[Taxonomy % Properties('tfa')]", gene_taxonomy
    )

    visualize = ctx.get_action("annotate", "_visualize_tfa")
    (visualization,) = visualize(
        feature_load=tfa, gene_taxonomy=gene_taxonomy, metadata=metadata
    )
    return tfa, gene_taxonomy, visualization
