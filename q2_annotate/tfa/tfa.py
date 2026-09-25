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

import biom
import numpy as np
import pandas as pd
import scipy.sparse as sp

from .types._format import _pair_id, _split_pair_id


def estimate_tfa(
    abundance_matrix: biom.Table,
    feature_inventory: biom.Table,
    taxonomy: pd.DataFrame,
    taxon_to_contig_map: dict,
) -> biom.Table:
    """Estimate sparse, sample-resolved taxon/function loads."""
    # Build a mapping from contig_id to taxon_id using taxon_to_contig_map and taxonomy
    contig_to_taxon_id = {}
    for taxon_id, contigs in taxon_to_contig_map.items():
        if taxon_id in taxonomy.index:
            taxon_name = taxonomy.loc[taxon_id, "Taxon"]
            if pd.notnull(taxon_name):
                for cid in contigs:
                    if (
                        cid in contig_to_taxon_id
                        and contig_to_taxon_id[cid] != taxon_id
                    ):
                        raise ValueError(
                            f"Contig {cid!r} is assigned to multiple taxa."
                        )
                    contig_to_taxon_id[cid] = taxon_id

    # Find the intersection of contig IDs across all inputs
    common_contigs = (
        set(abundance_matrix.ids(axis="observation"))
        & set(feature_inventory.ids(axis="sample"))
        & set(contig_to_taxon_id.keys())
    )

    feature_ids = feature_inventory.ids(axis="observation")
    sample_ids = abundance_matrix.ids(axis="sample")

    if not common_contigs:
        return biom.Table(
            sp.csr_matrix((0, len(sample_ids))),
            observation_ids=[],
            sample_ids=sample_ids,
        )

    # Sort common contigs to ensure deterministic ordering
    common_contigs = sorted(list(common_contigs))

    # Slice and align the abundance matrix (observations x samples)
    abund_contig_index = {
        cid: idx for idx, cid in enumerate(abundance_matrix.ids(axis="observation"))
    }
    abund_indices = [abund_contig_index[cid] for cid in common_contigs]
    abundance_matrix_sparse = abundance_matrix.matrix_data.tocsr()

    A = abundance_matrix_sparse[abund_indices, :]

    # Slice and align the feature inventory matrix. Its observations are
    # functions and its samples are contigs.
    feat_contig_index = {
        cid: idx for idx, cid in enumerate(feature_inventory.ids(axis="sample"))
    }
    feat_indices = [feat_contig_index[cid] for cid in common_contigs]
    feature_inventory_sparse = feature_inventory.matrix_data.tocsc()

    inventory = feature_inventory_sparse[:, feat_indices].T.tocsr()
    inventory.eliminate_zeros()

    # For each taxon, multiply its sparse function-by-contig inventory by
    # contig-by-sample abundances. Only observed taxon/function pairs become
    # rows; no dense taxon x function x sample cube is constructed.
    contig_rows_by_taxon = defaultdict(list)
    for row, contig_id in enumerate(common_contigs):
        contig_rows_by_taxon[contig_to_taxon_id[contig_id]].append(row)
    blocks = []
    pair_ids = []
    for taxon_id, contig_rows in sorted(contig_rows_by_taxon.items()):
        taxon_inventory = inventory[contig_rows, :]
        feature_columns = np.unique(taxon_inventory.indices)
        if not len(feature_columns):
            continue
        loads = (taxon_inventory[:, feature_columns].T @ A[contig_rows, :]).tocsr()
        loads.eliminate_zeros()
        blocks.append(loads)
        pair_ids.extend(
            _pair_id(taxon_id, feature_ids[column]) for column in feature_columns
        )

    result_matrix = (
        sp.vstack(blocks, format="csr")
        if blocks
        else sp.csr_matrix((0, len(sample_ids)))
    )
    return biom.Table(
        result_matrix,
        observation_ids=pair_ids,
        sample_ids=sample_ids,
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
    feature_load: biom.Table, identifiers: list[str], dimension: int
) -> dict[str, biom.Table]:
    label = "taxon" if dimension == 0 else "function"
    selected = {identifier: ([], []) for identifier in identifiers}
    if not selected:
        raise ValueError(f"Select at least one {label} ID.")

    for index, observation_id in enumerate(feature_load.ids(axis="observation")):
        pair = _split_pair_id(observation_id)
        if pair[dimension] in selected:
            rows, other_ids = selected[pair[dimension]]
            rows.append(index)
            other_ids.append(pair[1 - dimension])

    missing = [identifier for identifier, (rows, _) in selected.items() if not rows]
    if missing:
        raise ValueError(f"No TFA rows found for {label} IDs: {missing!r}.")

    matrix = feature_load.matrix_data.tocsr()
    sample_ids = feature_load.ids(axis="sample")
    return {
        _collection_key(identifier): biom.Table(
            matrix[rows, :], observation_ids=other_ids, sample_ids=sample_ids
        )
        for identifier, (rows, other_ids) in selected.items()
    }


def tfa_by_taxon(feature_load: biom.Table, taxon_ids: list[str]) -> biom.Table:
    """Return function-by-sample tables for the requested taxon IDs."""
    # QIIME 2 uses the collection member's view type as the return annotation.
    return _select_tfa(feature_load, taxon_ids, 0)


def tfa_by_function(feature_load: biom.Table, function_ids: list[str]) -> biom.Table:
    """Return taxon-by-sample tables for the requested function IDs."""
    return _select_tfa(feature_load, function_ids, 1)
