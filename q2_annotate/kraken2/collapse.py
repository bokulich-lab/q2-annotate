# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json
import os
import shutil
from importlib import resources
from typing import Optional

import biom
import numpy as np
import pandas as pd
import q2templates
import scipy.sparse as sp
from q2_types.kraken2 import Kraken2OutputDirectoryFormat

TEMPLATES = resources.files("q2_annotate") / "assets"


def _build_contig_map(kraken2_outputs: Kraken2OutputDirectoryFormat) -> dict:
    """
    Build a mapping from taxon IDs to lists of contig IDs.

    Args:
        kraken2_outputs: Kraken2OutputDirectoryFormat with classification results

    Returns:
        dict: Mapping of taxon ID to a list of contig IDs
    """
    taxon_to_contigs = {}

    for relpath, output_df in kraken2_outputs.outputs.iter_views(pd.DataFrame):
        contig_ids = output_df.iloc[:, 1].astype(str)
        taxon_ids = output_df.iloc[:, 2].astype(str)

        # Group by taxon ID and aggregate contig IDs into lists
        grouped = contig_ids.groupby(taxon_ids).apply(list)

        for taxon_id, contigs in grouped.items():
            if taxon_id not in taxon_to_contigs:
                taxon_to_contigs[taxon_id] = []
            taxon_to_contigs[taxon_id].extend(contigs)

    return taxon_to_contigs


def _average_by_count(
    collapsed_table: biom.Table, original_table: biom.Table, contig_map_rev: dict
) -> biom.Table:
    """
    Average abundances by dividing by the number of contigs per sample per taxon.

    Args:
        collapsed_table: Collapsed biom.Table with taxonomy IDs as observations
        original_table: Original biom.Table with contig IDs as observations
        contig_map_rev: Reverse mapping of contig IDs to taxonomy IDs

    Returns:
        biom.Table: Table with averaged abundances
    """

    # We want: for each (taxon, sample), divide the collapsed abundance by the
    # number of *contigs present* (abundance > 0) for that taxon in that sample.

    sample_ids = list(collapsed_table.ids(axis="sample"))
    taxon_ids = list(collapsed_table.ids(axis="observation"))

    taxon_index = {str(tid): i for i, tid in enumerate(taxon_ids)}

    # original table as COO so we can remap row indices efficiently
    orig = original_table.matrix_data.tocoo()

    # keep only present contigs (abundance > 0). In COO, data is non-zero by
    # definition, but we keep the explicit check for safety.
    present_mask = orig.data > 0
    contig_rows = orig.row[present_mask]
    sample_cols = orig.col[present_mask]

    # map each contig row index to a taxon row index
    obs_ids = np.array(list(original_table.ids(axis="observation")), dtype=object)
    contig_to_taxon = np.array(
        [str(contig_map_rev.get(str(cid), "0")) for cid in obs_ids], dtype=object
    )

    mapped_taxon_rows = np.fromiter(
        (taxon_index.get(tid, -1) for tid in contig_to_taxon[contig_rows]),
        dtype=np.int64,
        count=len(contig_rows),
    )

    # drop any contigs whose taxon isn't in the collapsed table (shouldn't happen)
    valid = mapped_taxon_rows >= 0
    mapped_taxon_rows = mapped_taxon_rows[valid]
    sample_cols = sample_cols[valid]

    # build sparse contig-counts per (taxon, sample)
    ones = np.ones(mapped_taxon_rows.shape[0], dtype=np.int32)
    count_matrix = sp.coo_matrix(
        (ones, (mapped_taxon_rows, sample_cols)),
        shape=(len(taxon_ids), len(sample_ids)),
    ).tocsr()

    # divide only where collapsed table has data (stay sparse)
    # if we divided the two matrices directly we would get a dense matrix
    collapsed = collapsed_table.matrix_data.tocoo()
    denom = count_matrix[collapsed.row, collapsed.col].A1

    with np.errstate(divide="ignore", invalid="ignore"):
        new_data = np.divide(
            collapsed.data,
            denom,
            out=np.zeros_like(collapsed.data, dtype=float),
            where=denom > 0,
        )

    averaged = sp.coo_matrix(
        (new_data, (collapsed.row, collapsed.col)),
        shape=collapsed_table.matrix_data.shape,
    ).tocsr()

    return biom.Table(averaged, taxon_ids, sample_ids)


def _df_to_json_per_sample(df: pd.DataFrame, output_dir: str):
    """
    Writes a pandas DataFrame to separate JSON files per sample.
    Each JSON contains an array of objects with taxon and abundance fields.

    Args:
        df (pd.DataFrame): The DataFrame to write (with taxon, sample,
                           abundances columns).
        output_dir (str): The directory where the 'data' subdirectory
                          will be located or created, and the JSON files
                          will be saved.
    """
    data_dir = os.path.join(output_dir, "data")

    # Group by sample and save each to a separate JSON
    for sample, group in df.groupby("sample"):
        # Explode the abundances column to get one row per abundance value
        sample_df = group[["taxon", "abundances"]].copy()
        sample_df = sample_df.explode("abundances")
        sample_df.columns = ["taxon", "abundance"]

        # Save to JSON file named after the sample
        json_filename = f"{sample}.json"
        json_path = os.path.join(data_dir, json_filename)
        sample_df.to_json(json_path, orient="records", indent=2)


def _table_to_json(
    table: biom.Table,
    contig_map_rev: dict,
    taxonomy: Optional[pd.Series],
    output_dir: str,
):
    """
    Save table data as JSON files with one file per sample.

    Args:
        table: biom.Table to save
        contig_map_rev: Dictionary mapping contig IDs to taxonomy IDs
        taxonomy: Optional pd.Series mapping taxonomy IDs to taxonomy strings
        output_dir: Directory where the 'data' subdirectory will be located
    """
    # Access sparse matrix directly (CSR format)
    matrix = table.matrix_data
    obs_ids = table.ids(axis="observation")
    sample_ids = table.ids(axis="sample")

    # Convert sparse matrix to COO format for efficient processing
    coo = matrix.tocoo()

    # Build arrays efficiently - use numpy array indexing for speed
    obs_ids_arr = np.array(obs_ids, dtype=object)
    sample_ids_arr = np.array(sample_ids, dtype=object)

    obs_id_array = obs_ids_arr[coo.row].astype(str)
    sample_id_array = sample_ids_arr[coo.col].astype(str)
    abundance_array = coo.data

    # Create initial DataFrame with contig_id, sample, abundance
    table_df = pd.DataFrame(
        {
            "contig_id": obs_id_array,
            "sample": sample_id_array,
            "abundance": abundance_array,
        }
    )

    # Filter out zero abundances (shouldn't be any in COO, but just in case)
    table_df = table_df[table_df["abundance"] > 0].copy()

    # Map contig_id to taxon_id
    table_df["taxon_id"] = table_df["contig_id"].map(
        lambda x: contig_map_rev.get(x, "0")
    )

    # Map taxon_id to taxon string if taxonomy provided
    if taxonomy is not None:
        table_df["taxon"] = (
            table_df["taxon_id"].astype(str).map(lambda x: taxonomy.get(x, x))
        )
    else:
        table_df["taxon"] = table_df["taxon_id"].astype(str)

    # Group by taxon and sample, aggregating abundances into arrays
    grouped = (
        table_df.groupby(["taxon", "sample"])["abundance"].apply(list).reset_index()
    )
    grouped.rename(columns={"abundance": "abundances"}, inplace=True)

    # Save to JSON files per sample
    _df_to_json_per_sample(grouped, output_dir)


def _extract_mean_abundances(
    collapsed_table: biom.Table, taxonomy: pd.Series = None
) -> dict:
    """
    Extract mean abundances per taxon per sample from the collapsed table.

    Args:
        collapsed_table: feature table with contig IDs collapsed
            to taxonomy IDs
        taxonomy: Optional taxonomy for taxonomy strings

    Returns:
        dict: Mapping of taxon -> sample -> mean abundance
    """
    collapsed_df = collapsed_table.to_dataframe(dense=True)

    # Map taxon IDs to taxonomy strings if taxonomy provided
    if taxonomy is not None:
        collapsed_df.index = collapsed_df.index.map(
            lambda x: taxonomy.get(str(x), str(x))
        )

    # Calculate mean abundances per taxon per sample
    mean_abundances = {}
    for taxon in collapsed_df.index:
        taxon_str = str(taxon)
        mean_abundances[taxon_str] = {}

        for sample in collapsed_df.columns:
            sample_str = str(sample)
            value = collapsed_df.loc[taxon, sample]
            if pd.notna(value) and value > 0:
                mean_abundances[taxon_str][sample_str] = float(value)

    return mean_abundances


def _visualize_collapsed_contigs(
    output_dir: str,
    table: biom.Table,
    collapsed_table: biom.Table,
    contig_map: dict,
    taxonomy: pd.Series = None,
):
    """
    Generate visualization for original (pre-collapsed) contig abundances.

    Args:
        output_dir: Directory to write visualization files
        table: FeatureTable[Frequency] artifact with contig IDs as observation IDs
        collapsed_table: FeatureTable[Frequency] artifact with contig IDs
            collapsed to taxonomy IDs
        contig_map: FeatureMap[TaxonomyToContigs] artifact
        taxonomy: Optional FeatureData[Taxonomy] artifact for taxonomy strings
    """
    # Extract biom.Table and build reverse mapping
    contig_map_rev = {
        contig: tax_id for tax_id, contigs in contig_map.items() for contig in contigs
    }

    os.makedirs(os.path.join(output_dir, "data"), exist_ok=True)

    # Save table data as JSON files with one file per sample
    _table_to_json(table, contig_map_rev, taxonomy, output_dir)

    # Get unique samples for dropdowns from the table directly
    sample_ids = [str(id_) for id_ in table.ids(axis="sample")]
    samples_list = sorted(sample_ids)

    # Extract mean abundances from collapsed_table for efficient sorting
    mean_abundances = _extract_mean_abundances(collapsed_table, taxonomy)

    templates = [
        TEMPLATES / "kraken2_collapse" / "index.html",
    ]

    vega_spec_path = (
        TEMPLATES / "kraken2_collapse" / "vega" / "abundance_histogram.json"
    )
    with open(vega_spec_path, "r") as f:
        vega_spec = json.load(f)

    context = {
        "samples": json.dumps(samples_list),
        "vega_abundance_histogram_spec": json.dumps(vega_spec),
        "mean_abundances": json.dumps(mean_abundances),
    }

    # Copy JS/CSS files
    for d in ("js", "css"):
        src_dir = TEMPLATES / "kraken2_collapse" / d
        dst_dir = os.path.join(output_dir, d)
        if src_dir.exists():
            shutil.copytree(src_dir, dst_dir, dirs_exist_ok=True)

    q2templates.render(templates, output_dir, context=context)


def collapse_contigs(
    ctx,
    contig_map,
    table,
    taxonomy=None,
):
    """
    Map contig IDs to taxonomy strings based on Kraken2 classifications.

    Args:
        contig_map: Taxon-to-contigs mapping.
        table: Feature table of contig abundances per sample.
        taxonomy: Optional taxonomy for taxonomy strings.

    Returns:
        collapsed_table: Feature table with averaged contig abundances
            collapsed per taxonomy.
        visualization: Visualization of abundance distributions.
    """
    contig_map_dict = contig_map.view(dict)
    contig_map_rev = {
        contig: tax_id
        for tax_id, contigs in contig_map_dict.items()
        for contig in contigs
    }

    contig_ft = table.view(biom.Table)

    # Collapse the table
    contig_ft_collapsed = contig_ft.collapse(
        f=lambda id_, md: contig_map_rev.get(id_, "0"), axis="observation", norm=False
    )
    # Average by contig counts per sample
    contig_ft_averaged = _average_by_count(
        contig_ft_collapsed, contig_ft, contig_map_rev
    )
    collapsed_table = ctx.make_artifact("FeatureTable[Frequency]", contig_ft_averaged)

    # Generate visualization using original table (before collapsing)
    _visualize = ctx.get_action("annotate", "_visualize_collapsed_contigs")
    (viz,) = _visualize(table, collapsed_table, contig_map, taxonomy)

    return collapsed_table, viz


def map_taxonomy_to_contigs(
    ctx,
    reports,
    outputs,
    coverage_threshold=0.1,
):
    """
    Map contig IDs to taxonomy strings based on Kraken2 classifications.

    Args:
        reports (Kraken2ReportDirectoryFormat): Directory containing Kraken2
            report files for contigs.
        outputs (Kraken2OutputDirectoryFormat): Directory containing Kraken2
            output files with contig classifications.
        coverage_threshold (float): Minimum percent coverage required to produce
            a feature. Default is 0.1.

    Returns:
        dict: Mapping of contig IDs to taxonomy strings.
    """

    _to_features = ctx.get_action("annotate", "kraken2_to_features")
    _, taxonomy = _to_features(reports, coverage_threshold)

    # add unclassified
    taxonomy = taxonomy.view(pd.Series)
    taxonomy["0"] = "d__Unclassified"
    feature_map = _build_contig_map(outputs.view(Kraken2OutputDirectoryFormat))

    if len(taxonomy) < len(feature_map.keys()):
        # we need to account for the taxa with the coverage under the threshold
        missing_taxa = set(feature_map.keys()) - set(taxonomy.index)
        for taxon in missing_taxa:
            feature_map["0"].extend(feature_map[taxon])
            del feature_map[taxon]

    taxonomy = ctx.make_artifact("FeatureData[Taxonomy]", taxonomy)
    feature_map = ctx.make_artifact("FeatureMap[TaxonomyToContigs]", feature_map)

    return feature_map, taxonomy
