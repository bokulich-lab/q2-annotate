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

import biom
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.ipc as ipc
import q2templates
from q2_types.kraken2 import Kraken2OutputDirectoryFormat, Kraken2ReportDirectoryFormat

TEMPLATES = resources.files("q2_annotate") / "assets"


def _build_contig_map(kraken2_outputs: Kraken2OutputDirectoryFormat):
    """
    Build a mapping from taxon IDs to lists of contig IDs.

    Args:
        kraken2_outputs: Kraken2OutputDirectoryFormat with classification results

    Returns:
        dict: Mapping of taxon ID to list of contig IDs
    """
    taxon_to_contigs = {}

    for relpath, output_df in kraken2_outputs.outputs.iter_views(pd.DataFrame):
        if output_df.empty:
            continue

        df = pd.DataFrame(
            {
                "contig_id": output_df.iloc[:, 1].astype(str),
                "taxon_id": output_df.iloc[:, 2].astype(str),
            }
        )

        # Group by taxon ID and aggregate contig IDs into lists
        grouped = df.groupby("taxon_id")["contig_id"].apply(list).to_dict()

        for taxon_id, contigs in grouped.items():
            if taxon_id not in taxon_to_contigs:
                taxon_to_contigs[taxon_id] = []
            taxon_to_contigs[taxon_id].extend(contigs)

    return taxon_to_contigs


def _average_by_count(table: biom.Table, contig_map: dict) -> biom.Table:
    """
    Average abundances by dividing by the number of collapsed contigs.

    Args:
        table: Collapsed biom.Table
        contig_map: Mapping of contig IDs to taxonomy IDs

    Returns:
        biom.Table: Table with averaged abundances
    """
    counts = {tax_id: len(contigs) for tax_id, contigs in contig_map.items()}

    def divide_by_count(data, id_, _):
        return data / counts[id_]

    return table.transform(divide_by_count, axis="observation")


def _df_to_arrow_with_arrays(df: pd.DataFrame, filename: str, output_dir: str):
    """
    Writes a pandas DataFrame with array columns to an Arrow IPC file.
    Arrow format supports list/array types natively.

    Args:
        df (pd.DataFrame): The DataFrame to write (may contain list/array columns).
        filename (str): The name of the Arrow file to create.
        output_dir (str): The directory where the 'data' subdirectory
                          will be located or created, and the Arrow file
                          will be saved.
    """
    # Convert DataFrame to Arrow Table
    # PyArrow automatically handles list columns
    table = pa.Table.from_pandas(df)
    fp = os.path.join(output_dir, "data", filename)
    with ipc.RecordBatchFileWriter(fp, table.schema) as writer:
        writer.write_table(table)


def _save_table_data_efficiently(
    table: biom.Table,
    contig_map_rev: dict,
    taxonomy: pd.Series,
    output_dir: str
):
    """
    Save table data as Parquet with arrays of abundances per taxonomy per sample.

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
    table_df = pd.DataFrame({
        "contig_id": obs_id_array,
        "sample": sample_id_array,
        "abundance": abundance_array
    })
    
    # Filter out zero abundances (shouldn't be any in COO, but just in case)
    table_df = table_df[table_df["abundance"] > 0].copy()
    
    # Map contig_id to taxon_id
    table_df["taxon_id"] = table_df["contig_id"].map(lambda x: contig_map_rev.get(x, "0"))
    
    # Map taxon_id to taxon string if taxonomy provided
    if taxonomy is not None:
        table_df["taxon"] = table_df["taxon_id"].astype(str).map(lambda x: taxonomy.get(x, x))
    else:
        table_df["taxon"] = table_df["taxon_id"].astype(str)
    
    # Group by taxon and sample, aggregating abundances into arrays
    grouped = table_df.groupby(["taxon", "sample"])["abundance"].apply(list).reset_index()
    grouped.rename(columns={"abundance": "abundances"}, inplace=True)
    
    # Keep as Python list - pyarrow will convert to Arrow list type automatically
    # Arrow natively supports list/array types
    # grouped["abundances"] is already a list column, which pyarrow handles correctly
    
    # Save to Arrow format (supports arrays and works with Vega-Lite)
    _df_to_arrow_with_arrays(grouped, "abundance_data.arrow", output_dir)


def _visualize_collapsed_contigs(
    output_dir: str,
    table: biom.Table,
    contig_map: dict,
    taxonomy: pd.Series = None,
):
    """
    Generate visualization for original (pre-collapsed) contig abundances.

    Args:
        output_dir: Directory to write visualization files
        table: FeatureTable[Frequency] artifact with contig IDs as observation IDs
        contig_map: FeatureMap[TaxonomyToContigs] artifact
        taxonomy: Optional FeatureData[Taxonomy] artifact for taxonomy strings
    """
    # Extract biom.Table and build reverse mapping
    contig_map_rev = {
        contig: tax_id
        for tax_id, contigs in contig_map.items()
        for contig in contigs
    }
    
    os.makedirs(os.path.join(output_dir, "data"), exist_ok=True)

    # Save table data as Parquet with arrays of abundances per taxonomy per sample
    _save_table_data_efficiently(table, contig_map_rev, taxonomy, output_dir)
    
    # Get unique samples for dropdowns from the table directly
    sample_ids = [str(id_) for id_ in table.ids(axis="sample")]
    samples_list = sorted(sample_ids)

    templates = [
        TEMPLATES / "kraken2_collapse" / "index.html",
    ]

    vega_spec_path = TEMPLATES / "kraken2_collapse" / "vega" / "abundance_histogram.json"
    with open(vega_spec_path, 'r') as f:
        vega_spec = json.load(f)
    
    context = {
        "samples": json.dumps(samples_list),
        "vega_abundance_histogram_spec": json.dumps(vega_spec),
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
        table: Feature table of contig abundances per sample.
        taxonomy: Optional FeatureData[Taxonomy] for taxonomy strings.

    Returns:
        collapse_table: Feature table with contig abundances collapsed per taxonomy.
        visualization: Visualization of abundance distributions.
    """
    contig_map_dict = contig_map.view(dict)
    contig_map_rev = {
        contig: tax_id
        for tax_id, contigs in contig_map_dict.items()
        for contig in contigs
    }

    contig_ft = table.view(biom.Table)
    
    # Generate visualization using original table (before collapsing)
    _visualize = ctx.get_action("annotate", "_visualize_collapsed_contigs")
    (viz,) = _visualize(table, contig_map, taxonomy)
    
    # Collapse the table
    contig_ft_collapsed = contig_ft.collapse(
        f=lambda id_, md: contig_map_rev.get(id_, "0"), axis="observation", norm=False
    )
    contig_ft_averaged = _average_by_count(contig_ft_collapsed, contig_map_dict)
    table = ctx.make_artifact("FeatureTable[Frequency]", contig_ft_averaged)

    return table, viz


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
    feature_map = _build_contig_map(
        outputs.view(Kraken2OutputDirectoryFormat)
    )

    # TODO: we should probably still account for the fact that "0" may not
    #  be present in the feature map sometimes

    if len(taxonomy) < len(feature_map.keys()):
        # we need to account for the taxa with the coverage under the threshold
        missing_taxa = set(feature_map.keys()) - set(taxonomy.index)
        for taxon in missing_taxa:
            feature_map["0"].extend(feature_map[taxon])
            del feature_map[taxon]

    taxonomy = ctx.make_artifact("FeatureData[Taxonomy]", taxonomy)
    feature_map = ctx.make_artifact("FeatureMap[TaxonomyToContigs]", feature_map)

    return feature_map, taxonomy
