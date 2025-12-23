# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import biom
import pandas as pd
from q2_types.kraken2 import Kraken2OutputDirectoryFormat, Kraken2ReportDirectoryFormat


def _build_taxonomy_to_contig_mapping(kraken2_outputs: Kraken2OutputDirectoryFormat):
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

        df = pd.DataFrame({
            'contig_id': output_df.iloc[:, 1].astype(str),
            'taxon_id': output_df.iloc[:, 2].astype(str)
        })

        # Group by taxon ID and aggregate contig IDs into lists
        grouped = df.groupby('taxon_id')['contig_id'].apply(list).to_dict()

        for taxon_id, contigs in grouped.items():
            if taxon_id not in taxon_to_contigs:
                taxon_to_contigs[taxon_id] = []
            taxon_to_contigs[taxon_id].extend(contigs)

    return taxon_to_contigs


def _calculate_collapsed_counts(contig_ft_collapsed):
    """
    Calculate the number of contigs collapsed into each taxonomy.

    Args:
        contig_ft_collapsed: Collapsed biom.Table

    Returns:
        dict: Mapping of taxonomy ID to count of collapsed contigs
    """
    counts = {}
    for obs_id in contig_ft_collapsed.ids(axis="observation"):
        metadata = contig_ft_collapsed.metadata(id=obs_id, axis="observation")
        counts[obs_id] = len(metadata["collapsed_ids"])
    return counts


def _average_by_count(contig_ft_collapsed):
    """
    Average abundances by dividing by the number of collapsed contigs.

    Args:
        contig_ft_collapsed: Collapsed biom.Table

    Returns:
        biom.Table: Table with averaged abundances
    """
    counts = _calculate_collapsed_counts(contig_ft_collapsed)

    def divide_by_count(data, id_, metadata):
        return data / counts[id_]

    return contig_ft_collapsed.transform(divide_by_count, axis="observation")


def kraken2_to_contig_taxonomy(
    ctx,
    reports,
    outputs,
    table,
    coverage_threshold=0.1,
):
    """
    Map contig IDs to taxonomy strings based on Kraken2 classifications.

    Args:
        reports (Kraken2ReportDirectoryFormat): Directory containing Kraken2
            report files for contigs.
        outputs (Kraken2OutputDirectoryFormat): Directory containing Kraken2
            output files with contig classifications.
        table: Feature table of contig abundances per sample.
        coverage_threshold (float): Minimum percent coverage required to produce
            a feature. Default is 0.1.

    Returns:
        tuple: A tuple containing:
            - Taxonomy artifact: Series mapping contig IDs to taxonomy strings
            - Taxonomy abundance table: Feature table with averaged abundances
    """
    # Step 1: Get taxonomy mapping from Kraken2 reports
    _to_features = ctx.get_action("annotate", "kraken2_to_features")
    _, taxonomy_df = _to_features(reports, coverage_threshold)
    taxonomy_df = taxonomy_df.view(pd.DataFrame)
    outputs = outputs.view(Kraken2OutputDirectoryFormat)
    id_to_taxonomy = taxonomy_df["Taxon"].to_dict()

    # Step 2: Build contig to taxonomy mapping
    contig_taxonomy = _build_taxonomy_to_contig_mapping(outputs)
    taxonomy_series = pd.Series(contig_taxonomy, name="Taxon")
    taxonomy_series.index.name = "Feature ID"

    # Step 3: Collapse feature table by taxonomy
    taxonomy_max_level = taxonomy_series.str.split(";").map(len).max()

    def _collapse(id_, md):
        try:
            tax = taxonomy_series.loc[id_]
        except KeyError:
            tax = "d__Unclassified"
        tax = [x.strip() for x in tax.split(";")]
        if len(tax) < taxonomy_max_level:
            padding = ["__"] * (taxonomy_max_level - len(tax))
            tax.extend(padding)
        return ";".join(tax)

    contig_ft = table.view(biom.Table)
    contig_ft_collapsed = contig_ft.collapse(
        f=_collapse, axis="observation", norm=False
    )

    # Step 4: Average abundances by number of collapsed contigs
    contig_ft_averaged = _average_by_count(contig_ft_collapsed)

    # Step 5: Create and return artifacts
    taxonomy_artifact = ctx.make_artifact("FeatureData[Taxonomy]", taxonomy_series)
    ft_artifact = ctx.make_artifact("FeatureTable[Frequency]", contig_ft_averaged)

    return taxonomy_artifact, ft_artifact

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
    taxonomy = ctx.make_artifact("FeatureData[Taxonomy]", taxonomy)

    outputs = outputs.view(Kraken2OutputDirectoryFormat)
    feature_map = ctx.make_artifact(
        "FeatureMap[TaxonomyToContigs]",
        _build_taxonomy_to_contig_mapping(outputs)
    )

    return feature_map, taxonomy
