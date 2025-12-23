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


def collapse_contigs(
    ctx,
    contig_map,
    table,
):
    """
    Map contig IDs to taxonomy strings based on Kraken2 classifications.

    Args:
        table: Feature table of contig abundances per sample.

    Returns:
        collapse_table: Feature table with contig abundances collapsed per taxonomy.
    """
    contig_map_dict = contig_map.view(dict)
    contig_map_rev = {
        contig: tax_id
        for tax_id, contigs in contig_map_dict.items()
        for contig in contigs
    }

    contig_ft = table.view(biom.Table)
    contig_ft_collapsed = contig_ft.collapse(
        f=lambda id_, md: contig_map_rev.get(id_, "0"), axis="observation", norm=False
    )
    contig_ft_averaged = _average_by_count(contig_ft_collapsed, contig_map_dict)
    table = ctx.make_artifact("FeatureTable[Frequency]", contig_ft_averaged)

    return table


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
        "FeatureMap[TaxonomyToContigs]", _build_taxonomy_to_contig_mapping(outputs)
    )

    return feature_map, taxonomy
