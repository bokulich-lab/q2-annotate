# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import biom
import glob
import os
from collections import deque
from typing import List

import numpy as np

from q2_annotate.kraken2.utils import _find_lca, _taxon_to_list, _join_ranks
from q2_types.kraken2 import Kraken2ReportDirectoryFormat, Kraken2OutputDirectoryFormat

import pandas as pd
import skbio

RANKS = "dkpcofgs"


def _fill_unclassified(row):
    if row.isnull().all():
        row[0] = "Unclassified"
    return row


def _find_lcas(taxa_list: List[pd.DataFrame], mode: str):
    """Find the least common ancestor in every DataFrame of taxa.

    Args:
        taxa_list (List[pd.DataFrame]): A list of taxonomy DataFrames.
        mode (str): The mode used to determine the least common ancestor.

    Returns:
        pd.DataFrame: A DataFrame containing the LCA of each feature (MAG).
    """
    methods = {
        "lca": _find_lca,
        # 'super': _find_super_lca,
        # 'majority': _find_lca_majority
    }
    func = methods[mode]
    taxa = pd.concat(taxa_list)

    # Convert taxonomies to list; optionally remove rank handle
    taxa["Taxon"] = taxa["Taxon"].apply(
        lambda x: _taxon_to_list(x, rank_handle=f"^[{RANKS[:-1]}]__|s1?__")
    )

    # Find LCA for every MAG
    results = {}
    for mag_id in taxa["mag_id"].unique():
        data = taxa[taxa["mag_id"] == mag_id]["Taxon"]
        result = func(data)
        results[mag_id] = result

    results = pd.DataFrame.from_dict(results, orient="index")
    results = results.apply(_fill_unclassified, axis=1)
    results = results.apply(lambda x: x.tolist(), axis=1).to_frame()
    results.columns = ["Taxon"]

    # Join ranks
    ranks = [*[f"{r}__" for r in RANKS], "ssp__"]
    results["Taxon"] = results["Taxon"].apply(lambda x: _join_ranks(x, ranks))

    results.index.name = "Feature ID"
    return results


def _add_unclassified_mags(
    table: pd.DataFrame,
    taxonomy: pd.DataFrame,
    reports: Kraken2ReportDirectoryFormat,
    coverage_threshold: float,
) -> pd.DataFrame:
    """
    Identify and process MAGs that are entirely unclassified based on
    Kraken 2 reports.

    Args:
        table (pd.DataFrame): A DataFrame representing the feature table with
                              MAGs as rows and assignments as columns.
        taxonomy (pd.DataFrame): A DataFrame containing taxonomy assignments
                                 for the features.
        reports (Kraken2ReportDirectoryFormat): Directory format containing
                                            Kraken 2 report files for each MAG.
        coverage_threshold (float): Minimum coverage threshold used to evaluate
                                    and classify MAGs from the Kraken 2 reports.

    Returns:
        pd.DataFrame: The updated taxonomy DataFrame with unclassified MAGs added.

    Raises:
        ValueError: If the unclassified fraction for a MAG is not 100%,
        indicating incomplete classification or inconsistent inputs.
    """
    samples = [
        os.path.basename(f).replace(".report.txt", "")
        for f in glob.glob(os.path.join(reports.path, "*.report.txt"))
    ]
    unclassified = set(samples) - set(table.index)

    # For each MAG missing from the table we will look at the original
    # Kraken 2 report and check whether the unclassified counts add up
    # by summing up the following fractions:
    # - line 1: fraction of unclassified
    # - line 2 (if present): fraction of classified as root (we treat
    #                        those as practically unclassified); if the
    #                        fraction is lower than the coverage threshold
    #                        we stop iterating as at this point we have
    #                        accounted for all unclassified sequences
    #                        (everything that follows should have been
    #                        excluded)
    # - line 3 on (if present): fraction lower than coverage threshold -
    #                           as soon as we hit this we stop for the same
    #                           reason as above
    # If the report was fine, at this point we should arrive at 100%. If not,
    # something must have been wrong with the report - we raise an error.
    for mag in unclassified:
        report_fp = os.path.join(reports.path, f"{mag}.report.txt")
        unclassified_fraction = 0.0
        with open(report_fp, "r") as f:
            for line in f:
                fraction = float(line.split("\t")[0])
                if "unclassified" in line:
                    unclassified_fraction += fraction
                elif "root" in line:
                    unclassified_fraction += fraction
                    if fraction < coverage_threshold:
                        break
                elif fraction < coverage_threshold:
                    unclassified_fraction += fraction
                    break

        if unclassified_fraction != 100.0:
            raise ValueError(
                f"Unclassified fraction for MAG '{mag}' is not "
                f"100.0%: {unclassified_fraction:.2f}%. "
                "Please check the Kraken 2 report."
            )

        taxonomy.loc[mag, "Taxon"] = "d__Unclassified"

    return taxonomy


def kraken2_to_contig_taxonomy(
    ctx,
    reports,
    outputs,
    table,
    coverage_threshold = 0.1,
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
            - Taxonomy artifact: Series mapping contig IDs (index) to taxonomy
              strings (values). Unclassified contigs are assigned "d__Unclassified".
            - Taxonomy abundance table: Feature table with contig IDs replaced by
              taxonomy strings. Contigs with the same taxonomy are summed.
    """
    # Get taxonomy mapping from reports (taxon ID -> taxonomy string)
    _to_features = ctx.get_action("annotate", "kraken2_to_features")
    _, taxonomy_df = _to_features(reports, coverage_threshold)

    taxonomy_df = taxonomy_df.view(pd.DataFrame)
    outputs = outputs.view(Kraken2OutputDirectoryFormat)

    # Create a dictionary mapping taxon IDs to taxonomy strings
    # taxonomy_df has "Feature ID" (taxon ID) as index and "Taxon" as column
    taxon_to_taxonomy = taxonomy_df["Taxon"].to_dict()

    # Dictionary to store contig ID -> taxonomy mapping
    contig_taxonomy = {}

    # Iterate over all output files
    for relpath, output_df in outputs.outputs.iter_views(pd.DataFrame):
        # Skip empty files
        if output_df.empty:
            continue

        # Map each contig to its taxonomy
        contig_ids = output_df.iloc[:, 1].astype(str)
        taxon_ids = output_df.iloc[:, 2].astype(str)

        # Map with fallback for unclassified/missing taxons
        taxonomy = taxon_ids.where(taxon_ids != "0").map(taxon_to_taxonomy).fillna("d__Unclassified")

        # Update dictionary (last occurrence wins)
        contig_taxonomy.update(dict(zip(contig_ids, taxonomy)))

    # Convert to pandas Series
    taxonomy_series = pd.Series(contig_taxonomy, name="Taxon")
    taxonomy_series.index.name = "Feature ID"

    # Convert contig abundance table to taxonomy abundance table using biom.Table
    def _collapse(id_, md):
        try:
            tax = taxonomy_series.loc[id_]
        except KeyError:
            tax = "d__Unclassified"
        tax = [x.strip() for x in tax.split(';')]
        if len(tax) < 9:
            padding = ['__'] * (9 - len(tax))
            tax.extend(padding)
        return ';'.join(tax)

    contig_ft = table.view(biom.Table)
    contig_ft_collapsed = contig_ft.collapse(
        f=_collapse, axis="observation", norm=False
    )

    counts = {}
    for obs_id in contig_ft_collapsed.ids(axis="observation"):
        metadata = contig_ft_collapsed.metadata(id=obs_id, axis="observation")
        counts[obs_id] = len(metadata['collapsed_ids'])

    def divide_by_count(data, id_, metadata):
        # Find the index of the current observation to get the correct count
        return data / counts[id_]

    contig_ft_collapsed = contig_ft_collapsed.transform(divide_by_count, axis='observation')

    # x = table.view(pd.DataFrame)
    # y = x.rename(columns=contig_taxonomy, inplace=False)
    # z = y.T.groupby(level=0).mean()

    # Create taxonomy artifact
    taxonomy_artifact = ctx.make_artifact("FeatureData[Taxonomy]", taxonomy_series)
    ft_artifact = ctx.make_artifact("FeatureTable[Frequency]", contig_ft_collapsed)

    return taxonomy_artifact, ft_artifact


def kraken2_to_mag_features(
    reports: Kraken2ReportDirectoryFormat,
    outputs: Kraken2OutputDirectoryFormat,
    coverage_threshold: float = 0.1,
    # lca_mode: str = 'lca'
) -> pd.DataFrame:
    table, taxonomy = kraken2_to_features(reports, coverage_threshold)

    taxa_list = []
    # convert IDs to match MAGs instead of taxids/db ids
    for mag_id in table.index:
        kraken_table_fp = outputs.path / f"{mag_id}.output.txt"
        hits_df = pd.read_csv(kraken_table_fp, sep="\t", header=None, dtype="str")
        MAG_COL = 1
        TAXA_COL = 2

        mag_series = table.loc[mag_id, :]
        mag_obs = mag_series[mag_series != 0]
        merged_df = hits_df.join(mag_obs, on=TAXA_COL, how="right")
        merged_df = merged_df.join(taxonomy, on=TAXA_COL, how="left")

        new_taxa = merged_df[[MAG_COL, "Taxon"]].set_index(MAG_COL)
        new_taxa.index.name = "Feature ID"
        new_taxa["mag_id"] = mag_id
        taxa_list.append(new_taxa)

    taxonomy = _find_lcas(taxa_list, mode="lca")
    return _add_unclassified_mags(table, taxonomy, reports, coverage_threshold)


def kraken2_to_features(
    reports: Kraken2ReportDirectoryFormat, coverage_threshold: float = 0.1
) -> (pd.DataFrame, pd.DataFrame):

    rows = []
    trees = []
    unclassified = {}
    for relpath, df in reports.reports.iter_views(pd.DataFrame):
        sample_id = os.path.basename(relpath).replace(".report.txt", "")

        filtered = df[df["perc_frags_covered"] >= coverage_threshold]
        unclassified[sample_id] = filtered[filtered["name"] == "unclassified"]
        tree = _kraken_to_ncbi_tree(filtered)
        tips = _ncbi_tree_to_tips(tree)
        if tips:
            table_row = pd.Series(True, index=_ncbi_tree_to_tips(tree))
            table_row.name = sample_id
            rows.append(table_row)
        trees.append(tree)

    full_tree = _combine_ncbi_trees(trees)

    table = pd.DataFrame(rows).fillna(False)
    taxonomy = _to_taxonomy(full_tree)
    # filter taxonomy to only IDs in table
    # use list to avoid index name change
    taxonomy = taxonomy.loc[list(table.columns)]

    if table.empty:
        raise ValueError(
            "The resulting feature table was empty. Please adjust the "
            "coverage threshold and try again."
        )

    return table, taxonomy


def _get_indentation(string, indent=2):
    return (len(string) - len(string.lstrip(" "))) // indent


def _kraken_to_ncbi_tree(df):
    tree = skbio.TreeNode()
    stack = deque([(0, tree)])
    for _, row in df.iterrows():
        r = row["rank"]
        label = row["name"]
        otu = str(row["taxon_id"])

        if r in ("U", "R"):
            continue  # unclassified or root

        indent = _get_indentation(label)
        name = f"{r.lower()}__{label.strip()}"
        node = skbio.TreeNode(name=name, length=0.0)

        # Don't include internal non-strain infra-clades as tips
        if len(r) == 1 or r.startswith("S"):
            id_node = skbio.TreeNode(name=otu, length=0.0)
            node.length = 1.0  # not infra-clade, so give it a length
            node.append(id_node)

        parent_indent, parent_node = stack[-1]
        if parent_indent >= indent and parent_node.children:
            parent_node.children[0].is_actual_tip = True

        while parent_indent >= indent:
            stack.pop()
            parent_indent, parent_node = stack[-1]

        parent_node.append(node)
        stack.append((indent, node))

    # last entry is always a tip
    _, parent_node = stack[-1]
    # It is possible for the last row to be an infra-clade tip,
    # so walk backwards up the stack until a standard node (with length)
    # is found
    while stack and parent_node.length == 0:
        _, parent_node = stack.pop()

    if parent_node.children:
        parent_node.children[0].is_actual_tip = True

    return tree


def _combine_ncbi_trees(trees):
    full_tree = trees[0]
    for tree in trees[1:]:
        tip_cache = {t.name: t for t in full_tree.tips()}
        for tip in list(tree.tips()):
            # check if taxid is already in this tree
            if tip.name in tip_cache:
                continue  # for clarity
            else:
                parents = list(tip.ancestors())[:-1]  # ignore unnamed root

                # check if node is a infra-clade (infra-clades have length 0).
                # then adds self to ancestor list, if it is an infra-clade.
                # this mimics what happens if the node isn't an infra-clade.
                # i.e node had an id_node and then self gets added to the
                # list of ancestors if you call .parent on an id_node.
                if tip.length == 0:
                    parents.insert(0, tip)
                matching = full_tree
                subtree_inserted = False
                while parents and not subtree_inserted:
                    node = parents.pop()
                    ancestor_found = False

                    for child in matching.children:
                        if child.name == node.name:
                            matching = child
                            ancestor_found = True
                            break
                    if not ancestor_found:
                        matching.append(node)
                        # This may be overkill but this checks to make sure
                        # that the tip is an infra clade (tip.length = 0)
                        # and doesn't have children. If thats these are both
                        # true then add this tip to tip cache.
                        if len(tip.children) == 0 and tip.length == 0:
                            tip_cache[tip.name] = tip
                        for t in node.tips():
                            tip_cache[t.name] = t
                        assert tip.name in tip_cache
                        subtree_inserted = True

                # did we exhaust our parents and still not find ourselves?
                if not subtree_inserted and not (tip.name == matching.name):
                    # should be impossible. Implies ancestor_found was always
                    # True but not in tip_cache AND the tip is not an
                    # intermediate node in the tree that already is present.
                    raise AssertionError(f"{tip.name} could not be inserted")
    return full_tree


def _ncbi_tree_to_tips(tree):
    return [n.name for n in tree.tips() if hasattr(n, "is_actual_tip")]


def _pad_ranks(ranks):
    order = ["d", "k", "p", "c", "o", "f", "g", "s"]
    available = {}
    taxonomy = []

    for rank in reversed(ranks):
        r, label = rank.split("__", 1)
        available[r] = label
        if len(r) > 1 and r.startswith("s"):
            taxonomy.append(rank)

    for r in reversed(order):
        label = available.get(r)

        if label is not None:
            taxonomy.append(f"{r}__{label}")
            last_good_label = f"{r}__{label}"
        elif taxonomy:
            if r == "k" and available.get("d") in ("Bacteria", "Archaea"):
                # it is too strange to propagate upwards for these 'kingdoms'
                taxonomy.append(f"k__{available['d']}")
            else:
                # smear the most specific label we have upwards
                taxonomy.append(f"{r}__containing {last_good_label}")

    return ";".join(reversed(taxonomy))


def _to_taxonomy(tree):
    rows = [(node.name, _pad_ranks(ranks)) for node, ranks in tree.to_taxonomy()]
    taxonomy = pd.DataFrame(rows, columns=["Feature ID", "Taxon"])
    taxonomy = taxonomy.set_index("Feature ID")

    return taxonomy
