# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os.path
import shutil
from typing import List, Dict

import pandas as pd
import skbio
from q2_types.feature_data import DNAFASTAFormat
from scipy.cluster.hierarchy import ward, fcluster

from q2_types_genomics.feature_data import MAGSequencesDirFmt
from q2_types_genomics.per_sample_data import MultiMAGSequencesDirFmt


# def find_similar_bins(
#     distance_matrix: pd.DataFrame, threshold: float
# ) -> List[List[str]]:
#     """
#     Finds similar bins based on a distance matrix and
#     a similarity threshold.
#
#     Args:
#         distance_matrix (pd.DataFrame): A pandas DataFrame containing the
#                                           pairwise distances between bins.
#         threshold (float): A float representing the maximum distance for
#                               two bins to be considered similar.
#
#     Returns:
#         A list where each element is a list  of similar bins.
#
#     Notes:
#         The output will look like:
#             [['sample1_bin1', 'sample3_bin3'],
#              ['sample2_bin4, 'sample3_bin1']]
#         meaning that sample1_bin1 is similar to sample3_bin3 and
#         sample2_bin4 is similar to sample3_bin1, according to the
#         provided threshold.
#     """
#     sample_names = list(distance_matrix.index)
#     similar_bins = {}
#     for i, bin1 in enumerate(sample_names):
#         for j, bin2 in enumerate(sample_names[i + 1:]):
#             distance = distance_matrix.loc[bin1, bin2]
#             if distance <= threshold:
#                 if bin1 not in similar_bins:
#                     similar_bins[bin1] = []
#                 similar_bins[bin1].append(bin2)
#     similar_bins = [[key, *value] for key, value in similar_bins.items()]
#     return similar_bins


def _find_similar_bins_fcluster(
        distance_matrix: pd.DataFrame, threshold: float
) -> List[List[str]]:
    """
    Group bins into clusters based on a distance threshold.

    Args:
        distance_matrix (pd.DataFrame): A distance matrix.
        threshold (float): The distance threshold for forming clusters.

    Returns:
         A list where each element is a list of similar bins.

    Notes:
         The output will look like:
             [['sample1_bin1', 'sample3_bin3'],
              ['sample2_bin4, 'sample3_bin1']]
         meaning that sample1_bin1 is similar to sample3_bin3 and
         sample2_bin4 is similar to sample3_bin1, according to
         the provided threshold.
    """

    # Perform hierarchical/agglomerative clustering
    tree = ward(distance_matrix.values)

    # Form flat clusters from the hierarchical clustering defined
    # by the given linkage matrix
    cluster_ids = fcluster(tree, t=threshold, criterion='distance')

    # Map each MAG to its corresponding cluster
    clusters = {i: [] for i in cluster_ids}
    for i, cluster in enumerate(cluster_ids):
        clusters[cluster].append(distance_matrix.index[i])

    return list(clusters.values())


def _get_bin_lengths(mags: MultiMAGSequencesDirFmt) -> pd.Series:
    """
    Calculates the length of each bin in a MultiMAGSequencesDirFmt object.

    Args:
        mags (MultiMAGSequencesDirFmt): An object containing all the
                                        original bins from all samples.

    Returns:
        A pandas Series where the index is the bin name and the value
        is the length of the bin.
    """
    bin_lengths = {}
    for path, seq in mags.sequences.iter_views(DNAFASTAFormat):
        tot = 0
        for _seq in skbio.io.read(str(seq.path), format="fasta"):
            tot += len(_seq)
        bin_lengths[str(path)] = tot
    ser = pd.Series(bin_lengths, name="length")
    ser.index = ser.index.map(
        lambda x: x.replace(".fasta", "").split("/")[-1]
    )
    return ser


def _remap_bins(
    bin_clusters: List[List[str]],
    longest_bins: List[str],
    distances: pd.DataFrame
) -> Dict[str, str]:
    """
    Maps duplicate bins to a single dereplicated bin and assigns
    new IDs to the unique bins.

    Args:
        bin_clusters (list): A list of lists, where each inner list contains
                                the IDs of similar bins.
        longest_bins (str): A list of longest bin for each cluster.
        distances (pd.DataFrame): The original bin distance matrix.

    Returns:
        A dictionary where the keys are the original bin names and the
        values are tuples of the old and new bin IDs.

    Notes:
        The output will look like:
            {'sample1_bin1': 'mag5',
             'sample1_bin2': 'mag2',
             'sample2_bin1': 'mag8',
             'sample2_bin2': 'mag1',
             'sample2_bin3': 'mag6',
             'sample2_bin4': 'mag3',
             'sample2_bin5': 'mag7',
             'sample3_bin1': 'mag3',
             'sample3_bin2': 'mag4',
             'sample3_bin3': 'mag5'}
        meaning that mag5 was assigned to sample1_bin1 and sample3_bin3.
    """
    final_bins = {}
    for i, similar_bins in enumerate(bin_clusters):
        for bin in similar_bins:
            final_bins[bin] = longest_bins[i]
    for bin in distances.index:
        if bin not in final_bins:
            final_bins[bin] = bin

    return final_bins


def _reassign_bins_to_samples(
    final_bins: Dict[str, str], manifest: pd.DataFrame
) -> Dict[str, Dict[str, int]]:
    """
    Assigns bins to samples based on the final bin mapping.

    Args:
        final_bins (dict): A dictionary where the keys are the original
                            bin names and the values are the
                            new bin IDs.
        manifest (pd.DataFrame): Manifest of the original sample set -
                                    required to recover information
                                    about original sample IDs.

    Returns:
        A dictionary where the keys are sample IDs and the values are
        dictionaries where the keys are MAG IDs and the values are
        the number of bins assigned to that MAG.

    Notes:
        The output will look like:
            {'sample1': {'mag1': 0, 'mag2': 1, 'mag3': 0, 'mag4': 0,
                         'mag5': 1, 'mag6': 0, 'mag7': 0, 'mag8': 0},),
             'sample2': {'mag1': 1, 'mag2': 0, 'mag3': 1, 'mag4': 0,
                         'mag5': 0, 'mag6': 1, 'mag7': 1, 'mag8': 1},),
             'sample3': {'mag1': 0, 'mag2': 0, 'mag3': 1, 'mag4': 1,
                         'mag5': 1, 'mag6': 0, 'mag7': 0, 'mag8': 0},)}.
    """
    all_samples = manifest.copy(deep=True) \
        .reset_index().replace({"mag-id": final_bins})
    all_derep_mags = set([mag_id for _, mag_id in final_bins.items()])

    samples_to_bins = {
        key: {mag: 0 for mag in all_derep_mags}
        for key in set(all_samples["sample-id"])
    }

    for i, row in all_samples.iterrows():
        sample_id = row["sample-id"]
        mag_id = row["mag-id"]
        samples_to_bins[sample_id][mag_id] += 1

    return samples_to_bins


def _write_unique_bins(
    all_bins: MultiMAGSequencesDirFmt, bins_remapped: Dict[str, str]
) -> MAGSequencesDirFmt:
    """
    Writes the unique bins to a new MAGSequencesDirFmt object
    based on the final bin mapping.

    Args:
        all_bins (MultiMAGSequencesDirFmt): An object with all the bins.
        bins_remapped (dict): A dictionary where the keys are the original
                                bin names and the values are tuples of
                                the old and new bin IDs.

    Returns:
        A MAGSequencesDirFmt object containing the unique bins.
    """
    derep_bins = MAGSequencesDirFmt()
    manifest = all_bins.manifest.view(pd.DataFrame)
    manifest.index = manifest.index.droplevel(0)
    for old_bin_id, new_bin_id in bins_remapped.items():
        dst_bin = os.path.join(str(derep_bins), f"{new_bin_id}.fasta")
        if not os.path.isfile(dst_bin):
            src_bin = manifest.loc[new_bin_id, "filename"]
            shutil.copy(src_bin, dst_bin)
    return derep_bins


def _generate_pa_table(
    unique_bins_per_sample: Dict[str, Dict[str, int]]
) -> pd.DataFrame:
    """
    Generates a presence-absence table from a dictionary of unique
    bins per sample.

    Args:
        unique_bins_per_sample: A dictionary where the keys are sample IDs
                                and the values are lists of unique bin IDs.

    Returns:
        A pandas DataFrame where the index is the sample ID and the columns
        are the unique bin IDs, with 1 indicating presence and 0 indicating
        absence.
    """
    presence_absence = pd.DataFrame.from_records(unique_bins_per_sample).T
    presence_absence = presence_absence.astype(bool).astype(int)
    presence_absence.index.name = "sample-id"
    presence_absence.sort_index(inplace=True)
    presence_absence = presence_absence[sorted(presence_absence.columns)]
    return presence_absence


def dereplicate_mags(
    mags: MultiMAGSequencesDirFmt,
    distance_matrix: skbio.DistanceMatrix,
    threshold: float = 0.99,
) -> (MAGSequencesDirFmt, pd.DataFrame):
    distances = distance_matrix.to_data_frame()

    # find similar bins, according to the threshold
    bin_clusters = _find_similar_bins_fcluster(distances, threshold)

    # find the longest bin in each cluster
    bin_lengths = _get_bin_lengths(mags)
    longest_bins = [bin_lengths[ids].idxmax() for ids in bin_clusters]

    # generate a map between the original bins and the dereplicated bins
    final_bins = _remap_bins(bin_clusters, longest_bins, distances)

    # generate dereplicated bin sequences
    unique_bin_seqs = _write_unique_bins(mags, final_bins)

    # generate a presence-absence table
    unique_bins_per_sample = _reassign_bins_to_samples(
        final_bins, mags.manifest.view(pd.DataFrame)
    )
    presence_absence = _generate_pa_table(unique_bins_per_sample)

    return unique_bin_seqs, presence_absence
