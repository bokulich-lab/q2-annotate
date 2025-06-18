import os
import glob
from qiime2 import Artifact
from q2_types.feature_data import DNAFASTAFormat
from q2_types.feature_data import DNAIterator

from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt, ContigSequencesDirFmt
# from skbio import read
# from skbio.sequence import DNA

from q2_types.feature_data import DNAIterator
from q2_types.per_sample_sequences import MultiFASTADirectoryFormat
import os
import shutil

import numpy as np
import pandas as pd
from qiime2.util import duplicate

from q2_types._util import _validate_num_partitions
# from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt
### new for filter_contigs
from q2_assembly.filter import filter_contigs
from qiime2 import Metadata
#filter_contigs, input dataframe->metadata that fits as the input of filter contigs (id:index-ids)
def partition_sample_data_mags(
    mags_qza, num_partitions: int = None
) -> MultiMAGSequencesDirFmt:
    """
    Returns a dictionary where each key is either the mag_id or an index, and
    values are the new objects with the mags.
    """
    mags = Artifact.load(mags_qza).view(MultiMAGSequencesDirFmt)
    partitioned_mags = {}
    mags_all = [{k: v} for k, v in mags.sample_dict().items()]

    num_partitions = _validate_num_partitions(
        len(mags_all), num_partitions, "sample"
    )

    arrays_of_mags = np.array_split(mags_all, num_partitions)

    for i, samples in enumerate(arrays_of_mags, 1):
        result = MultiMAGSequencesDirFmt()
        all_samples = set(k for d in samples for k in d.keys())
        manifest = pd.read_csv(mags.path / "MANIFEST", index_col=None)
        manifest = manifest[manifest["sample-id"].isin(all_samples)]
        manifest.to_csv(result.path / "MANIFEST", index=False)

        for sample_dict in samples:
            for sample_id, mag_dict in sample_dict.items():
                for mag_id, mag_fp in mag_dict.items():
                    os.makedirs(result.path / sample_id, exist_ok=True)
                    duplicate(
                        mag_fp,
                        result.path / sample_id / os.path.basename(mag_fp)
                    )

        # If num_partitions == num_samples we will only have gone through one
        # sample in the above loop and will use its id as a key. Otherwise we
        # may have gone through multiple MAGs in the above loop and will be
        # using indices for keys
        if num_partitions == len(mags_all):
            partitioned_mags[sample_id] = result
        else:
            partitioned_mags[i] = result
    
    print("\nFinal Partitioned MAGs Structure:")
    print(partitioned_mags)
    return partitioned_mags
def partition_sample_data_unbinned(
    unbinned_qza, mags_partitioned: MultiMAGSequencesDirFmt
) ->ContigSequencesDirFmt:
    """
    Partition unbinned contigs to match the partitioning of MAGs.
    
    Args:
        unbinned (ContigSequencesDirFmt): The unbinned contigs directory.
        mags_partitioned (MultiMAGSequencesDirFmt): The partitioned MAGs structure.
        
    Returns:
        Dict[str, ContigSequencesDirFmt]: Mapping of partitions to unbinned contigs.
    """
    unbinned = Artifact.load(unbinned_qza).view(ContigSequencesDirFmt)
    partitioned_unbinned = {}

    for partition_id, mag_partition in mags_partitioned.items():  # Loop over MAG partitions
        result = ContigSequencesDirFmt()  # Temporary storage for partitioned unbinned contigs
        for sample_id in mag_partition.sample_dict().keys():  # Extract sample names from MAGs
            expected_unbinned_filename = f"{sample_id}_contigs.fa"
            unbinned_sample_dict = unbinned.sample_dict()  # Get unbinned sample mapping
            
            if sample_id in unbinned_sample_dict:  # Check if unbinned contigs exist for this sample
                unbinned_fp = unbinned_sample_dict[sample_id]
                os.makedirs(result.path, exist_ok=True)
                duplicate(unbinned_fp, result.path / expected_unbinned_filename)
        
        partitioned_unbinned[partition_id] = result  # Store the partitioned unbinned contigs
    print("\nFinal Partitioned unbinned Structure:")
    print(partitioned_unbinned)
    return partitioned_unbinned

def partition_filtered_unbinned(
    unbinned_qza, mags_partitioned: MultiMAGSequencesDirFmt
) ->ContigSequencesDirFmt:
    """
    Partition unbinned contigs to match the partitioning of MAGs.
    
    Args:
        unbinned (ContigSequencesDirFmt): The unbinned contigs directory.
        mags_partitioned (MultiMAGSequencesDirFmt): The partitioned MAGs structure.
        
    Returns:
        Dict[str, ContigSequencesDirFmt]: Mapping of partitions to unbinned contigs.
    """
    unbinned = Artifact.load(unbinned_qza).view(ContigSequencesDirFmt)
    partitioned_unbinned = {}

    for partition_id, mag_partition in mags_partitioned.items():
        print("\npartition_id:", partition_id)
        print("\mag_partition:", mag_partition)
        print("\mag_partition.keys():", mags_partitioned.keys())
        print("\mag_partition.sample_dict().keys():", mag_partition.sample_dict().keys())
        sample_ids = list(mag_partition.sample_dict().keys())
        index = pd.Index(sample_ids, name="ID")
        metadata = Metadata(pd.DataFrame(index=index))
        # Build the WHERE string manually
        id_list = ", ".join([f"'{sid}'" for sid in sample_ids])
        where = f"ID IN ({id_list})"
        # Filter unbinned using metadata
        filtered_unbinned = filter_contigs(
            contigs=unbinned,
            metadata=metadata,
            where=where
        )
        partitioned_unbinned[partition_id] = filtered_unbinned
    # for unbinned_id, unbinned_value in filtered_unbinned.items():
    #     # print("\unbinned_id:", unbinned_id)
    #     print("\unbinned_value:", unbinned_value)
    #     print("\unbinned_value.sample_dict().keys():", unbinned_value.sample_dict().keys())
        # Save this filtered subset
    # partitioned_unbinned[partition_id] = filtered_unbinned
    print("filtered_unbinned:", filtered_unbinned)
    print("\nFinal Partitioned unbinned Structure:")
    print(partitioned_unbinned)
    return partitioned_unbinned
# def count_binned_contigs(bins: MultiFASTADirectoryFormat) -> int:
#     """Counts sequences in all FASTA files inside a MultiFASTADirectoryFormat using DNAIterator."""
    
#     total_sequences = 0

#     #Iterate over sequences directly using QIIME 2's DNAIterator
#     for fasta_fp, dna_iterator in bins.sequences.iter_views(DNAIterator):
#         sequence_count = sum(1 for _ in dna_iterator)  # Count sequences
#         # print(f"{fasta_fp} contains {sequence_count} sequences")  # Debugging output
#         total_sequences += sequence_count

#     return total_sequences

# def count_unbinned_contigs(unbinned: ContigSequencesDirFmt) -> int:
#     """Counts sequences in the unbinned FASTA file using DNAIterator."""
    
#     total_sequences = 0

#     for fasta_fp, dna_iterator in unbinned.sequences.iter_views(DNAIterator):
#         sequence_count = sum(1 for _ in dna_iterator)  # Count sequences
#         # print(f"{fasta_fp} contains {sequence_count} sequences")  # Debugging output
#         total_sequences += sequence_count

#     return total_sequences

# def calculate_unbinned_percentage(mags_qza, unbinned_qza):
#     """Calculates the percentage of unbinned contigs using skbio instead of manual counting."""
    
#     # Load QIIME2 artifacts and get correct directory formats
#     mags_dir = Artifact.load(mags_qza).view(MultiMAGSequencesDirFmt)  # Binned MAGs
#     unbinned_dir = Artifact.load(unbinned_qza).view(ContigSequencesDirFmt)  # Unbinned Contigs

#     # Count sequences using skbio
#     binned_contigs = count_binned_contigs(mags_dir)
#     unbinned_contigs_count = count_unbinned_contigs(unbinned_dir)

#     # Calculate percentage
#     total = binned_contigs + unbinned_contigs_count
#     percentage_unbinned = (unbinned_contigs_count / total) * 100 if total > 0 else 0

#     # Print results for verification
#     print(f"Total Contigs (Binned + Unbinned): {total}")
#     print(f"Binned Contigs: {binned_contigs}")
#     print(f"Unbinned Contigs: {unbinned_contigs_count}")
#     print(f"Percentage of Unbinned Contigs: {percentage_unbinned:.2f}%")

#     return percentage_unbinned 

# # Run test
# calculate_unbinned_percentage("./data/mags.qza", "./data/unbinned_contigs.qza")
mags_partitioned = partition_sample_data_mags("./data/mags.qza", 2)

# for partition_id, mag_partition in mags_partitioned.items():
#     output_path = f"./partitioned_output/partition_{partition_id}"
#     mag_partition.save(output_path)

# partition_filtered_unbinned("./data/unbinned_contigs.qza", mags_partitioned)