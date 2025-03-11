import os
import glob
from qiime2 import Artifact
from q2_types.feature_data import DNAFASTAFormat
from q2_types.feature_data import DNAIterator

from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt, ContigSequencesDirFmt
from skbio import read
from skbio.sequence import DNA

from q2_types.feature_data import DNAIterator
from q2_types.per_sample_sequences import MultiFASTADirectoryFormat

def count_binned_contigs(bins: MultiFASTADirectoryFormat) -> int:
    """Counts sequences in all FASTA files inside a MultiFASTADirectoryFormat using DNAIterator."""
    
    total_sequences = 0

    #Iterate over sequences directly using QIIME 2's DNAIterator
    for fasta_fp, dna_iterator in bins.sequences.iter_views(DNAIterator):
        sequence_count = sum(1 for _ in dna_iterator)  # Count sequences
        # print(f"{fasta_fp} contains {sequence_count} sequences")  # Debugging output
        total_sequences += sequence_count

    return total_sequences

def count_unbinned_contigs(unbinned: ContigSequencesDirFmt) -> int:
    """Counts sequences in the unbinned FASTA file using DNAIterator."""
    
    total_sequences = 0

    for fasta_fp, dna_iterator in unbinned.sequences.iter_views(DNAIterator):
        sequence_count = sum(1 for _ in dna_iterator)  # Count sequences
        # print(f"{fasta_fp} contains {sequence_count} sequences")  # Debugging output
        total_sequences += sequence_count

    return total_sequences

def calculate_unbinned_percentage(mags_qza, unbinned_qza):
    """Calculates the percentage of unbinned contigs using skbio instead of manual counting."""
    
    # Load QIIME2 artifacts and get correct directory formats
    mags_dir = Artifact.load(mags_qza).view(MultiMAGSequencesDirFmt)  # Binned MAGs
    unbinned_dir = Artifact.load(unbinned_qza).view(ContigSequencesDirFmt)  # Unbinned Contigs

    # Count sequences using skbio
    binned_contigs = count_binned_contigs(mags_dir)
    unbinned_contigs_count = count_unbinned_contigs(unbinned_dir)

    # Calculate percentage
    total = binned_contigs + unbinned_contigs_count
    percentage_unbinned = (unbinned_contigs_count / total) * 100 if total > 0 else 0

    # Print results for verification
    print(f"Total Contigs (Binned + Unbinned): {total}")
    print(f"Binned Contigs: {binned_contigs}")
    print(f"Unbinned Contigs: {unbinned_contigs_count}")
    print(f"Percentage of Unbinned Contigs: {percentage_unbinned:.2f}%")

    return percentage_unbinned 

# Run test
calculate_unbinned_percentage("./data/mags.qza", "./data/unbinned_contigs.qza")


# # Run test
# calculate_unbinned_percentage("./data/mags.qza", "./data/unbinned_contigs.qza")
