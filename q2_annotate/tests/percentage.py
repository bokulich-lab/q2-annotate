import os
import glob
from qiime2 import Artifact
from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt, ContigSequencesDirFmt

def count_fasta_sequences(directory):
    """Counts the number of sequences in all FASTA files inside a directory."""
    fasta_files = glob.glob(os.path.join(directory, "**", "*.fa*"), recursive=True)
    if not fasta_files:
        print(f"No FASTA files found in directory: {directory}")
    total_sequences = 0
    for file in fasta_files:
        # print(f"Found FASTA file: {file}")
        with open(file, "r") as f:
            sequence_count = sum(1 for line in f if line.startswith(">"))
            # print(f"{file} contains {sequence_count} contigs")
            total_sequences += sequence_count   # Count headers
    return total_sequences

def calculate_unbinned_percentage(mags_qza, unbinned_qza):
    #Load the artifacts and get correct directory formats
    mags_dir = Artifact.load(mags_qza).view(MultiMAGSequencesDirFmt)  # Get MAGs directory
    unbinned_dir = Artifact.load(unbinned_qza).view(ContigSequencesDirFmt)  # Get unbinned directory

    #Convert QIIME2 directory format to actual file system path
    mags_path = str(mags_dir)
    unbinned_path = str(unbinned_dir)
    # print(f"Binned MAGs directory: {mags_path}")
    # print(f"Unbinned Contigs directory: {unbinned_path}")
    #Count sequences
    total_contigs = count_fasta_sequences(mags_path)
    unbinned_contigs_count = count_fasta_sequences(unbinned_path)

    #Calculate percentage
    total = total_contigs + unbinned_contigs_count
    percentage_unbinned = (unbinned_contigs_count / total) * 100 if total > 0 else 0

    #Print results for verification
    print(f"Total Contigs (Binned + Unbinned): {total}")
    print(f"Binned Contigs: {total_contigs}")
    print(f"Unbinned Contigs: {unbinned_contigs_count}")
    print(f"Percentage of Unbinned Contigs: {percentage_unbinned:.2f}%")

    return percentage_unbinned 

#Run test
calculate_unbinned_percentage("./data/mags.qza", "./data/unbinned_contigs.qza")
