## kraken2-db/
A Kraken 2 database created using ~20 kbp nt of RefSeq genomes for the species
Bacillus anthracis, Mus musculus, Staphylococcus aureus, and Staphylococcus
epidermidis, the NCBI taxonomy, and all default parameters (Kraken 2 version 
2.1.3).

## reads/
Simulated paired-end reads taken from the above reference genomes. The nonsense
reads are randomly generated. The `ba-mm-mixed-*` reads are an approximately
75% Bacillus anthracis, 25% Mus musculus read mixture meant to mimic reads with
host contamination.

## contigs/
Simulated contigs taken from the above reference genomes. Contigs are not
overlapping. The `ba-mm-mixed_contigs` are also 75/25 ba/mm.  

## mags/
Simulated MAGs constructed from the contigs by copying most contigs for each
sample.

## coassembled-contigs/
Two coassembly artifacts, each containing a single `pooled.fasta`. Artifact 1
contains the first contig from `contigs/ba_contigs.fasta` and
`contigs/mm_contigs.fasta`; artifact 2 contains the first contig from
`contigs/sa_contigs.fasta` and `contigs/se_contigs.fasta`. Sequences are unchanged,
and IDs are prefixed with their source sample (`ba_contig0`, `mm_contig0`,
`sa_contig0`, `se_contig0`) to keep them unique across artifacts.
