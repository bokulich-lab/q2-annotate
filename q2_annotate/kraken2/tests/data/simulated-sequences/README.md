## kraken2-db/
A kraken2 database created using ~20k nt of refseq genomes for the species
Bacillus anthracis, Mus musculus, Staphylococcus aureus, and Staphylococcus
epidermidis, the ncbi taxonomy, and all default parameters (kraken2 version 
2.1.3).

## reads/
Simulated paired-end reads taken from the above reference genomes. The nonsense
reads are randomly generated. The ba-mm-mixed-* reads are an approximately
75% Bacillus anthracis, 25% Mus musculus read mixture meant to mimic reads with
host contamination.

## contigs/
Simulated contigs taken from the above reference genomes. Contigs are not
overlapping. The ba-mm-mixed_contigs are also 75/25 ba/mm.  

## mags/
Simulated MAGs constructed from the contigs by copying most contigs for each
sample. 