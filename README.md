# q2-moshpit
![CI](https://github.com/bokulich-lab/q2-moshpit/actions/workflows/ci.yaml/badge.svg)
[![codecov](https://codecov.io/gh/bokulich-lab/q2-moshpit/graph/badge.svg?token=PSCAYJUP01)](https://codecov.io/gh/bokulich-lab/q2-moshpit)

**MO**dular **SH**otgun metagenome **P**ipelines with **I**ntegrated provenance **T**racking
QIIME 2 plugin for functional annotation and taxonomic classification of shotgun metagenomes.

## Installation
_q2-moshpit_ is available as part of the QIIME 2 metagenome distribution. For installation and usage instructions please consult the official [QIIME 2 documentation](https://www.docs.qiime2.org). 

## Functionality
This QIIME 2 plugin contains actions used to annotate and classify (meta)genomes.
Below you will find an overview of actions available in the plugin.

| Action               | Description                                                | Underlying tool                                                    |
|----------------------|------------------------------------------------------------|--------------------------------------------------------------------|
| bin-contigs-metabat  | Bin contigs into MAGs using MetaBat 2.                     | [MetaBat 2](https://bitbucket.org/berkeleylab/metabat/src/master/) |
| build-custom-diamond-db | Create a DIAMOND reference database from a FASTA input file. | [Diamond](https://github.com/bbuchfink/diamond) |
| build-eggnog-diamond-db | Create a DIAMOND reference database for the specified taxon. | [Diamond](https://github.com/bbuchfink/diamond) |
| build-kraken-db      | Fetch an existing or build a custom Kraken 2 database.     | [Kraken 2](https://ccb.jhu.edu/software/kraken2/)                  |
| classify-kaiju | Classify reads using Kaiju. | [Kaiju](https://bioinformatics-centre.github.io/kaiju/) |
| classify-kraken2      | Classify reads/MAGs using Kraken 2.                        | [Kraken 2](https://ccb.jhu.edu/software/kraken2/)                  |
| dereplicate-mags | Dereplicate MAGs from multiple samples. | - |
| eggnog-annotate          | Annotate orthologs against eggNOG database. | [EggNOG mapper](https://github.com/eggnogdb/eggnog-mapper) |
| eggnog-diamond-search    | Run eggNOG search using diamond aligner. | [EggNOG mapper](https://github.com/eggnogdb/eggnog-mapper) |
| eggnog-hmmer-search      | Run eggNOG search using HMMER aligner. | [EggNOG mapper](https://github.com/eggnogdb/eggnog-mapper) |
| estimate-bracken         | Perform read abundance re-estimation using Bracken. | [Kraken 2](https://ccb.jhu.edu/software/bracken/) |
| estimate-mag-abundance   | Estimate MAG abundance. | - | 
| evaluate-busco           | Evaluate quality of the generated MAGs using BUSCO. | [BUSCO](https://busco.ezlab.org) |
| extract-annotations      | Extract annotation frequencies from all annotations. | - |
| fetch-busco-db           | Download BUSCO database. | [BUSCO](https://busco.ezlab.org) |
| fetch-diamond-db         | Fetch the complete Diamond database necessary to run the eggnog-diamond-search action. | [EggNOG mapper](https://github.com/eggnogdb/eggnog-mapper) |
| fetch-eggnog-db          | Fetch the databases necessary to run the eggnog-annotate action. | [EggNOG mapper](https://github.com/eggnogdb/eggnog-mapper) |
| fetch-eggnog-hmmer-db    | Fetch the taxon specific database necessary to run the eggnog-hmmer-search action. | [EggNOG mapper](https://github.com/eggnogdb/eggnog-mapper) |
| fetch-eggnog-proteins    | Fetch the databases necessary to run the build-eggnog-diamond-db action. | [EggNOG mapper](https://github.com/eggnogdb/eggnog-mapper) |
| fetch-kaiju-db           | Fetch Kaiju database. | [Kaiju](https://bioinformatics-centre.github.io/kaiju/) |
| fetch-ncbi-taxonomy      | Fetch NCBI reference taxonomy. | [EggNOG mapper](https://github.com/eggnogdb/eggnog-mapper) |
| filter-derep-mags        | Filter dereplicated MAGs. | - |
| filter-mags              | Filter MAGs. | - |
| filter-reads-pangenome   | Remove contaminating human reads. | [Bowtie 2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) |
| get-feature-lengths      | Get feature lengths. | - |
| inspect-kraken2-db       | Inspect a Kraken 2 database. |
| kraken2-to-features      | Select downstream features from Kraken 2. | - |
| kraken2-to-mag-features  | Select downstream MAG features from Kraken 2. | - |
| multiply-tables          | Multiply two feature tables. | - |
| predict-genes-prodigal   | Predict gene sequences from MAGs using Prodigal. | [Prodigal](https://github.com/hyattpd/Prodigal) |
