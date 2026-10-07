# q2-annotate
![CI](https://github.com/bokulich-lab/q2-annotate/actions/workflows/ci.yaml/badge.svg)
[![codecov](https://codecov.io/gh/bokulich-lab/q2-annotate/graph/badge.svg?token=PSCAYJUP01)](https://codecov.io/gh/bokulich-lab/q2-annotate)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

QIIME 2 plugin for functional annotation and taxonomic classification of shotgun metagenomes.

## Taxonomic functional attribution

`estimate-tfa` returns two artifacts:

- `tfa`: a sparse `FeatureTable[Frequency % Properties("tfa")]` with the original sample
  IDs as columns and one feature per observed taxon/gene pair.
- `gene_taxonomy`: `FeatureData[Taxonomy % Properties("tfa")]`, a standard
  taxonomy TSV with `Feature ID` and `Taxon` followed by `Taxon ID` and `Gene ID`.

Feature IDs concatenate the taxon and gene IDs with `|`, for example
`T1|bla_TEM`. Each component is percent-encoded, so a literal `|` becomes
`%7C`, `%` becomes `%25`, and a space becomes `%20`. This keeps IDs readable
and distinct even when an input ID contains the separator. IDs are stable
across runs and input ordering. The mapping retains the original, unescaped
IDs and taxonomy labels.

The `tfa` property marks the frequency table and its companion taxonomy as
representing taxon/gene pairs. They remain compatible with actions accepting
ordinary frequency tables and taxonomy artifacts. The taxonomy's semantic
validator requires nonempty IDs and labels and unique taxon/gene pairs.
TFA-specific actions require the property on both inputs.

For taxon *t*, gene *g*, and sample *s*, a load is the sum over contigs assigned
to *t* of contig abundance in *s* multiplied by the contig's count of *g*.
Only pairs with a nonzero load in at least one sample become rows; sample zeros
stay sparse. Estimation reports an error when no pairs remain, because the
standard taxonomy format requires at least one record.
Fractional abundance-weighted values are preserved without rounding,
following the `FeatureTable[Frequency]` convention used by `q2-mag
estimate-abundance` for RPKM/TPM. Both outputs must be retained to interpret
the table's feature IDs.

`explore-tfa` accepts the frequency table and its required gene-taxonomy
mapping. Taxonomy labels from the mapping appear in the taxon selector and
heatmap. Choose a taxonomy level to collapse taxa sharing that lineage
prefix. Loads are summed within each sample before the heatmap statistic is
computed; the taxon dropdown and all plots use the collapsed groups. Missing
ranks appear as Unclassified under their known lineage. Original taxa keeps
the uncollapsed view.
The full lineage appears in the selector and heatmap tooltip, while the
heatmap axis shows the terminal taxon name. Its heatmap ranks and colors
taxon/function links by the selected sum, mean, or median load across all
samples, including zero loads; selecting a link shows its per-sample
loads as a boxplot and a sample bar chart. Select a taxon or function alone
to sum its matching pairs within each sample. Optional sample metadata adds
categorical columns to the grouping control, so the boxplot can compare
sample groups. Samples with zero load and samples missing
metadata remain in the comparison. The heatmap displays up to 100 links, and
the bar chart shows the 100 highest-load samples; the boxplot uses all samples.
The Vega plots use
the same CDN scripts as the existing Kraken 2 visualizer, so viewing them
requires an internet connection.

A larger, reproducible synthetic dataset and a runnable command are in
[`examples/tfa_explorer/README.md`](examples/tfa_explorer/README.md).

The sample-resolved output replaces the older `estimate-tfa` layout, which
summed abundances across samples and used functional feature IDs as sample IDs.
Run `estimate-tfa` again from the original inputs to obtain per-sample loads;
the original layout cannot recover them.

## Installation
_q2-annotate_ is available as part of the QIIME 2 moshpit distribution. For installation and usage instructions please consult the official [QIIME 2 documentation](https://docs.qiime2.org).

## Functionality
This QIIME 2 plugin contains actions used to annotate and classify (meta)genomes:

| Action                    | Description                                                                            | Underlying tool       |
|---------------------------|----------------------------------------------------------------------------------------|-----------------------|
| build-custom-diamond-db   | Create a DIAMOND reference database from a FASTA input file.                           | [Diamond](https://github.com/bbuchfink/diamond) |
| build-eggnog-diamond-db   | Create a DIAMOND reference database for the specified taxon.                           | [Diamond](https://github.com/bbuchfink/diamond) |
| build-kraken-db           | Fetch an existing or build a custom Kraken 2 database.                                 | [Kraken 2](https://ccb.jhu.edu/software/kraken2/)      |
| classify-kaiju            | Classify reads using Kaiju.                                                            | [Kaiju](https://bioinformatics-centre.github.io/kaiju/) |
| classify-kraken2          | Classify reads/MAGs using Kraken 2.                                                    | [Kraken 2](https://ccb.jhu.edu/software/kraken2/)      |
| construct-pangenome-index | Construct the human pangenome Bowtie 2 index.                                          | [Bowtie 2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) |
| estimate-bracken          | Perform read abundance re-estimation using Bracken.                                    | [Kraken 2](https://ccb.jhu.edu/software/bracken/) |
| estimate-tfa              | Estimate a sparse load table and gene taxonomy mapping. | - |
| explore-tfa               | Explore TFA links and compare sample groups interactively. | [Vega-Lite](https://vega.github.io/vega-lite/) |
| extract-annotations       | Extract annotation frequencies from all annotations.                                   | - |
| fetch-diamond-db          | Fetch the complete Diamond database necessary to run the eggnog-diamond-search action. | [EggNOG mapper](https://github.com/eggnogdb/eggnog-mapper) |
| fetch-eggnog-db           | Fetch the databases necessary to run the eggnog-annotate action.                       | [EggNOG mapper](https://github.com/eggnogdb/eggnog-mapper) |
| fetch-eggnog-hmmer-db     | Fetch the taxon specific database necessary to run the eggnog-hmmer-search action.     | [EggNOG mapper](https://github.com/eggnogdb/eggnog-mapper) |
| fetch-eggnog-proteins     | Fetch the databases necessary to run the build-eggnog-diamond-db action.               | [EggNOG mapper](https://github.com/eggnogdb/eggnog-mapper) |
| fetch-kaiju-db            | Fetch Kaiju database.                                                                  | [Kaiju](https://bioinformatics-centre.github.io/kaiju/) |
| fetch-ncbi-taxonomy       | Fetch NCBI reference taxonomy.                                                         | [EggNOG mapper](https://github.com/eggnogdb/eggnog-mapper) |
| filter-kraken2-reads      | Filter Kraken 2-classified reads by taxonomy (for example, retain only bacterial reads). | - |
| filter-kraken2-results    | Filter Kraken 2 reports and outputs by metadata and/or abundance                       | - |
| filter-reads-pangenome    | Remove contaminating human reads.                                                      | [Bowtie 2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) |
| inspect-kraken2-db        | Inspect a Kraken 2 database.                                                           | - |
| kraken2-to-features       | Select downstream features from Kraken 2.                                              | - |
| kraken2-to-mag-features   | Select downstream MAG features from Kraken 2.                                          | - |
| map-eggnog                | Annotate orthologs against eggNOG database.                                            | [EggNOG mapper](https://github.com/eggnogdb/eggnog-mapper) |
| multiply-tables           | Multiply two feature tables.                                                           | - |
| predict-genes-prodigal    | Predict gene sequences from MAGs using Prodigal.                                       | [Prodigal](https://github.com/hyattpd/Prodigal) |
| search-orthologs-diamond  | Run eggNOG search using Diamond aligner.                                               | [EggNOG mapper](https://github.com/eggnogdb/eggnog-mapper) |
| search-orthologs-hmmer    | Run eggNOG search using HMMER aligner.                                                 | [EggNOG mapper](https://github.com/eggnogdb/eggnog-mapper) |
