# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

url = "https://scop.berkeley.edu/downloads/scopeseq-2.07/astral-scopedom-seqres"
"-gd-sel-gs-bib-40-2.07.fa"


def diamond_makedb(use):
    fasta_input = use.init_artifact_from_url("sequences", url)

    _ = use.action(
        use.UsageAction("annotate", "build_custom_diamond_db"),
        use.UsageInputs(
            seqs=fasta_input,
        ),
        use.UsageOutputNames(
            db="diamond_db",
        ),
    )


def build_kaiju_db_example(use):
    # This example shows how to build a custom Kaiju database from protein sequences
    # Note: This is a simplified example - in practice you would have real genome protein data
    
    proteins = use.init_artifact(
        'proteins',
        # Placeholder for GenomeData[Proteins] - would contain protein FASTA files
        # for each genome, e.g., genome1.fasta, genome2.fasta, etc.
    )
    
    metadata = use.init_metadata(
        'genome_metadata',
        # Metadata mapping genome IDs to NCBI taxon IDs
        # Format: genome_id -> taxon_id (e.g., genome1 -> 511145)
    )
    
    kaiju_db_result = use.action(
        use.UsageAction("annotate", "build_kaiju_db"),
        use.UsageInputs(
            proteins=proteins,
            metadata=metadata,
        ),
        use.UsageOutputNames(
            db="custom_kaiju_db",
        ),
    )
