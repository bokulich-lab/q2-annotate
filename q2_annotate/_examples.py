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
