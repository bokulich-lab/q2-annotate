# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ---------------------------------------------------------------------------
import glob
import os.path
import re
import shutil
import tempfile
from uuid import uuid4

from q2_types.per_sample_sequences import (
    BAMDirFmt,
    ContigSequencesDirFmt,
    MultiFASTADirectoryFormat,
)

from q2_annotate._utils import _process_common_input_params, run_command
from q2_annotate.metabat2.metabat2 import _assert_samples, _generate_contig_map
from q2_annotate.semibin2.utils import _process_semibin2_arg


def _run_semibin2(samp_name, samp_props, loc, mode, common_args):
    bins_dp = os.path.join(loc, samp_name)
    bins_prefix = os.path.join(bins_dp, "bin")
    os.makedirs(bins_dp)
    cmd = [
        "SemiBin2",
        mode,
        "--input-fasta",
        samp_props["contigs"],
        "--input-bam",
        samp_props["map"],
        "--output",
        bins_prefix,
        "--compression",
        "none",
        "--verbose",
    ]
    cmd.extend(common_args)
    run_command(cmd, verbose=True)
    return bins_dp


def _process_sample(samp_name, samp_props, mode, common_args, result_loc):
    with tempfile.TemporaryDirectory() as tmp:
        output_dp = _run_semibin2(samp_name, samp_props, tmp, mode, common_args)
        bins_dp = os.path.join(output_dp, "bin/output_bins/")

        all_outputs = glob.glob(os.path.join(bins_dp, "*.fa"))
        all_bins = [x for x in all_outputs if re.search(r"SemiBin\_[0-9]+\.fa$", x)]

        # rename using UUID v4
        bin_dest_dir = os.path.join(str(result_loc), samp_name)
        os.makedirs(bin_dest_dir, exist_ok=True)
        for old_bin in all_bins:
            new_bin = os.path.join(bin_dest_dir, f"{uuid4()}.fa")
            shutil.move(old_bin, new_bin)


def _bin_contigs_semibin2(
    contigs: ContigSequencesDirFmt,
    alignment_maps: BAMDirFmt,
    mode: str,
    common_args: list,
) -> (MultiFASTADirectoryFormat, dict):
    sample_set = _assert_samples(contigs, alignment_maps)

    bins = MultiFASTADirectoryFormat()
    for samp, props in sample_set.items():
        _process_sample(samp, props, mode, common_args, str(bins))

    if not glob.glob(os.path.join(str(bins), "*/*.fa")):
        raise ValueError(
            "No MAGs were formed during binning, please check your inputs."
        )

    contig_map = _generate_contig_map(bins)

    return bins, contig_map


def bin_contigs_semibin2(
    contigs: ContigSequencesDirFmt,
    alignment_maps: BAMDirFmt,
    # mode: str,
    training_type: str | None = None,
    orf_finder: str | None = None,
    environment: str | None = None,
    engine: str | None = None,
    sequencing_type: str | None = None,
    minfasta_kbs: int | None = None,
    no_recluster: bool | None = None,
    epochs: int | None = None,
    batch_size: int | None = None,
    max_node: int | None = None,
    max_edges: int | None = None,
    ratio: float | None = None,
    min_len: int | None = None,
    ml_threshold: int | None = None,
    threads: int | None = None,
    random_seed: int | None = None,
    debug: bool | None = None,
    verbose: bool | None = None,
) -> (MultiFASTADirectoryFormat, dict):
    kwargs = {
        k: v
        for k, v in locals().items()
        if k not in ["contigs", "alignment_maps", "mode", "training_type", "verbose"]
    }

    # mode = "single_easy_bin" if mode == "single" else "multi_easy_bin"
    mode = "single_easy_bin"

    common_args = _process_common_input_params(
        processing_func=_process_semibin2_arg, params=kwargs
    )

    return _bin_contigs_semibin2(
        contigs=contigs,
        alignment_maps=alignment_maps,
        mode=mode,
        common_args=common_args,
    )
