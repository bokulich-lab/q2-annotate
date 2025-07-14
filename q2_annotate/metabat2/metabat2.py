# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import re
from uuid import uuid4
from pathlib import Path

import os.path
import shutil
import tempfile
from copy import deepcopy

import skbio.io
from q2_types.feature_data import DNAIterator

from q2_types.per_sample_sequences import ContigSequencesDirFmt, BAMDirFmt
from q2_types.per_sample_sequences import MultiFASTADirectoryFormat

from q2_annotate._utils import run_command, _process_common_input_params
from q2_annotate.metabat2.utils import _process_metabat2_arg


def _assert_samples(contigs, maps) -> dict:
    contig_fps = contigs.sample_dict().values()
    map_fps = glob.glob(os.path.join(str(maps), "*.bam"))
    contig_fps, map_fps = sorted(contig_fps), sorted(map_fps)

    contig_samples = contigs.sample_dict().keys()
    map_samples = [Path(fp).stem.rsplit("_alignment", 1)[0] for fp in map_fps]
    if set(contig_samples) != set(map_samples):
        raise Exception(
            "Contigs and alignment maps should belong to the same sample set. "
            f'You provided contigs for samples: {",".join(contig_samples)} '
            f'but maps for samples: {",".join(map_samples)}. Please check '
            "your inputs and try again."
        )

    return {
        s: {"contigs": contig_fps[i], "map": map_fps[i]}
        for i, s in enumerate(contig_samples)
    }


def _sort_bams(samp_name, samp_props, loc):
    sorted_bam = os.path.join(loc, f"{samp_name}_alignment_sorted.bam")
    new_props = deepcopy(samp_props)
    run_command(["samtools", "sort", new_props["map"], "-o", sorted_bam], verbose=True)
    new_props["map"] = sorted_bam
    return new_props


def _estimate_depth(samp_name, samp_props, loc):
    depth_fp = os.path.join(str(loc), f"{samp_name}_depth.txt")
    run_command(
        [
            "jgi_summarize_bam_contig_depths",
            "--outputDepth",
            depth_fp,
            samp_props["map"],
        ],
        verbose=True,
    )
    return depth_fp


def _run_metabat2(samp_name, samp_props, loc, depth_fp, common_args):
    bins_dp = os.path.join(loc, samp_name)
    bins_prefix = os.path.join(bins_dp, "bin")
    os.makedirs(bins_dp)
    cmd = [
        "metabat2",
        "-i",
        samp_props["contigs"],
        "-a",
        depth_fp,
        "-o",
        bins_prefix,
        "--unbinned",
    ]
    cmd.extend(common_args)
    run_command(cmd, verbose=True)
    return bins_dp


def _process_sample(samp_name, samp_props, common_args, result_loc, unbinned_loc):
    with tempfile.TemporaryDirectory() as tmp:
        # sort alignment map
        props = _sort_bams(samp_name, samp_props, tmp)

        # calculate depth
        depth_fp = _estimate_depth(samp_name, props, tmp)

        # run metabat2
        bins_dp = _run_metabat2(samp_name, props, tmp, depth_fp, common_args)

        all_outputs = glob.glob(os.path.join(bins_dp, "*.fa"))
        all_bins = [x for x in all_outputs if re.search(r"bin\.[0-9]+\.fa$", x)]
        unbinned_fp = os.path.join(bins_dp, "bin.unbinned.fa")

        # rename using UUID v4
        bin_dest_dir = os.path.join(str(result_loc), samp_name)
        os.makedirs(bin_dest_dir, exist_ok=True)
        for old_bin in all_bins:
            new_bin = os.path.join(bin_dest_dir, f"{uuid4()}.fa")
            shutil.move(old_bin, new_bin)

        # move unbinned contigs
        unbinned_dest = os.path.join(str(unbinned_loc), f"{samp_name}_contigs.fa")
        if os.path.isfile(unbinned_fp):
            shutil.move(unbinned_fp, unbinned_dest)


def _generate_contig_map(bins: MultiFASTADirectoryFormat) -> dict:
    contig_map = {}
    for bin_fp, _ in bins.sequences.iter_views(DNAIterator):
        # bin_fp will look like /path/to/some/where/uuid4-bin-name.fa
        bin_id = os.path.splitext(os.path.basename(bin_fp))[0]
        seqs = skbio.read(
            os.path.join(str(bins), str(bin_fp)), format="fasta", verify=False
        )
        contigs = [x.metadata["id"] for x in seqs]
        contig_map[bin_id] = contigs
    return contig_map


def _bin_contigs_metabat(
    contigs: ContigSequencesDirFmt, alignment_maps: BAMDirFmt, common_args: list
) -> (MultiFASTADirectoryFormat, dict, ContigSequencesDirFmt):
    sample_set = _assert_samples(contigs, alignment_maps)

    bins = MultiFASTADirectoryFormat()
    unbinned = ContigSequencesDirFmt()
    for samp, props in sample_set.items():
        _process_sample(samp, props, common_args, str(bins), str(unbinned))

    if not glob.glob(os.path.join(str(bins), "*/*.fa")):
        raise ValueError(
            "No MAGs were formed during binning, please check your inputs."
        )

    contig_map = _generate_contig_map(bins)

    return bins, contig_map, unbinned


def bin_contigs_metabat(
    contigs: ContigSequencesDirFmt,
    alignment_maps: BAMDirFmt,
    min_contig: int = None,
    max_p: int = None,
    min_s: int = None,
    max_edges: int = None,
    p_tnf: int = None,
    no_add: bool = None,
    min_cv: int = None,
    min_cv_sum: int = None,
    min_cls_size: int = None,
    num_threads: int = None,
    seed: int = None,
    debug: bool = None,
    verbose: bool = None,
) -> (MultiFASTADirectoryFormat, dict, ContigSequencesDirFmt):

    kwargs = {
        k: v for k, v in locals().items() if k not in ["contigs", "alignment_maps"]
    }
    common_args = _process_common_input_params(
        processing_func=_process_metabat2_arg, params=kwargs
    )

    return _bin_contigs_metabat(
        contigs=contigs, alignment_maps=alignment_maps, common_args=common_args
    )
