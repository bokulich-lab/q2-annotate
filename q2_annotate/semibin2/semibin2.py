# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import gzip
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
from q2_annotate.semibin2.utils import _process_semibin2_arg


def _assert_samples(contigs, maps) -> dict:
    """Validate that contigs and alignment maps belong to the same sample set."""
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
    """Sort BAM files for SemiBin2 processing."""
    sorted_bam = os.path.join(loc, f"{samp_name}_alignment_sorted.bam")
    new_props = deepcopy(samp_props)
    run_command(["samtools", "sort", new_props["map"], "-o", sorted_bam], verbose=True)
    new_props["map"] = sorted_bam
    return new_props


def _run_semibin2_single(samp_name, samp_props, loc, common_args, multi_sample=False):
    """Run SemiBin2 in single-sample mode."""
    bins_dp = os.path.join(loc, samp_name)
    os.makedirs(bins_dp, exist_ok=True)
    
    cmd = [
        "SemiBin2",
        "single_easy_bin" if not multi_sample else "multi_easy_bin",
        "-i", samp_props["contigs"],
        "-b", samp_props["map"],
        "-o", bins_dp,
        "--verbose",  # Always run in verbose mode as specified in requirements
    ]
    cmd.extend(common_args)
    run_command(cmd, verbose=True)
    return bins_dp


def _concatenate_contigs_with_semibin2(sample_set, loc):
    """Concatenate contigs using SemiBin2's concatenate_fasta command."""
    # Collect all contig file paths
    contig_files = []
    for samp_name, props in sample_set.items():
        if os.path.exists(props["contigs"]):
            contig_files.append(props["contigs"])
    
    # Use SemiBin2 concatenate_fasta to combine and rename contigs
    combined_contigs = os.path.join(loc, "combined_contigs.fa")
    cmd = [
        "SemiBin2",
        "concatenate_fasta",
        "-i",
    ]
    cmd.extend(contig_files)  # Add all contig file paths
    cmd.extend(["-o", combined_contigs])
    run_command(cmd, verbose=True)
    return combined_contigs


def _run_semibin2_multi(sample_set, loc, common_args):
    """Run SemiBin2 in multi-sample mode."""
    bins_dp = os.path.join(loc, "multi_sample_bins")
    os.makedirs(bins_dp, exist_ok=True)
    
    # Use SemiBin2's concatenate_fasta command to combine contigs with proper naming
    combined_contigs = _concatenate_contigs_with_semibin2(sample_set, loc)
    
    # Collect all BAM files
    bam_files = [props["map"] for props in sample_set.values()]
    
    cmd = [
        "SemiBin2",
        "multi_easy_bin",
        "-i", combined_contigs,
        "-b", *bam_files,
        "-o", bins_dp,
        "--verbose",  # Always run in verbose mode as specified in requirements
    ]
    cmd.extend(common_args)
    run_command(cmd, verbose=True)
    return bins_dp


def _process_semibin2_outputs(bins_dp, samp_name, result_loc):
    """Process SemiBin2 outputs and organize them like MetaBAT2."""
    # SemiBin2 outputs bins to output_recluster_bins/ directory
    semibin_bins_dir = os.path.join(bins_dp, "output_recluster_bins")
    
    if not os.path.exists(semibin_bins_dir):
        # Fallback to output_bins if recluster_bins doesn't exist
        semibin_bins_dir = os.path.join(bins_dp, "output_bins")
        print(f"WARNING: output_recluster_bins/ directory not found, using output_bins/ instead.")
    
    if not os.path.exists(semibin_bins_dir):
        print(f"WARNING: No bins found for sample {samp_name}.")
        return  # No bins generated
    
    all_outputs = glob.glob(os.path.join(semibin_bins_dir, "*.fa.gz"))
    
    # Rename bins using UUID v4 (following MetaBAT2 pattern)
    bin_dest_dir = os.path.join(str(result_loc), samp_name)
    os.makedirs(bin_dest_dir, exist_ok=True)
    for old_bin in all_outputs:
        print(f"Moving {old_bin} to {bin_dest_dir}")
        new_bin = os.path.join(bin_dest_dir, f"{uuid4()}.fa")
        with gzip.open(old_bin, 'rb') as f_in:
            with open(new_bin, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)


def _process_sample(samp_name, samp_props, common_args, result_loc, multi_sample=False):
    """Process a single sample through the SemiBin2 pipeline."""
    with tempfile.TemporaryDirectory() as tmp:
        # Sort alignment map
        props = _sort_bams(samp_name, samp_props, tmp)

        # Run SemiBin2
        bins_dp = _run_semibin2_single(samp_name, props, tmp, common_args, multi_sample)

        # Process outputs
        _process_semibin2_outputs(bins_dp, samp_name, result_loc)


def _generate_contig_map(bins: MultiFASTADirectoryFormat) -> dict:
    """Generate a mapping from bin IDs to contig IDs."""
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


def _bin_contigs_semibin2(
    contigs: ContigSequencesDirFmt, 
    alignment_maps: BAMDirFmt, 
    common_args: list,
    multi_sample: bool = False
) -> (MultiFASTADirectoryFormat, dict):
    """Main function to bin contigs using SemiBin2."""
    sample_set = _assert_samples(contigs, alignment_maps)

    bins = MultiFASTADirectoryFormat()
    
    if multi_sample and len(sample_set) > 1:
        # Run multi-sample binning
        with tempfile.TemporaryDirectory() as tmp:
            # Sort all BAMs first
            sorted_sample_set = {}
            for samp, props in sample_set.items():
                sorted_sample_set[samp] = _sort_bams(samp, props, tmp)
            
            bins_dp = _run_semibin2_multi(sorted_sample_set, tmp, common_args)
            
            # Process outputs for multi-sample mode
            # In multi-sample mode, we need to distribute bins back to samples
            # This is more complex and may require parsing SemiBin2's contig-to-bin mapping
            for samp in sample_set:
                _process_semibin2_outputs(bins_dp, samp, str(bins))
    else:
        # Run single-sample binning for each sample
        for samp, props in sample_set.items():
            _process_sample(samp, props, common_args, str(bins), False)

    if not glob.glob(os.path.join(str(bins), "*/*.fa")):
        raise ValueError(
            "No MAGs were formed during binning, please check your inputs."
        )

    contig_map = _generate_contig_map(bins)

    return bins, contig_map


def bin_contigs_semibin2(
    contigs: ContigSequencesDirFmt,
    alignment_maps: BAMDirFmt,
    multi_sample: bool = False,
    threads: int = None,
    min_len: int = None,
    batch_size: int = None,
    epochs: int = None,
    random_seed: int = None,
    sequencing_type: str = None,
) -> (MultiFASTADirectoryFormat, dict):
    """
    Bin contigs into MAGs using SemiBin2.
    
    Parameters
    ----------
    contigs : ContigSequencesDirFmt
        Contigs to be binned.
    alignment_maps : BAMDirFmt
        Reads-to-contig alignment maps.
    multi_sample : bool, optional
        Whether to perform multi-sample binning. Default is False (single-sample).
    threads : int, optional
        Number of threads to use.
    min_len : int, optional
        Minimum length of contigs for binning.
    batch_size : int, optional
        Batch size for neural network training.
    epochs : int, optional
        Number of epochs for neural network training.
    random_seed : int, optional
        Random seed for reproducibility.
    sequencing_type : str, optional
        Sequencing data type (e.g., 'human_gut', 'ocean', 'soil', etc.).
        
    Returns
    -------
    bins : MultiFASTADirectoryFormat
        The resulting MAGs.
    contig_map : dict
        Mapping of MAG identifiers to contig identifiers.
    """
    kwargs = {
        k: v for k, v in locals().items() 
        if k not in ["contigs", "alignment_maps", "multi_sample"]
    }
    common_args = _process_common_input_params(
        processing_func=_process_semibin2_arg, params=kwargs
    )

    return _bin_contigs_semibin2(
        contigs=contigs, 
        alignment_maps=alignment_maps, 
        common_args=common_args,
        multi_sample=multi_sample
    )