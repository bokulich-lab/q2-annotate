# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json

import os
import tempfile
from copy import deepcopy
from shutil import copytree
from typing import List, Dict, Union

import pandas as pd
import q2templates

from q2_annotate.busco.plots_detailed import _draw_detailed_plots
from q2_annotate.busco.plots_summary import _draw_marker_summary_histograms, \
    _draw_selectable_summary_histograms, _draw_selectable_unbinned_histograms

from q2_annotate.busco.utils import (
    _parse_busco_params, _collect_summaries, _rename_columns,
    _parse_df_columns, _partition_dataframe, _calculate_summary_stats,
    _get_feature_table, _cleanup_bootstrap, _get_mag_lengths,
    _validate_lineage_dataset_input
)
from q2_annotate._utils import _process_common_input_params, run_command
from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt
from q2_annotate.busco.types import BuscoDatabaseDirFmt
from q2_types.feature_data_mag import MAGSequencesDirFmt
### NEW
from q2_types.feature_data import DNAIterator
from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt, ContigSequencesDirFmt

# from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt, ContigSequencesDirFmt
from q2_types.per_sample_sequences import MultiFASTADirectoryFormat

#####for partition
# import shutil
import numpy as np
# import pandas as pd
from qiime2.util import duplicate
# from q2_types._util import _validate_num_partitions
### new for filter_contigs
from q2_assembly.filter import filter_contigs
from qiime2 import Metadata
import warnings
from skbio import DNA
from skbio.io import read
from typing import List
from pathlib import Path
from qiime2.core.type import SemanticType
from qiime2 import Artifact
#####
def filter_unbinned_for_partition(unbinned_contigs, mag_partition, _filter_contigs):
    """
    Filters the unbinned contigs to match the sample IDs in a MAG partition.

    Args:
        unbinned_contigs (ContigSequencesDirFmt): The full unbinned contigs.
        mag_partition (MultiMAGSequencesDirFmt): One partition of MAGs.
        _filter_contigs (Action): QIIME 2 action to filter contigs.

    Returns:
        ContigSequencesDirFmt: Filtered unbinned contigs matching the partition samples.
    """

    # sample_ids = list(mag_partition.sample_dict().keys())
    sample_ids = list(mag_partition.view(MultiMAGSequencesDirFmt).sample_dict().keys())
    metadata = Metadata(pd.DataFrame(index=pd.Index(sample_ids, name="ID")))
    id_list = ", ".join([f"'{sid}'" for sid in sample_ids])
    where = f"ID IN ({id_list})"
    (filtered_unbinned,) = _filter_contigs(
        contigs=unbinned_contigs,
        metadata=metadata,
        where=where
    )
    return filtered_unbinned
def get_fasta_files_from_dir(directory: Path) -> list:
    # Only match common FASTA extensions starting with '.fa' or '.fna'
    return [f for f in directory.glob('*') if f.suffix in {'.fa', '.fasta'}]
# def count_contigs(dirfmt) -> int:
#     """
#     Count the number of DNA sequences in a QIIME 2 directory format.

#     Parameters
#     ----------
#     dirfmt : DirectoryFormat
#         A QIIME 2 directory format containing one or more FASTA files.
    
#     Returns
#     -------
#     int
#         The total number of sequences.
#     """
    
#     total_sequences = 0
#     #accepts a list of file paths 
#     # instead of  DNAIterator for seq in read(str(fp), format="fasta", constructor=DNA):
#     #Iterate over sequences directly using QIIME 2's DNAIterator
#     for fasta_fp, dna_iterator in dirfmt.sequences.iter_views(DNAIterator): #one element
#         # print(f"{fasta_fp} contains {sequence_count} sequences")  # Debugging output
#         total_sequences += sum(1 for _ in dna_iterat

def count_contigs(file_paths: List[Path]) -> int:
    """
    Count the number of DNA sequences across a list of FASTA files.

    Parameters
    ----------
    file_paths : list of Path
        List of FASTA file paths (.fa, .fasta, .fna).

    Returns
    -------
    int
        Total number of sequences across all files.
    """
    total_sequences = 0

    for fp in file_paths:
        total_sequences += sum(1 for _ in read(str(fp), format="fasta", constructor=DNA))

    return total_sequences
# def calculate_unbinned_percentage(all_binned_mags, sample_unbinned_contigs) -> tuple[float, int]:
#     """
#     Calculate the percentage and absolute count of unbinned contigs for a single sample.

#     Parameters
#     ----------
#     all_binned_mags : MultiMAGSequencesDirFmt
#         Binned contigs (MAGs) from all samples.
    
#     sample_unbinned_contigs : ContigSequencesDirFmt
#         Unbinned contigs from one specific sample.

#     Returns
#     -------
#     percentage_unbinned : float
#         The percentage of unbinned contigs relative to the total number of contigs 
#         (binned + unbinned) for this sample.

#     unbinned_contigs_count : int
#         The number of unbinned contigs in this sample.
#     """
#     # Count sequences
#     # binned_contigs = count_binned_contigs(all_binned_mags)
#     unbinned_contigs_count = count_contigs(sample_unbinned_contigs)

#     # Calculate percentage
#     total = all_binned_mags + unbinned_contigs_count
#     percentage_unbinned = (unbinned_contigs_count / total) * 100 if total > 0 else 0

#     return percentage_unbinned, unbinned_contigs_count

def calculate_unbinned_percentage(mags_per_sample, unbinned_contigs_per_sample) -> tuple[float, int]:
    """
    Calculate the percentage and absolute count of unbinned contigs for a single sample.

    Parameters
    ----------
    mags_per_sample : MultiMAGSequencesDirFmt
        Binned contigs (MAGs) from one specific sample.
    
    sample_unbinned_contigs : ContigSequencesDirFmt
        Unbinned contigs from one specific sample.

    Returns
    -------
    percentage_unbinned : float
        The percentage of unbinned contigs relative to the total number of contigs 
        (binned + unbinned) for this sample.

    unbinned_contigs_count : int
        The number of unbinned contigs in this sample.
    """
    #change input
    # Count sequences
    binned_contigs = count_contigs(mags_per_sample)
    unbinned_contigs_count = count_contigs(unbinned_contigs_per_sample)

    # Calculate percentage
    total = binned_contigs + unbinned_contigs_count
    percentage_unbinned = (unbinned_contigs_count / total) * 100 if total > 0 else 0

    return percentage_unbinned, unbinned_contigs_count
### NEW ENDS

def _run_busco(
    output_dir: str,
    mags: Union[MultiMAGSequencesDirFmt, MAGSequencesDirFmt],
    params: List[str]
) -> Dict[str, str]:
    """Evaluates bins for all samples using BUSCO.

    Args:
        output_dir (str): Location where the final results should be stored.
        mags (MultiMAGSequencesDirFmt): The mags to be analyzed.
        params (List[str]): List of parsed arguments to pass to BUSCO.

    Returns:
        dict: Dictionary where keys are sample IDs and values are the paths
            to the `batch_summary.txt` generated by BUSCO, e.g.
            `tmp/busco_output/<sample_id>/batch_summary.txt`.
    """
    base_cmd = ["busco", *params]
    # Get samples directories from MAGs
    if isinstance(mags, MultiMAGSequencesDirFmt):
        manifest: pd.DataFrame = mags.manifest.view(pd.DataFrame)
        manifest["sample_dir"] = manifest.filename.apply(
            lambda x: os.path.dirname(x)
        )
        sample_dirs = manifest["sample_dir"].unique()

    elif isinstance(mags, MAGSequencesDirFmt):
        sample_dirs = [str(mags)]

    path_to_run_summaries = {}

    # For every unique sample dir run busco
    for sample_dir in sample_dirs:
        sample = os.path.split(sample_dir)[-1]

        cmd = deepcopy(base_cmd)
        cmd.extend([
            "--in",
            sample_dir,
            "--out_path",
            output_dir,
            "-o",
            sample
        ])
        run_command(cmd,  cwd=os.path.dirname(output_dir))

        path_to_run_summary = os.path.join(
            output_dir, sample, "batch_summary.txt"
        )
        if os.path.isfile(path_to_run_summary):
            path_to_run_summaries[sample] = path_to_run_summary
        else:
            raise FileNotFoundError(
                f"BUSCO batch summary file {path_to_run_summary} not found."
            )

    return path_to_run_summaries


def _busco_helper(bins, common_args):
    with tempfile.TemporaryDirectory() as tmp:
        path_to_run_summaries = _run_busco(
            output_dir=os.path.join(tmp, "busco_output"),
            mags=bins,
            params=common_args,
        )
        # Returns pd.Dataframe
        all_summaries = _collect_summaries(
            run_summaries_fp_map=path_to_run_summaries,
        )
    all_summaries = _rename_columns(all_summaries)

    lengths = _get_mag_lengths(bins)
    all_summaries = all_summaries.join(lengths, on="mag_id")

    return all_summaries


def _evaluate_busco(
    bins: Union[MultiMAGSequencesDirFmt, MAGSequencesDirFmt],
    busco_db: BuscoDatabaseDirFmt,
    unbinned_contigs: ContigSequencesDirFmt = None, ### NEW unbinned
    mode: str = "genome",
    lineage_dataset: str = None,
    augustus: bool = False,
    augustus_parameters: str = None,
    augustus_species: str = None,
    auto_lineage: bool = False,
    auto_lineage_euk: bool = False,
    auto_lineage_prok: bool = False,
    cpu: int = 1,
    config: str = None,
    contig_break: int = 10,
    evalue: float = 1e-03,
    force: bool = False,
    limit: int = 3,
    long: bool = False,
    metaeuk_parameters: str = None,
    metaeuk_rerun_parameters: str = None,
    miniprot: bool = False,
    scaffold_composition: bool = False,
) -> pd.DataFrame:
    kwargs = {
        k: v for k, v in locals().items() if k not in ["bins", "unbinned_contigs", "busco_db"]#exclude unbinned
    }
    kwargs["offline"] = True
    kwargs["download_path"] = f"{str(busco_db)}/busco_downloads"

    if lineage_dataset is not None:
        _validate_lineage_dataset_input(
            lineage_dataset, auto_lineage, auto_lineage_euk, auto_lineage_prok,
            busco_db, kwargs  # kwargs may be modified inside this function
        )

    # Filter out all kwargs that are None, False or 0.0
    common_args = _process_common_input_params(
        processing_func=_parse_busco_params, params=kwargs
    )
    print("xixo")
    # if isinstance(bins.view(MultiMAGSequencesDirFmt), MultiMAGSequencesDirFmt):
    if isinstance(bins, MultiMAGSequencesDirFmt):
    # if issubclass(type(bins), MultiMAGSequencesDirFmt):
        print("xoxo")
        busco_results = _busco_helper(bins, common_args)
        busco_results["unbinned_contigs"] = pd.NA
        busco_results["unbinned_contigs_count"] = pd.NA
        # binned_contigs = count_contigs(bins)
        for unbinned_id, unbinned_path in unbinned_contigs.sample_dict().items(): 
            # percentage, count =  calculate_unbinned_percentage(bins, unbinned_values)
            # binned_sample_dir = MultiMAGSequencesDirFmt(path=bins.path / unbinned_id, mode='r')
            binned_dir = bins.path / unbinned_id
            binned_fasta_paths = get_fasta_files_from_dir(binned_dir)
            # unbinned_fasta_paths = get_fasta_files_from_dir(unbinned_path)
            # percentage, count = calculate_unbinned_percentage(MultiMAGSequencesDirFmt(path=bins.path / unbinned_id, mode='r'), ContigSequencesDirFmt(unbinned_path, mode='r'))
            percentage, count = calculate_unbinned_percentage(binned_fasta_paths, [unbinned_path])
            print(percentage)
            busco_results.loc[busco_results["sample_id"] == unbinned_id, "unbinned_contigs"] = float(percentage)
            busco_results.loc[busco_results["sample_id"] == unbinned_id, "unbinned_contigs_count"] = int(count)
            
            # busco_results.loc[busco_results["sample_id"] == unbinned_id, "unbinned_contigs"] = percentage
            # busco_results.loc[busco_results["sample_id"] == unbinned_id, "unbinned_contigs_count"] = count
            # ## loc -no dictionary
        return busco_results 
    
    return _busco_helper(bins, common_args)


def _visualize_busco(output_dir: str, busco_results: pd.DataFrame) -> None:
    busco_results.to_csv(
        os.path.join(output_dir, "busco_results.csv"),
        index=False
    )
    # Outputs different df for sample and feature data
    busco_results = _parse_df_columns(busco_results)
    max_rows = 100

    # Partition data frames
    if len(busco_results["sample_id"].unique()) >= 2:
        counter_col = "sample_id"
        assets_subdir = "sample_data"
        tab_title = ["Sample details", "Feature details"]
        is_sample_data = True

        # Draw selectable histograms (only for sample data mags)
        tabbed_context = {
            "vega_summary_selectable_json":
            json.dumps(_draw_selectable_summary_histograms(busco_results)),
            "vega_selectable_unbinned_json":
            json.dumps(_draw_selectable_unbinned_histograms(busco_results)),
            
        }
    else:
        counter_col = "mag_id"
        tab_title = ["BUSCO plots", "BUSCO table"]
        assets_subdir = "feature_data"
        is_sample_data = False
        tabbed_context = {}  # Init as empty bc we update it below

    dfs = _partition_dataframe(busco_results, max_rows, is_sample_data)

    # Copy BUSCO results from tmp dir to output_dir
    TEMPLATES = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        "assets",
        "busco"
    )
    templates = [
        os.path.join(TEMPLATES, assets_subdir, file_name)
        for file_name in ["index.html", "detailed_view.html", "table.html"]
    ]
    copytree(
        src=os.path.join(TEMPLATES, assets_subdir),
        dst=output_dir,
        dirs_exist_ok=True
    )
    for folder in ["css", "js"]:
        os.makedirs(os.path.join(output_dir, folder))
        copytree(
            src=os.path.join(TEMPLATES, folder),
            dst=os.path.join(output_dir, folder),
            dirs_exist_ok=True
        )

    # Partition data frames and draw detailed plots
    context = {}
    counter_left = 1
    for i, df in enumerate(dfs):
        count = df[counter_col].nunique()
        counter_right = counter_left + count - 1
        counters = {"from": counter_left, "to": counter_right}
        counter_left += count
        subcontext = _draw_detailed_plots(
            df,
            is_sample_data,
            width=600,
            height=30,
            title_font_size=20,
            label_font_size=17,
            spacing=20,
        )
        context.update({
            f"partition_{i}":
            {
                "subcontext": subcontext,
                "counters": counters,
                "ids": df[counter_col].unique().tolist(),
            }
        })

    # Render
    vega_json = json.dumps(context)
    vega_json_summary = json.dumps(
        _draw_marker_summary_histograms(busco_results)
    )
    table_json = _get_feature_table(busco_results)
    # print("BUSCO DF TYPES:\n", df.dtypes)
    # print("BUSCO DF unbinned_contigs values:\n", df["unbinned_contigs"].tolist())
    stats_json = _calculate_summary_stats(busco_results)
    tabbed_context.update({
        "tabs": [
            {"title": "QC overview", "url": "index.html"},
            {"title": tab_title[0], "url": "detailed_view.html"},
            {"title": tab_title[1], "url": "table.html"}
        ],
        "vega_json": vega_json,
        "vega_summary_json": vega_json_summary,
        "table": table_json,
        "summary_stats_json": stats_json,
        "page_size": 100
    })
    q2templates.render(templates, output_dir, context=tabbed_context)

    # Final cleanup, needed until we fully migrate to Bootstrap 5
    _cleanup_bootstrap(output_dir)


def evaluate_busco(
    ctx,
    bins,
    unbinned_contigs,
    busco_db,
    mode="genome",
    lineage_dataset=None,
    augustus=False,
    augustus_parameters=None,
    augustus_species=None,
    auto_lineage=False,
    auto_lineage_euk=False,
    auto_lineage_prok=False,
    cpu=1,
    config=None,
    contig_break=10,
    evalue=1e-03,
    force=False,
    limit=3,
    long=False,
    metaeuk_parameters=None,
    metaeuk_rerun_parameters=None,
    miniprot=False,
    scaffold_composition=False,
    num_partitions=None
):

    kwargs = {
        k: v for k, v in locals().items()
        if k not in ["bins", "unbinned_contigs", "ctx", "busco_db", "num_partitions"] # NEW exlude unbinned, why busco_db not excluded ??? "num_partitions"
    }
    _evaluate_busco = ctx.get_action("annotate", "_evaluate_busco")
    collate_busco_results = ctx.get_action("annotate", "collate_busco_results")
    _visualize_busco = ctx.get_action("annotate", "_visualize_busco")
    _filter_contigs = ctx.get_action("assembly", "filter_contigs")
    #add partition cases for unbinned this is possible/ makes sense to have unbinned only when we partition one mag per sample
    # print(bins.type)
    # Get samples directories from MAGs
    # print("TYPE:", type(bins))
    # print("IS INSTANCE OF Artifact?", isinstance(bins, Artifact))
    # print("SEMANTIC TYPE:", bins.type)
    # print("FORMAT:", bins.format)
    # print(str(bins.type))
    if str(bins.type) == 'SampleData[MAGs]':
    # if isinstance(bins.view(MultiMAGSequencesDirFmt), MultiMAGSequencesDirFmt):
    # if bins.type == SemanticType('SampleData[MAGs]'):
    # if issubclass(type(bins), MultiMAGSequencesDirFmt):
        partition_action = "partition_sample_data_mags"
    else:
        partition_action = "partition_feature_data_mags" #here unbined does not make sense! CREATE A WARNING
        warnings.warn("FeatureData[MAG] artifact was provided - unbinned contigs will be ignored.")

    partition_mags = ctx.get_action("types", partition_action)
    #####NEW 
    # if unbinned_action is not None:
    #     partition_unbinned = ctx.get_action("types", unbinned_action)
    #     (partitioned_unbinned,) = partition_unbinned(unbinned_contigs, num_partitions)
    # else:
        
    #     warnings.warn("Unbinned contigs will not be partitioned because partitioning by feature is active.")
    #     partitioned_unbinned = None  # Unused later, or handled conditionally   
    #####NEW ENDS 

    (partitioned_mags, ) = partition_mags(bins, num_partitions)
   
    #####new go through mags but also partition unbinned contigs
    results = []
    # print("TYPE:", type(bins))
    # print("IS INSTANCE OF Artifact?", isinstance(bins, Artifact))
    if str(bins.type) == 'SampleData[MAGs]':
    # if isinstance(bins.view(MultiMAGSequencesDirFmt), MultiMAGSequencesDirFmt):
    # if isinstance(bins, MultiMAGSequencesDirFmt):
    # if bins.type == SemanticType('SampleData[MAGs]'):
    # if issubclass(type(bins), MultiMAGSequencesDirFmt):
        for partition_id, mag_partition in partitioned_mags.items():
            # # Get the sample IDs in this partition
            # sample_ids = list(mag_partition.view(MultiMAGSequencesDirFmt).sample_dict().keys())
            # # Filter the unbinned contigs for this partition
            # metadata = Metadata(pd.DataFrame(index=pd.Index(sample_ids, name="ID")))
            # id_list = ", ".join([f"'{sid}'" for sid in sample_ids])
            # where = f"ID IN ({id_list})"
            # filtered_unbinned = _filter_contigs(
            #     contigs=unbinned_contigs,
            #     metadata=metadata,
            #     where=where
            # )
            filtered_unbinned = filter_unbinned_for_partition(unbinned_contigs, mag_partition, _filter_contigs)
            # filtered_unbinned, = _filter_contigs(
            #     contigs=unbinned_contigs, metadata=metadata
            # )
            # Run BUSCO for this partition of MAGs (with filtered unbinned if needed)
            # If you're calculating unbinned percentages inside `_evaluate_busco`, pass it in
            print("FILTERED UNBINNED TYPE:", type(filtered_unbinned))

            (busco_result,) = _evaluate_busco(mag_partition, busco_db, filtered_unbinned, **kwargs)
            # df = busco_result.view(pd.DataFrame)
            # columns = df.columns

            # # Extract the first and last two column names
            # selected_columns = [columns[0], columns[-2], columns[-1]]

            # # Subset the DataFrame using those columns
            # subset_df = df[selected_columns]
            # print(subset_df.to_string(index=False))
            results.append(busco_result)
    else:
        for mag in partitioned_mags.values(): # each mag is a subset of bins
            (busco_result, ) = _evaluate_busco(mag, busco_db, unbinned_contigs, **kwargs)
            results.append(busco_result)
    for i, result in enumerate(results):
        print(f"\nResult {i}: Type = {type(result)}")
        if isinstance(result, Artifact):
            try:
                df = result.view(pd.DataFrame)
                print(f"Artifact {i} DataFrame:\n", df.head())
            except Exception as e:
                print(f"Failed to view Artifact {i} as DataFrame: {e}")
        else:
            print("Not a QIIME2 Artifact.")
    collated_results, = collate_busco_results(results)
    visualization, = _visualize_busco(collated_results)

    return collated_results, visualization
