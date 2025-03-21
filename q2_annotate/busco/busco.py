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
    _draw_selectable_summary_histograms

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
        k: v for k, v in locals().items() if k not in ["bins", "busco_db"]
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
            json.dumps(_draw_selectable_summary_histograms(busco_results))
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
        if k not in ["bins", "ctx", "num_partitions"]
    }

    _evaluate_busco = ctx.get_action("annotate", "_evaluate_busco")
    collate_busco_results = ctx.get_action("annotate", "collate_busco_results")
    _visualize_busco = ctx.get_action("annotate", "_visualize_busco")

    if issubclass(bins.format, MultiMAGSequencesDirFmt):
        partition_action = "partition_sample_data_mags"
    else:
        partition_action = "partition_feature_data_mags"
    partition_mags = ctx.get_action("types", partition_action)

    (partitioned_mags, ) = partition_mags(bins, num_partitions)
    results = []
    for mag in partitioned_mags.values():
        (busco_result, ) = _evaluate_busco(mag, **kwargs)
        results.append(busco_result)

    collated_results, = collate_busco_results(results)
    visualization, = _visualize_busco(collated_results)

    return collated_results, visualization
