# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import json
import os
import tempfile
from copy import deepcopy
from importlib import resources
from shutil import copytree
from typing import List, Union

import pandas as pd
import q2templates

from q2_annotate.busco.utils import (
    _parse_busco_params,
    _parse_df_columns,
    _calculate_summary_stats,
    _get_feature_table,
    _cleanup_bootstrap,
    _validate_lineage_dataset_input,
    _extract_json_data,
    _process_busco_results,
    _add_unbinned_metrics,
)

from q2_annotate._utils import _process_common_input_params, run_command
from q2_annotate.busco.types import BuscoDatabaseDirFmt
from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_types.per_sample_sequences import (
    MultiMAGSequencesDirFmt,
    ContigSequencesDirFmt,
)
import warnings

TEMPLATES = resources.files("q2_annotate") / "assets" / "busco"


def _load_vega_spec(spec_name: str) -> dict:
    """Load a Vega spec from JSON file."""
    spec_path = TEMPLATES / "vega" / f"{spec_name}.json"
    return json.loads(spec_path.read_text())


def _prepare_histogram_data(results: pd.DataFrame) -> str:
    """Prepare melted data for histogram plots."""
    cols = [
        ["single", "duplicated", "fragmented", "missing", "complete"],
        ["completeness", "contamination", "contigs_n50", "length"],
    ]

    if not ("completeness" in results.columns and "contamination" in results.columns):
        cols[0].remove("completeness")
        cols[1].remove("contamination")

    melted = pd.melt(
        results,
        id_vars=["sample_id", "mag_id", "dataset", "n_markers"],
        value_vars=[*cols[0], *cols[1]],
        value_name="metric",
        var_name="category",
    )
    return json.dumps(melted.to_dict("records")).replace("NaN", "null")


def _prepare_box_plot_data(results: pd.DataFrame) -> dict:
    """Prepare data for box plots, one per metric."""
    metrics = [
        "single",
        "duplicated",
        "fragmented",
        "missing",
        "complete",
        "completeness",
        "contamination",
        "contigs_n50",
        "length",
    ]

    if not ("completeness" in results.columns and "contamination" in results.columns):
        metrics.remove("completeness")
        metrics.remove("contamination")

    data_by_metric = {}
    for metric in metrics:
        if metric in results.columns:
            df = results[["sample_id", "mag_id", metric]].copy()
            df.columns = ["sample_id", "mag_id", "value"]
            data_by_metric[metric] = df.to_dict("records")

    return data_by_metric


def _prepare_scatter_data(results: pd.DataFrame) -> tuple:
    """Prepare data for completeness vs contamination scatter plot."""
    if "completeness" not in results.columns or "contamination" not in results.columns:
        return None, False, 110, 110

    # Calculate axis bounds
    max_comp = pd.to_numeric(results["completeness"], errors="coerce").max(skipna=True)
    max_cont = pd.to_numeric(results["contamination"], errors="coerce").max(skipna=True)
    max_comp = 0 if pd.isna(max_comp) else float(max_comp)
    max_cont = 0 if pd.isna(max_cont) else float(max_cont)
    upper_x = max(5.0, min(101.0, round(max_comp * 1.01, 1)))
    upper_y = max(5.0, min(101.0, round(max_cont * 1.01, 1)))

    data = results.to_dict("records")
    return json.dumps(data).replace("NaN", "null"), True, upper_x, upper_y


def _prepare_detailed_data(results: pd.DataFrame) -> str:
    """Prepare melted data for detailed BUSCO plots."""
    busco_plot_data = pd.melt(
        results,
        id_vars=["sample_id", "mag_id", "dataset", "n_markers"],
        value_vars=["single", "duplicated", "fragmented", "missing"],
        value_name="BUSCO_percentage",
        var_name="category",
    )

    # Estimate fraction of sequences in each BUSCO category
    busco_plot_data["frac_markers"] = (
        "~"
        + (busco_plot_data["BUSCO_percentage"] * busco_plot_data["n_markers"] / 100)
        .round()
        .astype(int)
        .astype(str)
        + "/"
        + busco_plot_data["n_markers"].astype(str)
    )

    return json.dumps(busco_plot_data.to_dict("records")).replace("NaN", "null")


def _prepare_assembly_data(results: pd.DataFrame) -> str:
    """Prepare data for assembly metrics plots."""
    cols = [
        "sample_id",
        "mag_id",
        "scaffold_n50",
        "contigs_n50",
        "percent_gaps",
        "scaffolds",
    ]
    # Only include columns that exist
    cols = [c for c in cols if c in results.columns]
    data = results[cols].to_dict("records")
    return json.dumps(data).replace("NaN", "null")


def _run_busco(input_dir: str, output_dir: str, sample_id: str, params: List[str]):
    """Runs BUSCO on one (sample) directory

    Args:
        input_dir (str): Location where the MAG FASTA files are stored.
        output_dir (str): Location where the final results should be stored.
        sample_id (str): The sample ID.
        params (List[str]): List of parsed arguments to pass to BUSCO.
    """
    base_cmd = ["busco", *params]

    cmd = deepcopy(base_cmd)
    cmd.extend(["--in", input_dir, "--out_path", output_dir, "-o", sample_id])
    run_command(cmd, cwd=os.path.dirname(output_dir))


def _busco_helper(mags, common_args, additional_metrics):
    results_all = []
    # Get samples directories from MAGs
    if isinstance(mags, MultiMAGSequencesDirFmt):
        sample_dir = mags.sample_dict()
    elif isinstance(mags, MAGSequencesDirFmt):
        sample_dir = {"feature_data": mags.feature_dict()}

    with tempfile.TemporaryDirectory() as tmp:
        for sample_id, feature_dict in sample_dir.items():

            _run_busco(
                input_dir=os.path.join(
                    str(mags), "" if sample_id == "feature_data" else sample_id
                ),
                output_dir=str(tmp),
                sample_id=sample_id,
                params=common_args,
            )
            # Extract and process results from JSON files for one sample
            for mag_id, mag_fp in feature_dict.items():

                json_path = glob.glob(
                    os.path.join(
                        str(tmp), sample_id, os.path.basename(mag_fp), "*.json"
                    )
                )[0]

                results = _process_busco_results(
                    _extract_json_data(json_path),
                    sample_id,
                    mag_id,
                    os.path.basename(mag_fp),
                    additional_metrics,
                )
                results_all.append(results)

    return pd.DataFrame(results_all)


def _evaluate_busco(
    mags: Union[MultiMAGSequencesDirFmt, MAGSequencesDirFmt],
    db: BuscoDatabaseDirFmt,
    unbinned_contigs: ContigSequencesDirFmt = None,  # NEW unbinned
    mode: str = "genome",
    lineage_dataset: str = None,
    augustus: bool = False,
    augustus_parameters: str = None,
    augustus_species: str = None,
    cpu: int = 1,
    contig_break: int = 10,
    evalue: float = 1e-03,
    limit: int = 3,
    long: bool = False,
    metaeuk_parameters: str = None,
    metaeuk_rerun_parameters: str = None,
    miniprot: bool = False,
    additional_metrics: bool = False,
) -> pd.DataFrame:
    kwargs = {
        k: v
        for k, v in locals().items()
        if k not in ["mags", "db", "additional_metrics", "unbinned_contigs"]
    }
    kwargs["offline"] = True
    kwargs["download_path"] = str(db)

    if lineage_dataset is not None:
        _validate_lineage_dataset_input(
            lineage_dataset,
            db,
        )

    common_args = _process_common_input_params(
        processing_func=_parse_busco_params, params=kwargs
    )

    busco_results = _busco_helper(mags, common_args, additional_metrics)

    return busco_results


def _visualize_busco(
    output_dir: str,
    results: pd.DataFrame,
    unbinned_contigs: ContigSequencesDirFmt = None,
) -> None:
    results.to_csv(os.path.join(output_dir, "busco_results.csv"), index=False)

    # Add unbinned metrics if unbinned_contigs are provided
    if unbinned_contigs and "unbinned_contigs_count" not in results.columns:
        results = _add_unbinned_metrics(results, unbinned_contigs)

    results = _parse_df_columns(results)

    if len(results["sample_id"].unique()) >= 2:
        assets_subdir = "sample_data"
        tab_title = ["Sample details", "Feature details"]
        is_sample_data = True
    else:
        tab_title = ["Per-MAG metrics", "Feature details"]
        assets_subdir = "feature_data"
        is_sample_data = False

    templates = [
        str(TEMPLATES / assets_subdir / filename)
        for filename in ["index.html", "detailed_view.html", "table.html"]
    ]
    copytree(src=str(TEMPLATES / assets_subdir), dst=output_dir, dirs_exist_ok=True)
    for folder in ["css", "js", "vega"]:
        folder_dst = os.path.join(output_dir, folder)
        os.makedirs(folder_dst, exist_ok=True)
        copytree(
            src=str(TEMPLATES / folder),
            dst=folder_dst,
            dirs_exist_ok=True,
        )

    # Load Vega specs from JSON files
    vega_histogram_spec = json.dumps(_load_vega_spec("histogram"))
    vega_unbinned_spec = json.dumps(_load_vega_spec("unbinned"))
    vega_boxplot_spec = json.dumps(_load_vega_spec("box_plot"))
    vega_completeness_spec = json.dumps(_load_vega_spec("completeness"))
    vega_busco_detailed_spec = json.dumps(_load_vega_spec("busco_detailed"))
    vega_assembly_detailed_spec = json.dumps(_load_vega_spec("assembly_detailed"))

    # Prepare data for plots
    histogram_data = _prepare_histogram_data(results)
    box_plot_data = json.dumps(_prepare_box_plot_data(results)).replace("NaN", "null")
    scatter_data, comp_cont, upper_x, upper_y = _prepare_scatter_data(results)
    detailed_data = _prepare_detailed_data(results)
    assembly_data = _prepare_assembly_data(results)

    # Get sorted list of MAG IDs for consistent ordering
    mag_ids_sorted = json.dumps(sorted(results["mag_id"].unique().tolist()))

    # Provide sample IDs for coordinated filtering in templates
    sample_ids = []
    if is_sample_data:
        sample_ids = sorted([sid for sid in results["sample_id"].unique() if sid])

    # Prepare table and summary stats
    table_json = _get_feature_table(results)
    summary_stats_json = _calculate_summary_stats(results)

    # Check for unbinned contigs data
    unbinned, unbinned_data = False, None
    if (
        "unbinned_contigs" in results.columns
        and "unbinned_contigs_count" in results.columns
    ):
        unbinned = True
        unbinned_df = results.drop_duplicates(subset=["sample_id"])[
            ["sample_id", "unbinned_contigs_count"]
        ]
        unbinned_data = json.dumps(unbinned_df.to_dict("records")).replace(
            "NaN", "null"
        )

    # Available metrics for histograms/box plots
    metrics = ["single", "duplicated", "fragmented", "missing", "complete"]
    if "completeness" in results.columns:
        metrics.extend(["completeness", "contamination"])
    metrics.extend(["contigs_n50", "length"])

    # Assembly metrics for detailed view
    assembly_metrics = {
        "contigs_n50": "Contig N50 (bp)",
        "percent_gaps": "Percent Gaps (%)",
        "scaffolds": "Number of Scaffolds",
    }

    tabbed_context = {
        "tabs": [
            {"title": "QC overview", "url": "index.html"},
            {"title": tab_title[0], "url": "detailed_view.html"},
            {"title": tab_title[1], "url": "table.html"},
        ],
        # Vega specs
        "vega_histogram_spec": vega_histogram_spec,
        "vega_unbinned_spec": vega_unbinned_spec,
        "vega_box_plot_spec": vega_boxplot_spec,
        "vega_scatter_spec": vega_completeness_spec,
        "vega_busco_detailed_spec": vega_busco_detailed_spec,
        "vega_assembly_detailed_spec": vega_assembly_detailed_spec,
        # Data for plots
        "histogram_data": histogram_data,
        "box_plot_data": box_plot_data,
        "scatter_data": scatter_data,
        "detailed_data": detailed_data,
        "assembly_data": assembly_data,
        "mag_ids_sorted": mag_ids_sorted,
        # Metadata
        "metrics_json": json.dumps(metrics),
        "assembly_metrics_json": json.dumps(assembly_metrics),
        "sample_ids_json": json.dumps(sample_ids),
        "is_sample_data": is_sample_data,
        "comp_cont": comp_cont,
        "upper_x": upper_x,
        "upper_y": upper_y,
        "unbinned": unbinned,
        "unbinned_data": unbinned_data,
        # Table data
        "table": table_json,
        "summary_stats_json": summary_stats_json,
        "page_size": 100,
    }
    q2templates.render(templates, output_dir, context=tabbed_context)

    # Final cleanup, needed until we fully migrate to Bootstrap 5
    _cleanup_bootstrap(output_dir)


def evaluate_busco(
    ctx,
    mags,
    db,
    unbinned_contigs=None,
    mode="genome",
    lineage_dataset=None,
    augustus=False,
    augustus_parameters=None,
    augustus_species=None,
    cpu=1,
    contig_break=10,
    evalue=1e-03,
    limit=3,
    long=False,
    metaeuk_parameters=None,
    metaeuk_rerun_parameters=None,
    miniprot=False,
    additional_metrics=True,
    num_partitions=None,
):
    if not lineage_dataset:
        raise ValueError("'lineage-dataset' is required as a parameter")

    kwargs = {
        k: v
        for k, v in locals().items()
        if k not in ["mags", "unbinned_contigs", "ctx", "db", "num_partitions"]
    }
    _evaluate_busco = ctx.get_action("annotate", "_evaluate_busco")
    collate_busco_results = ctx.get_action("annotate", "collate_busco_results")
    _visualize_busco = ctx.get_action("annotate", "_visualize_busco")

    if issubclass(mags.format, MultiMAGSequencesDirFmt):
        partition_action = "partition_sample_data_mags"
    else:
        partition_action = "partition_feature_data_mags"
        if unbinned_contigs is not None:
            warnings.warn(
                "FeatureData[MAG] artifact was provided - "
                "unbinned contigs will be ignored."
            )

    partition_mags = ctx.get_action("types", partition_action)

    (partitioned_mags,) = partition_mags(mags, num_partitions)
    results = []

    # Run BUSCO evaluation on all partitions without unbinned contigs
    for mag_partition in partitioned_mags.values():
        (busco_result,) = _evaluate_busco(mag_partition, db, **kwargs)
        results.append(busco_result)

    (collated_results,) = collate_busco_results(results)

    # Pass unbinned_contigs to visualization step
    (visualization,) = _visualize_busco(collated_results, unbinned_contigs)

    return collated_results, visualization
