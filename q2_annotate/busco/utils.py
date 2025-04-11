# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json
import os
import warnings
import pandas as pd
from typing import List, Union
import skbio.io
from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt
from q2_annotate.busco.types import BuscoDatabaseDirFmt
from q2_types.feature_data_mag import MAGSequencesDirFmt

arguments_with_hyphens = {
    "auto_lineage": "auto-lineage",
    "auto_lineage_euk": "auto-lineage-euk",
    "auto_lineage_prok": "auto-lineage-prok",
    "list_datasets": "list-datasets",
    "update_data": "update-data",
}

MARKER_COLS = ["single", "duplicated", "fragmented", "missing", "complete"]


def _validate_lineage_dataset_input(
        lineage_dataset: str,
        auto_lineage: bool,
        auto_lineage_euk: bool,
        auto_lineage_prok: bool,
        busco_db: BuscoDatabaseDirFmt,
        kwargs
        ) -> None:
    # When lineage_dataset is specified all other lineage flags are ignored
    if any([auto_lineage, auto_lineage_euk, auto_lineage_prok]):
        warnings.warn(
            f"`--p-lineage-dataset` was specified as '{lineage_dataset}'. "
            "`--p-auto-lineage` flags will be ignored."
        )
        kwargs["auto_lineage"] = False
        kwargs["auto_lineage_euk"] = False
        kwargs["auto_lineage_prok"] = False

    # Check that lineage in deed exits inside Busco DB (if provided)
    if busco_db is not None:
        if not os.path.exists(
            f"{str(busco_db)}/busco_downloads/lineages/{lineage_dataset}"
        ):
            present_lineages = os.listdir(
                os.path.join(str(busco_db), "busco_downloads/lineages/")
            )
            raise ValueError(
                f"The specified lineage_dataset ({lineage_dataset}) "
                "is not present in input database. "
                "The datasets present in the input database are: "
                f"{present_lineages}"
            )


def _parse_busco_params(arg_key, arg_val) -> List[str]:
    """Creates a list with argument and its value to be consumed by BUSCO.
    Argument names will be converted to command line parameters by
    appending a '--' prefix and in some cases replacing "_" for "-"
    (only for e.g. `arguments_with_hyphens`)

    Args:
        arg_key (str): Argument name.
        arg_val: Argument value.
    Returns:
        [converted_arg, arg_value]: List containing a prepared command line
            parameter and, optionally, its value.
    """

    # If the key is in arguments_with_hyphens, modify key
    if arg_key in arguments_with_hyphens.keys():
        arg_key = arguments_with_hyphens[arg_key]

    if isinstance(arg_val, bool):
        return [f"--{arg_key}"]
    else:
        return [f"--{arg_key}", str(arg_val)]


def _partition_dataframe(
    df: pd.DataFrame, max_rows: int, is_sample_data: bool
) -> list:
    """
    Partitions a DataFrame into smaller DataFrames based on
    a maximum row limit.

    If is_sample_data = True:
    This function groups the DataFrame by 'sample_id' and then partitions
    these groups into smaller DataFrames. Each partition will have a total
    row count less than or equal to the max_rows parameter (unless a single
    partition exceeds the max_rows, in which case it will have all the
    MAGs included). The last group in a partition can exceed the max_rows
    limit.

    If is_sample_data = False:
    Partitions a DataFrame into smaller DataFrames based on
    a maximum row limit. Each partition will have a total
    row count less than or equal to the `max_rows` parameter.

    Args:
        df (pd.DataFrame): The DataFrame to partition. It should have a
            'sample_id' column.
        max_rows (int): The maximum number of rows that each partitioned
            DataFrame should have.

    Returns:
        list: A list of partitioned DataFrames. Each DataFrame in the
            list is a partition of the original DataFrame.
    """
    if is_sample_data:
        groups = [group for _, group in df.groupby('sample_id')]
        partitions = []
        temp = []
        total_rows = 0

        for group in groups:
            if total_rows + len(group) > max_rows:
                if temp:
                    partitions.append(pd.concat(temp))
                temp = [group]
                total_rows = len(group)
            else:
                temp.append(group)
                total_rows += len(group)

        if temp:
            partitions.append(pd.concat(temp))

        return partitions
    else:
        return [df[i:i+max_rows] for i in range(0, len(df), max_rows)]


def _collect_summaries(run_summaries_fp_map: dict) -> pd.DataFrame:
    """
    Reads-in the sample-wise summaries and concatenates them in one
    pd.DataFrame, which is saved to file.

    Args:
        run_summaries_fp_map (dict): dict where key is sample id
            and value is path "tmp/sample_id/batch_summary.txt"

    Returns:
        all_summaries (pd.DataFrame): DataFrame composed of the individual
            run summaries.
    """

    all_summaries = []
    for sample_id, path_to_summary in run_summaries_fp_map.items():
        df = pd.read_csv(filepath_or_buffer=path_to_summary, sep="\t")
        df["sample_id"] = sample_id
        all_summaries.append(df)

    return pd.concat(all_summaries, ignore_index=True)


def _get_feature_table(busco_results: pd.DataFrame):
    df = busco_results.reset_index(inplace=False, drop=False)

    new_cols = {
        "mag_id": "MAG", "sample_id": "Sample", "dataset": "Dataset",
        "single": "% single", "duplicated": "% duplicated",
        "fragmented": "% fragmented", "missing": "% missing",
        "complete": "% complete", "n_markers": "Total markers",
        "contigs_n50": "N50 contigs", "percent_gaps": "Percent gaps",
        "scaffolds": "Contigs", "length": "Length (bp)"
    }

    if len(busco_results["sample_id"].unique()) < 2:
        del new_cols["sample_id"]

    df = df[list(new_cols.keys())].rename(columns=new_cols, inplace=False)
    return df.to_json(orient='split')


def _parse_df_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds several columns required for generation of downloadable
    BUSCO plots.

    Args:
        df (pd.DataFrame): Unformatted DataFrame

    Returns:
        df (pd.DataFrame): Formatted DataFrame
    """
    df = df.reset_index(drop=False, inplace=False)
    df = df.rename(columns={"id": "mag_id"}, inplace=False)

    # fix data types
    df["percent_gaps"] = df["percent_gaps"].str.split(
        '%', expand=True
    )[0].map(float)
    for col in MARKER_COLS:
        df[col] = df[col].map(float)
    df["n_markers"] = df["n_markers"].map(int)

    return df


def _rename_columns(df):
    cols = {
        "Input_file": "input_file", "Dataset": "dataset",
        "Complete": "complete", "Single": "single",
        "Duplicated": "duplicated", "Fragmented": "fragmented",
        "Missing": "missing", "n_markers": "n_markers",
        "Scaffold N50": "scaffold_n50", "Contigs N50": "contigs_n50",
        "Percent gaps": "percent_gaps", "Number of scaffolds": "scaffolds",
        "sample_id": "sample_id"
    }

    cols_reshuffled = [
        "mag_id", "sample_id", "input_file", "dataset", "complete",
        "single", "duplicated", "fragmented", "missing", "n_markers",
        "scaffold_n50", "contigs_n50", "percent_gaps", "scaffolds",
    ]

    df = df.rename(columns=cols, inplace=False)
    df["mag_id"] = df["input_file"].str.split(".", expand=True)[0]

    return df[cols_reshuffled]


def _cleanup_bootstrap(output_dir):
    # Remove unwanted files
    # until Bootstrap 3 is replaced with v5, remove the v3 scripts as
    # the HTML files are adjusted to work with v5
    os.remove(
        os.path.join(
            output_dir, "q2templateassets", "css", "bootstrap.min.css"
        )
    )
    os.remove(
        os.path.join(
            output_dir, "q2templateassets", "js", "bootstrap.min.js"
        )
    )


def _calculate_summary_stats(df: pd.DataFrame) -> json:
    stats = pd.DataFrame({
        "min": df[MARKER_COLS].min(),
        "median": df[MARKER_COLS].median(),
        "mean": df[MARKER_COLS].mean(),
        "max": df[MARKER_COLS].max(),
        "count": df[MARKER_COLS].count()
    })
    return stats.T.to_json(orient='table')


def _get_mag_lengths(bins: Union[MultiMAGSequencesDirFmt, MAGSequencesDirFmt]):
    lengths = {}
    if isinstance(bins, MultiMAGSequencesDirFmt):
        for sample, mags in bins.sample_dict().items():
            for mag_id, mag_fp in mags.items():
                seq = skbio.io.read(mag_fp, format="fasta")
                lengths[mag_id] = sum([len(s) for s in seq])
        return pd.Series(lengths, name="length")
    else:
        for mag_id, mag_fp in bins.feature_dict().items():
            seq = skbio.io.read(mag_fp, format="fasta")
            lengths[mag_id] = sum([len(s) for s in seq])
        return pd.Series(lengths, name="length")


def _validate_parameters(lineage_dataset, auto_lineage,
                         auto_lineage_euk, auto_lineage_prok):
    if not any([lineage_dataset, auto_lineage, auto_lineage_euk, auto_lineage_prok]):
        raise ValueError(
            "At least one of these parameters must be provided/set to True: "
            "'lineage-dataset', 'auto-lineage', 'auto-lineage-euk', "
            "'auto-lineage-prok'."
        )
    if lineage_dataset and any([auto_lineage, auto_lineage_euk, auto_lineage_prok]):
        raise ValueError(
            "If 'lineage-dataset' is provided, all the parameters 'auto-lineage', "
            "'auto-lineage-euk' and 'auto-lineage-prok' must be set to False."
        )
