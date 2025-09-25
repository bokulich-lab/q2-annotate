# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import base64
import glob
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Union

import pandas as pd
from q2_types.per_sample_sequences import (
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
    SequencesWithQuality,
    PairedEndSequencesWithQuality,
    ContigSequencesDirFmt,
    Contigs,
    JoinedSequencesWithQuality,
)

from q2_annotate._utils import run_command
from q2_types.kaiju import KaijuDBDirectoryFormat
from q2_types.sample_data import SampleData

DEFAULT_PREFIXES = ["d__", "p__", "c__", "o__", "f__", "g__", "s__", "ssp__"]
RANKS = [x.replace("__", "") for x in DEFAULT_PREFIXES]


def _get_sample_paths(df_index, df_row, paired):
    """
    Args:
        df_index: The index of the DataFrame.
        df_row: A row from the dataframe.
        paired: A boolean flag indicating whether the data is paired.

    Returns:
        sample_name: The name of the sample.
        fps: A list of file paths.

    """
    if paired:
        sample_name, fps = df_index, df_row.tolist()
    else:
        sample_name, fps = df_index, [df_row["forward"], ""]
    return sample_name, fps


def _rename_taxon(_id: str, id_to_taxon: dict) -> str:
    """
    Args:
        _id (str): The ID of the taxon to be renamed.
        id_to_taxon (dict): A dictionary mapping taxon IDs to taxon names.

    Returns:
        str: The renamed taxon as a string.

    """
    _id = id_to_taxon[_id].rstrip(";").split(";")
    x = zip(_id, DEFAULT_PREFIXES)
    return ";".join([f"{prefix}{taxon}" for taxon, prefix in x])


def _encode_unclassified_ids(table: pd.DataFrame, text: str) -> pd.DataFrame:
    """
    Args:
        table (pd.DataFrame): The input DataFrame with the taxon information.
        text (str): The text to filter the taxon_name column in the table.

    Returns:
        pd.DataFrame: The updated DataFrame with encoded unclassified ids.

    """
    taxon = table.loc[table["taxon_name"].str.startswith(text), "taxon_name"].iloc[0]
    encoded = base64.b64encode(taxon.encode()).decode()
    table.loc[table["taxon_name"].str.startswith(text), "taxon_id"] = encoded[:8]
    return table


def _fix_id_types(table: pd.DataFrame) -> pd.DataFrame:
    """
    Fixes the data types of the "taxon_id" column in a given DataFrame.

    Args:
        table: A pandas DataFrame containing the data.

    Returns:
        A new pandas DataFrame with the "taxon_id" column modified by:
            - Filling any missing values with 0.
            - Converting the "taxon_id" column to integer type.
            - Encoding the values "cannot be assigned" and "unclassified"
                as integers.

    """
    table["taxon_id"].fillna(0, inplace=True)
    table["taxon_id"] = table["taxon_id"].astype(int)
    table = _encode_unclassified_ids(table, "cannot be assigned")
    table = _encode_unclassified_ids(table, "unclassified")
    return table


def _construct_feature_table(table_fp: str) -> (pd.DataFrame, pd.DataFrame):
    """
    Args:
        table_fp (str): The file path of the table.

    Returns:
        pd.DataFrame, pd.DataFrame: A tuple containing two pandas DataFrames.
            The first DataFrame is the feature table, where the rows
            represent samples and the columns represent features (taxon IDs).
            The second DataFrame contains the taxonomy, where the index
            contains taxon IDs and the only column represents taxon names.

    """
    table = pd.read_csv(table_fp, sep="\t")

    # clean up taxon IDs
    table = _fix_id_types(table)

    # extract sample name from the file path
    table["sample"] = table["file"].map(lambda x: Path(x).stem)

    # clean up all the NAs
    table["taxon_name"] = table["taxon_name"].str.replace("NA", "Unspecified")

    # rename taxon IDs to taxon names and clean up
    taxa = table.set_index("taxon_id")["taxon_name"].to_dict()
    table["taxon_name"] = table["taxon_id"].map(lambda x: _rename_taxon(x, taxa))

    # create taxonomy table
    taxonomy = table[["taxon_id", "taxon_name"]].drop_duplicates()
    taxonomy.set_index("taxon_id", inplace=True, drop=True)
    taxonomy.index = taxonomy.index.astype(str)
    taxonomy.index.name = "Feature ID"
    taxonomy.columns = ["Taxon"]

    # convert to sample x feature format
    table = table.groupby(["taxon_id", "sample"], as_index=False)["reads"].sum()
    table = table.pivot(index="sample", columns="taxon_id", values="reads")

    # convert column names to strings
    table.columns = table.columns.astype(str)

    return table, taxonomy


def _process_kaiju_reports(tmpdir, all_args):
    """
    Args:
        tmpdir (str): The temporary directory where the Kaiju report files
            are located.
        all_args (dict): A dictionary containing the original arguments
            passed to the classification action.

    Returns:
        pd.DataFrame: The feature table constructed from all reports.

    """
    table_args = [
        "-r",
        all_args["r"],
        "-t",
        os.path.join(str(all_args["db"].path), "nodes.dmp"),
        "-n",
        os.path.join(str(all_args["db"].path), "names.dmp"),
        "-l",
        "superkingdom,phylum,class,order,family,genus,species",
    ]
    if all_args["exp"]:
        table_args.append("-e")
    if all_args["u"]:
        table_args.append("-u")
    if 0 <= all_args["c"] <= 1:
        table_args.extend(["-m", str(all_args["c"])])
    else:
        table_args.extend(["-c", str(all_args["c"])])

    report_fps = sorted(glob.glob(os.path.join(tmpdir, "*.out")))

    table_fp = os.path.join(tmpdir, "results.tsv")

    cmd = ["kaiju2table", "-v", "-o", table_fp, *table_args, *report_fps]
    try:
        run_command(cmd=cmd, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running Kaiju2Table, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )

    return _construct_feature_table(table_fp)


def _classify_kaiju_helper(seqs, all_args: dict) -> (pd.DataFrame, pd.DataFrame):
    """
    Args:
        seqs: Sequences object, can be of class SingleLanePerSampleSingleEndFastqDirFmt,
              SingleLanePerSamplePairedEndFastqDirFmt, ContigSequencesDirFmt,
              MAGSequencesDirFmt or MultiFASTADirectoryFormat.
        all_args (dict): A dictionary containing arguments for running Kaiju.

    Returns:
        (pd.DataFrame, pd.DataFrame): A tuple containing two DataFrames -
            the feature table and taxonomy generated by Kaiju.

    Raises:
        Exception: If an error is encountered while running Kaiju.

    """
    kaiju_args = [
        "-z",
        str(all_args["z"]),
        "-a",
        str(all_args["a"]),
        "-e",
        str(all_args["e"]),
        "-m",
        str(all_args["m"]),
        "-s",
        str(all_args["s"]),
        "-E",
        str(all_args["evalue"]),
        "-x" if all_args["x"] else "-X",
        "-t",
        os.path.join(str(all_args["db"].path), "nodes.dmp"),
        "-f",
        glob.glob(os.path.join(str(all_args["db"].path), "kaiju_*.fmi"))[0],
    ]

    base_cmd = ["kaiju-multi", "-v", *kaiju_args]
    read_types = (
        SingleLanePerSampleSingleEndFastqDirFmt,
        SingleLanePerSamplePairedEndFastqDirFmt,
    )
    files_fwd, files_rev, output_fps = [], [], []
    paired = False

    with tempfile.TemporaryDirectory() as tmpdir:
        if isinstance(seqs, read_types):
            manifest: pd.DataFrame = seqs.manifest.view(pd.DataFrame)
            paired = "reverse" in manifest.columns

            for index, row in manifest.iterrows():
                sample_name, fps = _get_sample_paths(index, row, paired)
                files_fwd.append(fps[0])
                files_rev.append(fps[1])
                output_fps.append(f"{os.path.join(tmpdir, sample_name)}.out")

        elif isinstance(seqs, ContigSequencesDirFmt):
            for sample_id, fp in seqs.sample_dict().items():
                files_fwd.append(fp)
                output_fps.append(f"{os.path.join(tmpdir, sample_id)}.out")

        base_cmd.extend(["-i", ",".join(files_fwd)])
        if paired:
            base_cmd.extend(["-j", ",".join(files_rev)])
        base_cmd.extend(["-o", ",".join(output_fps)])

        try:
            run_command(cmd=base_cmd, verbose=True)
        except subprocess.CalledProcessError as e:
            raise Exception(
                "An error was encountered while running Kaiju, "
                f"(return code {e.returncode}), please inspect "
                "stdout and stderr to learn more."
            )

        table, taxonomy = _process_kaiju_reports(tmpdir, all_args)

    return table, taxonomy


def _classify_kaiju(
    seqs: Union[
        SingleLanePerSamplePairedEndFastqDirFmt,
        SingleLanePerSampleSingleEndFastqDirFmt,
        ContigSequencesDirFmt,
    ],
    db: KaijuDBDirectoryFormat,
    z: int = 1,
    a: str = "greedy",
    e: int = 3,
    m: int = 11,
    s: int = 65,
    evalue: float = 0.01,
    x: bool = True,
    r: str = "species",
    c: float = 0.0,
    exp: bool = False,
    u: bool = False,
) -> (pd.DataFrame, pd.DataFrame):

    return _classify_kaiju_helper(seqs, dict(locals().items()))


def classify_kaiju(
    ctx,
    seqs,
    db,
    z=1,
    a="greedy",
    e=3,
    m=11,
    s=65,
    evalue=0.01,
    x=True,
    r="species",
    c=0.0,
    exp=False,
    u=False,
    num_partitions=None,
):
    kwargs = {
        k: v
        for k, v in locals().items()
        if k not in ["seqs", "db", "ctx", "num_partitions"]
    }

    _classify_kaiju = ctx.get_action("annotate", "_classify_kaiju")
    collate_feature_tables = ctx.get_action("feature_table", "merge")
    collate_taxonomies = ctx.get_action("feature_table", "merge_taxa")

    if seqs.type <= SampleData[SequencesWithQuality | JoinedSequencesWithQuality]:
        partition_method = ctx.get_action("demux", "partition_samples_single")
    elif seqs.type <= SampleData[PairedEndSequencesWithQuality]:
        partition_method = ctx.get_action("demux", "partition_samples_paired")
    elif seqs.type <= SampleData[Contigs]:
        partition_method = ctx.get_action("types", "partition_contigs")
    else:
        raise NotImplementedError()

    (partitioned_seqs,) = partition_method(seqs, num_partitions)

    tables = []
    taxonomies = []
    for seq in partitioned_seqs.values():
        (table, taxonomy) = _classify_kaiju(seq, db, **kwargs)
        tables.append(table)
        taxonomies.append(taxonomy)

    (combined_table,) = collate_feature_tables(tables)
    (collated_taxonomy,) = collate_taxonomies(taxonomies)

    return combined_table, collated_taxonomy
