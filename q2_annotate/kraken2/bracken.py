# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import platform
import re
import subprocess
import tempfile

import pandas as pd

from q2_annotate._utils import run_command
from q2_annotate.kraken2.select import kraken2_to_features
from q2_types.kraken2 import (
    Kraken2ReportDirectoryFormat,
    BrackenDBDirectoryFormat,
)


def _run_bracken_one_sample(
    bracken_db: str,
    kraken2_report_fp: str,
    bracken_report_dir: str,
    tmp_dir: str,
    threshold: int,
    read_len: int,
    level: str,
) -> pd.DataFrame:
    sample_id = os.path.basename(kraken2_report_fp).replace(".report.txt", "")
    bracken_output_fp = os.path.join(tmp_dir, f"{sample_id}.bracken.output.txt")
    bracken_report_fp = os.path.join(bracken_report_dir, f"{sample_id}.report.txt")
    cmd = [
        "bracken",
        "-d",
        bracken_db,
        "-i",
        str(kraken2_report_fp),
        "-o",
        bracken_output_fp,
        "-w",
        bracken_report_fp,
        "-t",
        str(threshold),
        "-r",
        str(read_len),
        "-l",
        level,
    ]
    try:
        run_command(cmd=cmd, verbose=True)
    except subprocess.CalledProcessError as e:
        # TODO: what should be the behaviour when no reads could be classified?
        raise Exception(
            "An error was encountered while running Bracken, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )
    bracken_table = pd.read_csv(bracken_output_fp, sep="\t", index_col=0)
    bracken_table["sample_id"] = sample_id
    bracken_table["taxonomy_id"] = bracken_table["taxonomy_id"].astype(str)

    return bracken_table


def _estimate_bracken(
    kraken2_reports: Kraken2ReportDirectoryFormat,
    db: BrackenDBDirectoryFormat,
    threshold: int,
    read_len: int,
    level: str,
) -> (pd.DataFrame, Kraken2ReportDirectoryFormat):
    bracken_tables = []
    bracken_reports = Kraken2ReportDirectoryFormat()

    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            for report_fp in kraken2_reports.path.iterdir():
                bracken_table = _run_bracken_one_sample(
                    bracken_db=str(db),
                    kraken2_report_fp=report_fp,
                    bracken_report_dir=str(bracken_reports),
                    tmp_dir=tmpdir,
                    threshold=threshold,
                    read_len=read_len,
                    level=level,
                )
                bracken_tables.append(bracken_table)
        except subprocess.CalledProcessError as e:
            raise Exception(
                "An error was encountered while running Bracken, "
                f"(return code {e.returncode}), please inspect "
                "stdout and stderr to learn more."
            )

    bracken_table = pd.concat(bracken_tables).reset_index()
    bracken_table = bracken_table.pivot(
        index="sample_id", columns="taxonomy_id", values="new_est_reads"
    )
    bracken_table.fillna(0, inplace=True)

    return bracken_table, bracken_reports


def _assert_read_lens_available(bracken_db: BrackenDBDirectoryFormat, read_len: int):
    pattern = r".+database(\d{2,})mers\.kmer_distrib$"
    lengths = []
    for db in bracken_db.path.iterdir():
        lengths.extend(re.findall(pattern, str(db)))
    lengths = sorted([int(x) for x in lengths])
    if read_len not in lengths:
        raise ValueError(
            f"Provided read length ({read_len}) is not available in the "
            f"Bracken DB. The available values are: "
            f"{', '.join([str(x) for x in lengths])}."
        )


def _add_unclassified(
    table: pd.DataFrame, taxonomy: pd.Series, reports: Kraken2ReportDirectoryFormat
) -> (pd.DataFrame, pd.Series):
    """
    Adds unclassified read counts from Kraken2 reports to
    the given feature table and taxonomy.

    This function reads the Kraken2 report files for each sample and extracts
    the number of unclassified reads. It adds these counts to the table under
    the taxonomy ID '0'. If the '0' taxonomy ID does not exist in the
    taxonomy series, it is added with the description 'd__Unclassified'.
    Samples without unclassified reads will have a count of 0 added to the
    table for taxonomy ID '0'.

    Args:
        table (pd.DataFrame): DataFrame with sample read counts.
        taxonomy (pd.Series): Series mapping taxonomy IDs to names.
        reports (Kraken2ReportDirectoryFormat): Directory containing
            Kraken2 report files.

    Returns:
        pd.DataFrame: Updated table with unclassified read counts.
        pd.Series: Updated taxonomy with '0' as 'd__Unclassified'
            if not already present.
    """
    unclassified_counts = {}
    found_unclassified = False

    for sample in table.index:
        report_fp = os.path.join(reports.path, f"{sample}.report.txt")
        with open(report_fp, "r") as f:
            line = f.readline()
            if "unclassified" in line:
                unclassified_counts[sample] = int(line.split("\t")[1])
                found_unclassified = True
            else:
                unclassified_counts[sample] = 0

    if found_unclassified and "0" not in taxonomy:
        taxonomy.loc["0"] = "d__Unclassified"

    table["0"] = table.index.map(unclassified_counts).fillna(0)

    return table, taxonomy


def estimate_bracken(
    kraken2_reports: Kraken2ReportDirectoryFormat,
    db: BrackenDBDirectoryFormat,
    threshold: int = 0,
    read_len: int = 100,
    level: str = "S",
    include_unclassified: bool = True,
) -> (Kraken2ReportDirectoryFormat, pd.DataFrame, pd.DataFrame):

    # Check for macOS and raise an error if detected
    if platform.system() == "Darwin":
        raise RuntimeError(
            "The estimate_bracken action is currently not supported on macOS "
            "due to the unavailability of the Bracken package. Please try it on Linux."
        )

    _assert_read_lens_available(db, read_len)

    table, reports = _estimate_bracken(
        kraken2_reports=kraken2_reports,
        db=db,
        threshold=threshold,
        read_len=read_len,
        level=level,
    )

    _, taxonomy = kraken2_to_features(reports=reports, coverage_threshold=0.0)

    if include_unclassified:
        # Bracken does not report unclassified reads in its output table,
        # so we need to re-add them
        table, taxonomy = _add_unclassified(table, taxonomy, kraken2_reports)

    return reports, taxonomy, table
