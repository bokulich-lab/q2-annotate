# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import subprocess
from copy import deepcopy
from typing import Union

import pandas as pd
from q2_types.per_sample_sequences import (
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt
)

from q2_moshpit._utils import run_command, _process_common_input_params
from q2_moshpit.kraken2.utils import _process_kraken2_arg
from q2_types_genomics.kraken2 import (
    Kraken2ReportDirectoryFormat,
    Kraken2OutputDirectoryFormat,
    Kraken2DBDirectoryFormat
)
from q2_types_genomics.per_sample_data import MultiMAGSequencesDirFmt


def _get_seq_paths(df_index, df_row, df_columns):
    if "filename" in df_columns:
        _sample, _bin, fn = df_index[0], df_index[1], [df_row["filename"]]
    elif "reverse" in df_columns:
        _sample, _bin, fn = df_index, df_index, df_row.tolist()
    else:
        _sample, _bin, fn = df_index, df_index, [df_row["forward"]]
    return _sample, _bin, fn


def _construct_output_paths(
    _sample, _bin, kraken2_outputs_dir, kraken2_reports_dir
):
    sample_dir_report = os.path.join(kraken2_reports_dir.path, _sample)
    sample_dir_output = os.path.join(kraken2_outputs_dir.path, _sample)
    for s in [sample_dir_report, sample_dir_output]:
        os.makedirs(s, exist_ok=True)
    report_fp = os.path.join(sample_dir_report, f"{_bin}.report.txt")
    output_fp = os.path.join(sample_dir_output, f"{_bin}.output.txt")
    return output_fp, report_fp


def _classify_kraken(
    manifest, common_args
) -> (Kraken2ReportDirectoryFormat, Kraken2OutputDirectoryFormat):
    base_cmd = ["kraken2", *common_args]
    base_cmd.append("--paired") if "reverse" in manifest.columns else False

    kraken2_reports_dir = Kraken2ReportDirectoryFormat()
    kraken2_outputs_dir = Kraken2OutputDirectoryFormat()

    try:
        for index, row in manifest.iterrows():
            _sample, _bin, fn = _get_seq_paths(index, row, manifest.columns)
            output_fp, report_fp = _construct_output_paths(
                _sample, _bin, kraken2_outputs_dir, kraken2_reports_dir
            )
            cmd = deepcopy(base_cmd)
            cmd.extend(
                ["--report", report_fp, "--output", output_fp,
                 "--use-names", *fn]
            )
            run_command(cmd=cmd, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running Kraken 2, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )

    return kraken2_reports_dir, kraken2_outputs_dir


def classify_kraken(
    seqs: Union[
        SingleLanePerSamplePairedEndFastqDirFmt,
        SingleLanePerSampleSingleEndFastqDirFmt,
        MultiMAGSequencesDirFmt,
    ],
    db: Kraken2DBDirectoryFormat,
    threads: int = 1,
    confidence: float = 0.0,
    minimum_base_quality: int = 0,
    memory_mapping: bool = False,
    minimum_hit_groups: int = 2,
    quick: bool = False,
) -> (Kraken2ReportDirectoryFormat, Kraken2OutputDirectoryFormat):
    kwargs = {k: v for k, v in locals().items() if k not in ["seqs", "db"]}
    common_args = _process_common_input_params(
        processing_func=_process_kraken2_arg, params=kwargs
    )
    common_args.extend(["--db", str(db.path)])
    manifest: pd.DataFrame = seqs.manifest.view(pd.DataFrame)

    return _classify_kraken(manifest, common_args)
