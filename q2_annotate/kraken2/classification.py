# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import itertools
import os
import subprocess
from copy import deepcopy
from typing import Union, Optional

import pandas as pd
from q2_types.per_sample_sequences import (
    SequencesWithQuality,
    PairedEndSequencesWithQuality,
    JoinedSequencesWithQuality,
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
    ContigSequencesDirFmt, Contigs,
    MultiFASTADirectoryFormat, MAGs
)
from q2_types.sample_data import SampleData
from q2_types.feature_data import FeatureData
from q2_annotate._utils import run_command, _process_common_input_params
from q2_annotate.kraken2.utils import _process_kraken2_arg
from q2_types.feature_data_mag import MAGSequencesDirFmt, MAG
from q2_types.kraken2 import (
    Kraken2ReportDirectoryFormat,
    Kraken2OutputDirectoryFormat,
    Kraken2DBDirectoryFormat,
    Kraken2ReportFormat,
    Kraken2OutputFormat,
)


def _get_seq_paths(df_index, df_row, df_columns):
    if "reverse" in df_columns:
        _sample, fn = df_index, df_row.tolist()
    else:
        _sample, fn = df_index, [df_row["forward"]]
    return _sample, fn


def _construct_output_paths(
        _sample, kraken2_outputs_dir, kraken2_reports_dir
):
    report_fp = os.path.join(
        kraken2_reports_dir.path, f"{_sample}.report.txt"
    )
    output_fp = os.path.join(
        kraken2_outputs_dir.path, f"{_sample}.output.txt"
    )
    return output_fp, report_fp


def classify_kraken2(
    ctx,
    seqs,
    kraken2_db,
    threads=1,
    confidence=0.0,
    minimum_base_quality=0,
    memory_mapping=False,
    minimum_hit_groups=2,
    quick=False,
    report_minimizer_data=False,
    num_partitions=None
):
    '''
    '''
    kwargs = {
        k: v for k, v in locals().items()
        if k not in ["seqs", "kraken2_db", "ctx", "num_partitions"]
    }

    # classify sequences
    reports = []
    outputs = []
    for seqs_artifact in seqs:
        single_art_reports, single_art_outputs = _classify_single_artifact(
            ctx, seqs, kraken2_db, num_partitions, kwargs
        )
        reports.append(single_art_reports)
        outputs.append(single_art_outputs)

    # match and merge collated reports, outputs

    # return matched reports, outputs


def _classify_single_artifact(ctx, seqs, kraken2_db, num_partitions, kwargs):
    '''
    Runs the kraken2 software on the contents of a single artifact.

    Parameters
    ----------
    ctx : qiime2.sdk.Context
        The pipeline context object.
    seqs : qiime2.Artifact
        An artifact of type SampleData[SequencesWithQualtiy],
        SampleData[PairedEndSequencesWithQuality],
        SampleData[JoinedSequencesWithQuality], SampleData[Contigs],
        or SampleData[MAGs].
    kraken2_db : Kraken2DBDirectoryFormat
        The kraken2 database.
    num_partitions : int | None
        The number of partitions to create for parallel execution.
    kwargs : dict
        The remaining keyword arguments from the `classify_kraken2` pipeline.

    Returns
    -------
    tuple[Kraken2ReportDirectoryFormat, Kraken2OutputDirectoryFormat]
        The kraken2 report and output files.
    '''
    _classify_kraken2 = ctx.get_action("annotate", "_classify_kraken2")
    collate_kraken2_reports = ctx.get_action(
        "annotate", "collate_kraken2_reports"
    )
    collate_kraken2_outputs = ctx.get_action(
        "annotate", "collate_kraken2_outputs"
    )

    if seqs.type <= FeatureData[MAG]:
        # FeatureData[MAG] is not parallelized
        return _classify_kraken2(seqs, kraken2_db, **kwargs)
    else:
        partition_action = _get_partition_action(ctx, seqs)
        (partitioned_seqs,) = partition_action(seqs, num_partitions)

        all_reports = []
        all_outputs = []
        for seq in partitioned_seqs.values():
            reports, outputs = _classify_kraken2(seq, kraken2_db, **kwargs)
            all_reports.append(reports)
            all_outputs.append(outputs)

        collated_reports = collate_kraken2_reports(all_reports)
        collated_outputs = collate_kraken2_outputs(all_outputs)

        return collated_reports, collated_outputs


def _merge_kraken2_results(
    reports: list[Kraken2ReportDirectoryFormat],
    outputs: list[Kraken2OutputDirectoryFormat]
) -> tuple[Kraken2ReportDirectoryFormat, Kraken2OutputDirectoryFormat]:
    '''
    Merges kraken2 reports and outputs into a single report and format per
    unique sample id.

    Parameters
    ----------
    reports : list[Kraken2ReportDirectoryFormat]
        The kraken2 reports.
    outputs : list[Kraken2OutputDirectoryFromat]
        The kraken2 outputs.

    Returns
    -------
    tuple[Kraken2ReportDirectoryFormat, Kraken2OutputDirectoryFormat]
        The merged reports and formats.
    '''
    merged_reports = Kraken2ReportDirectoryFormat()
    merged_outputs = Kraken2OutputDirectoryFormat()

    mags = False
    if isinstance(next(iter(reports[0].file_dict().values(), dict))):
        mags = True

    report_mapping, output_mapping = _condense_formats(reports, outputs, mags)

    def _merge_formats(sample_id, merged_formats, mapping, merger, Format):
        for sample_id in report_mapping:
            if mags:
                for filename, formats in mapping[sample_id].items():
                    merged_format = merger(formats)
                    merged_reports.write_data(
                        merged_format,
                        Format,
                        sample_id=sample_id,
                        mag_id=filename
                    )
            else:
                formats = format[sample_id]
                merged_format = merger(formats)
                merged_formats.write_data(
                    merged_format, Format, sample_id=sample_id
                )

    for sample_id in report_mapping:
        _merge_formats(
            sample_id,
            merged_reports,
            report_mapping,
            _merge_reports,
            Kraken2ReportFormat
        )
        _merge_formats(
            sample_id,
            merged_outputs,
            output_mapping,
            _merge_outputs,
            Kraken2OutputFormat
        )

    return merged_reports, merged_outputs


def _condense_formats(
    reports: list[Kraken2ReportDirectoryFormat],
    outputs: list[Kraken2OutputDirectoryFormat],
    mags: bool
) -> tuple[dict, dict]:
    '''
    Condenses multiple report and output directory formats into a single
    mapping each. The structure is sample_id -> list[format] for reads/contigs
    and sample_id -> {filename -> list[format]} for mags.

    Parameters
    ----------
    reports : list[Kraken2ReportDirectoryFormat]
        The kraken2 reports.
    outputs : list[Kraken2OutputDirectoryFormat]
        The kraken2 outputs.
    mags : bool
        Whether the directory formats represent MAG results or not.

    Returns
    -------
    tuple[dict, dict]
        A tuple of mappings as described above. The first contains the kraken2
        reports and the second the kraken2 outputs.
    '''
    chained_reports = itertools.chain(
        *[report.file_dict() for report in reports]
    )
    chained_outputs = itertools.chain(
        *[output.file_dict() for output in outputs]
    )

    def _update_mapping(sample_id, chain, mapping, Format):
        if sample_id not in mapping:
            if mags:
                mapping[sample_id] = {}
                for filename, filepath in chain[sample_id].items():
                    format = Format(filepath, mode='r')
                    mapping[sample_id][filename] = [format]
            else:
                format = Format(chain[sample_id], mode='r')
                mapping[sample_id] = [format]
        else:
            if mags:
                for filename, filepath in chain[sample_id].items():
                    format = Format(filepath, mode='r')
                    mapping[sample_id][filename].append(format)
            else:
                format = Format(chain[sample_id], mode='r')
                mapping[sample_id].append(format)

    report_mapping = {}
    output_mapping = {}

    for sample_id in chained_reports:
        _update_mapping(
            sample_id, chained_reports, report_mapping, Kraken2ReportFormat
        )
        _update_mapping(
            sample_id, chained_outputs, output_mapping, Kraken2OutputFormat
        )

    return report_mapping, output_mapping


def _merge_reports(reports: list[Kraken2ReportFormat]) -> Kraken2ReportFormat:
    # make dataframe for each
    trees = []
    for report in reports:
        df = report.view(pd.DataFrame)
        trees.append(_kraken_to_ncbi_tree(df))

    tree = _combine_ncbi_trees(trees)

    merged_report = Kraken2ReportFormat()

    # dump tree to report



def _merge_outputs(outputs: list[Kraken2OutputFormat]) -> Kraken2OutputFormat:
    pass


def _get_partition_action(ctx, seqs):
    '''
    Returns the proper partition action for the given type of `seqs`.

    Parameters
    ----------
    ctx : qiime2.sdk.Context
        The pipeline context object.
    seqs : qiime2.Artifact
        An artifact of type SampleData[SequencesWithQualtiy],
        SampleData[PairedEndSequencesWithQuality],
        SampleData[JoinedSequencesWithQuality], SampleData[Contigs],
        or SampleData[MAGs].

    Returns
    -------
    qiime2.sdk.Action
        The partition action.
    '''
    if seqs.type <= SampleData[
        SequencesWithQuality | JoinedSequencesWithQuality
    ]:
        return ctx.get_action("demux", "partition_samples_single")
    elif seqs.type <= SampleData[PairedEndSequencesWithQuality]:
        return ctx.get_action("demux", "partition_samples_paired")
    elif seqs.type <= SampleData[Contigs]:
        return ctx.get_action("assembly", "partition_contigs")
    elif seqs.type <= SampleData[MAGs]:
        return ctx.get_action("types", "partition_sample_data_mags")
    else:
        raise NotImplementedError()


def _classify_kraken2(
        seqs: Union[
            SingleLanePerSamplePairedEndFastqDirFmt,
            SingleLanePerSampleSingleEndFastqDirFmt,
            ContigSequencesDirFmt,
            MAGSequencesDirFmt,
            MultiFASTADirectoryFormat
        ],
        kraken2_db: Kraken2DBDirectoryFormat,
        threads: int = 1,
        confidence: float = 0.0,
        minimum_base_quality: int = 0,
        memory_mapping: bool = False,
        minimum_hit_groups: int = 2,
        quick: bool = False,
        report_minimizer_data: bool = False
) -> (
        Kraken2ReportDirectoryFormat,
        Kraken2OutputDirectoryFormat,
):
    kwargs = {k: v for k, v in locals().items()
              if k not in ["seqs", "kraken2_db", "ctx"]}

    common_args = _process_common_input_params(
        processing_func=_process_kraken2_arg, params=kwargs
    )
    common_args.extend(["--db", str(kraken2_db.path)])
    return classify_kraken2_helper(seqs, common_args)


def classify_kraken2_helper(
        seqs, common_args
) -> (Kraken2ReportDirectoryFormat, Kraken2OutputDirectoryFormat):
    base_cmd = ["kraken2", *common_args]

    read_types = (
        SingleLanePerSampleSingleEndFastqDirFmt,
        SingleLanePerSamplePairedEndFastqDirFmt
    )

    if isinstance(seqs, read_types):
        manifest: Optional[pd.DataFrame] = seqs.manifest.view(pd.DataFrame)
        if manifest is not None and "reverse" in manifest.columns:
            base_cmd.append("--paired")

        iterate_over = manifest.iterrows()

        def get_paths_for_reads(index, row):
            return _get_seq_paths(index, row, list(manifest.columns))

        path_function = get_paths_for_reads

    elif isinstance(seqs, ContigSequencesDirFmt):
        iterate_over = seqs.sample_dict().items()

    elif isinstance(seqs, MAGSequencesDirFmt):
        iterate_over = seqs.feature_dict().items()

    elif isinstance(seqs, MultiFASTADirectoryFormat):
        iterate_over = (
            (sample_id, mag_id, mag_fp)
            for sample_id, mags in seqs.sample_dict().items()
            for mag_id, mag_fp in mags.items()
        )

    kraken2_reports_dir = Kraken2ReportDirectoryFormat()
    kraken2_outputs_dir = Kraken2OutputDirectoryFormat()

    try:
        for args in iterate_over:
            if isinstance(seqs, read_types):
                _sample, fps = path_function(*args)
            elif isinstance(seqs, MultiFASTADirectoryFormat):
                sample_id, mag_id, fps = args
                for p in (kraken2_reports_dir.path, kraken2_outputs_dir.path):
                    os.makedirs(os.path.join(p, sample_id), exist_ok=True)
                _sample = f"{sample_id}/{mag_id}"
                fps = [fps]
            else:
                _sample, fps = args
                fps = [fps]

            output_fp, report_fp = _construct_output_paths(
                _sample, kraken2_outputs_dir, kraken2_reports_dir
            )
            cmd = deepcopy(base_cmd)
            cmd.extend(
                ["--report", report_fp, "--output", output_fp, *fps]
            )
            run_command(cmd=cmd, verbose=True)

    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running Kraken 2, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )

    return kraken2_reports_dir, kraken2_outputs_dir
