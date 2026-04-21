# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import gzip
import os
from itertools import zip_longest
from pathlib import Path

import pandas as pd

from q2_annotate.kraken2.select import _get_indentation
from q2_types.kraken2 import (
    Kraken2OutputDirectoryFormat,
    Kraken2ReportDirectoryFormat,
    Kraken2ReportFormat,
)
from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt,
)


def _normalize_taxon_id(taxon_id) -> str:
    taxon_id = str(taxon_id).strip()
    if taxon_id.endswith(".0"):
        taxon_id = taxon_id[:-2]
    return taxon_id


def _normalize_read_id(read_id: str) -> str:
    read_id = str(read_id).strip()
    if read_id.startswith("@"):
        read_id = read_id[1:]
    read_id = read_id.split()[0]

    if read_id.endswith("/1") or read_id.endswith("/2"):
        read_id = read_id[:-2]

    return read_id


def _is_gzip_file(filepath: Path) -> bool:
    with open(filepath, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def _open_fastq_for_read(filepath: Path):
    if _is_gzip_file(filepath):
        return gzip.open(filepath, mode="rt")
    return open(filepath, mode="r")


def _open_fastq_for_write(filepath: Path, gzip_output: bool):
    if gzip_output:
        return gzip.open(filepath, mode="wt")
    return open(filepath, mode="w")


def _assert_distinct_input_output_paths(input_fp: Path, output_fp: Path):
    if input_fp.resolve() == output_fp.resolve():
        raise ValueError(
            "Input and output FASTQ paths are identical. "
            "Refusing in-place overwrite to prevent data loss."
        )

    # `duplicate` may create hardlinks. Writing to a hardlinked destination
    # would also mutate/truncate the input file.
    if output_fp.exists() and os.path.samefile(str(input_fp), str(output_fp)):
        raise ValueError(
            "Input and output FASTQ paths refer to the same inode (hardlink). "
            "Refusing in-place overwrite to prevent data loss."
        )


def _iter_fastq_records(fh):
    while True:
        header = fh.readline()
        if header == "":
            return

        sequence = fh.readline()
        separator = fh.readline()
        quality = fh.readline()
        if sequence == "" or separator == "" or quality == "":
            raise ValueError("Malformed FASTQ file: unexpected end-of-file.")

        yield header, sequence, separator, quality


def _extract_matching_read_ids_from_output(
    output_fp: Path, taxon_ids: set[str]
) -> set[str]:
    if not taxon_ids:
        return set()

    matched_read_ids = set()
    with open(output_fp, "r") as fh:
        for line in fh:
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 3:
                continue

            taxon_id = _normalize_taxon_id(columns[2])
            if taxon_id in taxon_ids:
                matched_read_ids.add(_normalize_read_id(columns[1]))

    return matched_read_ids


def _collect_matching_taxon_ids(
    report: Path,
    taxonomy: str,
    include_descendants: bool = True,
    contains: bool = False,
) -> set[str]:
    report_df = Kraken2ReportFormat(report, "r").view(pd.DataFrame)
    names = report_df["name"].astype(str)
    names_clean = names.str.strip()
    taxon_ids = report_df["taxon_id"].map(_normalize_taxon_id)

    if taxonomy.isdigit():
        match_mask = taxon_ids == _normalize_taxon_id(taxonomy)
    else:
        normalized_query = taxonomy.casefold()
        if contains:
            match_mask = names_clean.str.casefold().str.contains(
                normalized_query, regex=False
            )
        else:
            match_mask = names_clean.str.casefold() == normalized_query

    matched_positions = list(report_df.index[match_mask])
    if not matched_positions:
        return set()

    matched_taxa = set()
    indentations = names.map(_get_indentation)

    for position in matched_positions:
        matched_taxa.add(taxon_ids.iloc[position])

        if not include_descendants:
            continue

        parent_indent = indentations.iloc[position]
        for child_position in range(position + 1, len(report_df)):
            if indentations.iloc[child_position] <= parent_indent:
                break
            matched_taxa.add(taxon_ids.iloc[child_position])

    return matched_taxa


def _filter_single_end_fastq(
    input_fp: Path,
    output_fp: Path,
    matched_read_ids: set[str],
    exclude: bool = False,
):
    _assert_distinct_input_output_paths(input_fp, output_fp)

    gzip_input = _is_gzip_file(input_fp)
    with _open_fastq_for_read(input_fp) as in_fh, _open_fastq_for_write(
        output_fp, gzip_output=gzip_input
    ) as out_fh:
        for record in _iter_fastq_records(in_fh):
            read_id = _normalize_read_id(record[0])
            should_keep = read_id in matched_read_ids
            if exclude:
                should_keep = not should_keep

            if should_keep:
                out_fh.writelines(record)


def _filter_paired_end_fastq(
    forward_input_fp: Path,
    reverse_input_fp: Path,
    forward_output_fp: Path,
    reverse_output_fp: Path,
    matched_read_ids: set[str],
    exclude: bool = False,
):
    _assert_distinct_input_output_paths(forward_input_fp, forward_output_fp)
    _assert_distinct_input_output_paths(reverse_input_fp, reverse_output_fp)

    gzip_forward = _is_gzip_file(forward_input_fp)
    gzip_reverse = _is_gzip_file(reverse_input_fp)

    with (
        _open_fastq_for_read(forward_input_fp) as fwd_in,
        _open_fastq_for_read(reverse_input_fp) as rev_in,
        _open_fastq_for_write(forward_output_fp, gzip_output=gzip_forward) as fwd_out,
        _open_fastq_for_write(reverse_output_fp, gzip_output=gzip_reverse) as rev_out,
    ):
        for fwd_record, rev_record in zip_longest(
            _iter_fastq_records(fwd_in), _iter_fastq_records(rev_in)
        ):
            if fwd_record is None or rev_record is None:
                raise ValueError(
                    "Forward and reverse FASTQ files have different record counts."
                )

            fwd_id = _normalize_read_id(fwd_record[0])
            rev_id = _normalize_read_id(rev_record[0])
            if fwd_id != rev_id:
                raise ValueError(
                    "Forward and reverse FASTQ files are out of sync. "
                    f"Found '{fwd_id}' and '{rev_id}'."
                )

            should_keep = fwd_id in matched_read_ids
            if exclude:
                should_keep = not should_keep

            if should_keep:
                fwd_out.writelines(fwd_record)
                rev_out.writelines(rev_record)


def _validate_read_sample_ids(
    read_sample_ids: set[str],
    report_sample_ids: set[str],
    output_sample_ids: set[str],
) -> None:
    if report_sample_ids != output_sample_ids:
        missing_in_reports = sorted(output_sample_ids - report_sample_ids)
        missing_in_outputs = sorted(report_sample_ids - output_sample_ids)
        raise ValueError(
            "Sample IDs in Kraken2 reports and outputs do not match. "
            f"Missing in reports: {missing_in_reports}. "
            f"Missing in outputs: {missing_in_outputs}."
        )

    missing_in_reads = sorted(report_sample_ids - read_sample_ids)
    if missing_in_reads:
        raise ValueError(
            "Some Kraken2 classification sample IDs are missing from reads. "
            f"Missing in reads: {missing_in_reads}. "
        )

    missing_in_classifications = sorted(read_sample_ids - report_sample_ids)
    if missing_in_classifications:
        raise ValueError(
            "Some read sample IDs are missing from Kraken2 classifications. "
            f"Missing in Kraken2 classifications: {missing_in_classifications}. "
        )


def filter_kraken2_reads(
    reads: CasavaOneEightSingleLanePerSampleDirFmt,
    reports: Kraken2ReportDirectoryFormat,
    outputs: Kraken2OutputDirectoryFormat,
    taxonomy: str,
    include_descendants: bool = True,
    contains: bool = False,
    exclude: bool = False,
) -> CasavaOneEightSingleLanePerSampleDirFmt:
    """
    Filter reads based on Kraken2 classifications that match a taxonomy query.

    Parameters
    ----------
    reads : CasavaOneEightSingleLanePerSampleDirFmt
        Reads used as input for Kraken2 classification.
    reports : Kraken2ReportDirectoryFormat
        Kraken2 reports generated for `reads`.
    outputs : Kraken2OutputDirectoryFormat
        Kraken2 output files generated for `reads`.
    taxonomy : str
        Taxonomy query. This can be a Kraken2 taxon name (e.g., "Bacteria")
        or a taxon ID (e.g., "2").
    include_descendants : bool, optional
        If True, include all descendant taxa under each matching taxon.
    contains : bool, optional
        If True, perform case-insensitive substring matching on taxon names
        instead of exact matching.
    exclude : bool, optional
        If False (default), retain only reads matching the taxonomy query.
        If True, discard matching reads and retain everything else.
    """
    taxonomy = taxonomy.strip()
    if not taxonomy:
        raise ValueError("`taxonomy` cannot be an empty string.")

    sample_fastqs = reads.manifest.to_dict(orient="index")

    report_map = reports.file_dict()
    output_map = outputs.file_dict()

    read_sample_ids = set(sample_fastqs.keys())
    report_sample_ids = set(report_map.keys())
    output_sample_ids = set(output_map.keys())
    _validate_read_sample_ids(read_sample_ids, report_sample_ids, output_sample_ids)

    filtered_reads = CasavaOneEightSingleLanePerSampleDirFmt()

    taxonomy_found = False
    paired = all(fastqs.get("reverse") for fastqs in sample_fastqs.values())

    for sample_id, sample_fps in sample_fastqs.items():
        forward_input = Path(sample_fps["forward"])
        forward_output = filtered_reads.path / forward_input.name

        matched_taxon_ids = _collect_matching_taxon_ids(
            report_map[sample_id],
            taxonomy=taxonomy,
            include_descendants=include_descendants,
            contains=contains,
        )
        if matched_taxon_ids:
            taxonomy_found = True

        matched_read_ids = _extract_matching_read_ids_from_output(
            output_map[sample_id], matched_taxon_ids
        )

        if paired:
            reverse_input = Path(sample_fps["reverse"])
            reverse_output = filtered_reads.path / reverse_input.name
            _filter_paired_end_fastq(
                forward_input_fp=forward_input,
                reverse_input_fp=reverse_input,
                forward_output_fp=forward_output,
                reverse_output_fp=reverse_output,
                matched_read_ids=matched_read_ids,
                exclude=exclude,
            )
        else:
            _filter_single_end_fastq(
                input_fp=forward_input,
                output_fp=forward_output,
                matched_read_ids=matched_read_ids,
                exclude=exclude,
            )

    if not taxonomy_found:
        raise ValueError(
            f"Taxonomy query '{taxonomy}' was not found in any Kraken2 report."
        )

    return filtered_reads
