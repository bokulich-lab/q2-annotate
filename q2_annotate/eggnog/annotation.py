# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import shutil
import subprocess
import warnings
from pathlib import Path
from typing import Union

import pandas as pd

from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_types.feature_map import MAGtoContigsDirFmt
from q2_types.genome_data import (
    OrthologAnnotationDirFmt,
    Orthologs,
    SeedOrthologDirFmt,
    OrthologFileFmt,
)
from q2_types.reference_db import EggnogRefDirFmt
from q2_types.sample_data import SampleData


def _annotate_seed_orthologs_runner(
    seed_ortholog, eggnog_db, sample_label, output_loc, db_in_memory, num_cpus
):
    # at this point instead of being able to specify the type of target
    # orthologs, we want to annotate _all_.

    cmds = [
        "emapper.py",
        "-m",
        "no_search",
        "--annotate_hits_table",
        str(seed_ortholog),
        "--data_dir",
        str(eggnog_db),
        "-o",
        str(sample_label),
        "--output_dir",
        str(output_loc),
        "--cpu",
        str(num_cpus),
    ]
    if db_in_memory:
        cmds.append("--dbmem")

    subprocess.run(cmds, check=True, cwd=str(output_loc))


def _eggnog_annotate(
    eggnog_hits: SeedOrthologDirFmt,
    db: EggnogRefDirFmt,
    db_in_memory: bool = False,
    num_cpus: int = 1,
) -> OrthologAnnotationDirFmt:

    eggnog_db_fp = db.path

    result = OrthologAnnotationDirFmt()

    # run analysis
    for relpath, obj_path in eggnog_hits.seed_orthologs.iter_views(OrthologFileFmt):
        sample_label = str(relpath).rsplit(r".", 2)[0]

        _annotate_seed_orthologs_runner(
            seed_ortholog=obj_path,
            eggnog_db=eggnog_db_fp,
            sample_label=sample_label,
            output_loc=result,
            db_in_memory=db_in_memory,
            num_cpus=num_cpus,
        )

    return result


def map_eggnog(
    ctx, eggnog_hits, db, db_in_memory=False, num_cpus=1, num_partitions=None
):
    _eggnog_annotate = ctx.get_action("annotate", "_eggnog_annotate")
    collate_annotations = ctx.get_action("types", "collate_ortholog_annotations")

    if eggnog_hits.type <= SampleData[Orthologs]:
        partition_method = ctx.get_action("types", "partition_orthologs")
    else:
        raise NotImplementedError()

    (partitioned_orthologs,) = partition_method(eggnog_hits, num_partitions)

    annotations = []
    for orthologs in partitioned_orthologs.values():
        (annotation,) = _eggnog_annotate(orthologs, db, db_in_memory, num_cpus)
        annotations.append(annotation)

    (collated_annotations,) = collate_annotations(annotations)
    return collated_annotations


def _extract_generic(
    data: pd.DataFrame, column: str, transform_func: callable
) -> pd.Series:
    """
    Extracts a series from the annotation DataFrame and applies a
    transformation function to the specified column (annotation).

    Args:
        data (pd.DataFrame): The input DataFrame.
        column (str): The annotation column in the DataFrame to extract.
        transform_func (callable): The transformation function to apply
            to the extracted column.

    Returns:
        pd.Series: The extracted series with the transformation applied.
    """
    data = data.set_index("seed_ortholog", inplace=False)[column]
    data = data.apply(transform_func if data.notna().any() else lambda x: x)
    data = data.stack().reset_index(level=1, drop=True)
    data = data.value_counts()
    for char in ("", "-"):
        if char in data.index:
            data = data.drop(char, axis=0, inplace=False)
    return data


# this dictionary contains all the supported annotation types
# each value represents a tuple of:
# 1. original annotation column name (as it appears in the annotation table)
# 2. lambda function which will process values of that column and expand them
#   into a new series of values
extraction_methods = {
    "cog": ("COG_category", lambda x: pd.Series(list(x))),
    "kegg_ko": ("KEGG_ko", lambda x: pd.Series([i[3:] for i in x.split(",")])),
    "kegg_pathway": (
        "KEGG_Pathway",
        lambda x: pd.Series([i for i in x.split(",") if i.startswith("map")]),
    ),
    "kegg_module": ("KEGG_Module", lambda x: pd.Series(x.split(","))),
    "kegg_reaction": ("KEGG_Reaction", lambda x: pd.Series(x.split(","))),
    "brite": ("BRITE", lambda x: pd.Series(x.split(","))),
    "caz": ("CAZy", lambda x: pd.Series(x.split(","))),
    "ec": ("EC", lambda x: pd.Series(x.split(","))),
}


def _filter(data: pd.DataFrame, max_evalue: float, min_score: float) -> pd.DataFrame:
    data = data[(data["evalue"] <= max_evalue) & (data["score"] >= min_score)]
    if len(data) == 0:
        raise ValueError(
            "E-value/score filtering resulted in an empty table - "
            "please adjust your thresholds and try again."
        )
    return data


def extract_annotations(
    ortholog_annotations: OrthologAnnotationDirFmt,
    annotation: str,
    max_evalue: float = 1.0,
    min_score: float = 0.0,
) -> pd.DataFrame:
    extract_method = extraction_methods.get(annotation)
    if not extract_method:
        raise NotImplementedError(f"Annotation '{annotation}' not supported.")
    else:
        col, func = extract_method

    annotations = []
    for _id, fp in ortholog_annotations.annotation_dict().items():
        annot_df = pd.read_csv(
            fp, sep="\t", skiprows=4, index_col=0
        )  # skip the first 4 rows as they contain comments
        annot_df = annot_df.iloc[:-3, :]  # remove the last 3 comment rows
        annot_df = _filter(annot_df, max_evalue, min_score)
        annot_df = _extract_generic(annot_df, col, func)
        annot_df.name = _id
        annotations.append(annot_df)

    result = pd.concat(annotations, axis=1).fillna(0).T
    result.index.name = "id"
    return result


def _get_mag_ids_from_feature_data(mags: MAGSequencesDirFmt) -> set:
    """Extract MAG UUIDs from a FeatureData[MAG] artifact."""
    return set(mags.feature_dict().keys())


def _copy_annotation_files(
    source_annotations: OrthologAnnotationDirFmt,
    mag_ids: set,
    result: OrthologAnnotationDirFmt,
):
    """Copy annotation files from source to result for the given MAG IDs.

    """
    annotation_dict = source_annotations.annotation_dict()
    matched = 0
    for mag_id, src_path in annotation_dict.items():
        if mag_id in mag_ids:
            shutil.copy2(src_path, str(result.path / Path(src_path).name))
            matched += 1

    missing = mag_ids - set(annotation_dict.keys())
    if missing:
        warnings.warn(
            f"{len(missing)} MAG(s) in the destination had no matching "
            f"annotation file in the source and will be skipped: "
            f"{', '.join(sorted(missing))}",
            UserWarning,
        )
    if matched == 0:
        raise ValueError(
            "No annotation files matched the destination MAG IDs. "
            "Make sure the source annotations were derived from the same "
            "set of sequences as the destination."
        )


def transfer_annotations(
    ortholog_annotations: OrthologAnnotationDirFmt,
    reference: Union[MAGSequencesDirFmt, MAGtoContigsDirFmt],
) -> OrthologAnnotationDirFmt:
    """Transfer or aggregate eggNOG annotations based on the reference type.

    If *reference* is a ``FeatureData[MAG]`` (``MAGSequencesDirFmt``), copies
    annotation files for the MAGs present in the reference — useful for
    subsetting annotations after dereplication.

    If *reference* is a ``FeatureMap[MAGtoContigs]`` (``MAGtoContigsDirFmt``),
    aggregates contig-level annotations into MAG-level annotations using the
    contig-to-MAG mapping.
    """
    if isinstance(reference, MAGSequencesDirFmt):
        result = OrthologAnnotationDirFmt()
        mag_ids = _get_mag_ids_from_feature_data(reference)
        _copy_annotation_files(ortholog_annotations, mag_ids, result)
        return result
    else:
        return _annotate_mags_from_contigs(ortholog_annotations, reference)


def _annotate_mags_from_contigs(
    ortholog_annotations: OrthologAnnotationDirFmt,
    contig_map: MAGtoContigsDirFmt,
) -> OrthologAnnotationDirFmt:
    """Aggregate contig-level eggNOG annotations -> MAG-level annotations."""
    # contig_map: {mag_uuid: [contig_id, ...]}
    contig_map_dict = contig_map.view(dict)

    # reverse map: contig_id -> mag_uuid
    contig_to_mag = {
        contig_id: mag_uuid
        for mag_uuid, contig_ids in contig_map_dict.items()
        for contig_id in contig_ids
    }

    # Read all annotation files into one DataFrame
    frames = []
    for _id, fp in ortholog_annotations.annotation_dict().items():
        df = pd.read_csv(fp, sep="\t", skiprows=4, index_col=False)
        df = df.iloc[:-3]
        frames.append(df)

    all_annotations = pd.concat(frames, ignore_index=True)

    # Strip ORF suffix (e.g. contig_id_1 -> contig_id), then map to MAG UUID
    query_col = all_annotations.columns[0]
    all_annotations["mag_uuid"] = (
        all_annotations[query_col]
        .str.replace(r"_\d+$", "", regex=True)
        .map(contig_to_mag)
    )

    unmatched = all_annotations["mag_uuid"].isna().sum()
    if unmatched > 0:
        warnings.warn(
            f"{unmatched} annotation row(s) could not be matched to any "
            "MAG in the contig map and were skipped.",
            UserWarning,
        )

    matched = all_annotations.dropna(subset=["mag_uuid"])
    if matched.empty:
        raise ValueError(
            "No annotation rows could be matched to any MAG. "
            "Make sure the contig IDs in the contig map match those used "
            "during annotation (e.g., run eggNOG on the same contigs that "
            "were used for binning)."
        )

    result = OrthologAnnotationDirFmt()
    for mag_uuid, group in matched.groupby("mag_uuid"):
        out_fp = result.path / f"{mag_uuid}.emapper.annotations"
        n_contigs = len(contig_map_dict.get(mag_uuid, []))
        n_rows = len(group)
        with open(out_fp, "w") as fh:
            fh.write(
                f"## Annotations derived using annotate_mags_from_contigs\n"
                f"## MAG {mag_uuid}: {n_contigs} contig(s) in map, "
                f"{n_rows} annotation row(s) transferred\n"
            )
        group.drop(columns=["mag_uuid"]).to_csv(out_fp, sep="\t", index=False, mode="a")

    return result
