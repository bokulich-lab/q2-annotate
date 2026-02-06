# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import subprocess
from collections import defaultdict

import pandas as pd

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


def _extract_generic(
    data: pd.DataFrame, column: str, func: callable
) -> tuple[dict, pd.Series, pd.DataFrame]:
    """
    Converts annotation data to a feature map and counts of annotation values.

    Processes the annotation DataFrame by extracting and expanding annotation
    values from the specified column, grouping them by contig ID, and counting
    their occurrences. The function applies a transformation function to expand
    annotation values and removes empty or dash values.

    Args:
        data (pd.DataFrame): The input annotation DataFrame.
        column (str): The annotation column in the DataFrame to extract.
        func (callable): The transformation function to apply to expand
            annotation values into lists.

    Returns:
        tuple[dict, pd.Series, pd.DataFrame]: A tuple containing:
            - dict: Feature map with annotation values as keys and lists of
              unique contig IDs as values.
            - pd.Series: Value counts of all annotation values.
            - pd.DataFrame: Counts of each annotation value per contig.
    """
    data = data.copy()
    data["contig_id"] = data.index.astype(str).str.rsplit("_", n=1).str[0]

    def _expand(x):
        s = func(x)
        if isinstance(s, pd.Series):
            return s.tolist()
        # be tolerant if a future extractor returns a scalar or list
        if isinstance(s, (list, tuple, set)):
            return list(s)
        return [s]

    tmp = data[["contig_id", column]].dropna(subset=[column]).copy()
    tmp["annotation_value"] = tmp[column].map(_expand)
    tmp = tmp.explode("annotation_value")

    # drop empties/dashes and any remaining missing values
    tmp = tmp[
        tmp["annotation_value"].notna() & ~tmp["annotation_value"].isin(("", "-"))
    ]

    # count annotations per contig
    contig_annotation_counts = (
        tmp.groupby("contig_id")["annotation_value"]
        .value_counts()
        .unstack(fill_value=0)
    )

    # count all the annotation values
    annotation_counts = tmp["annotation_value"].value_counts()

    # for the feature map, ensure each (annotation_value, contig_id) pair appears once
    deduplicated = tmp.drop_duplicates(subset=["annotation_value", "contig_id"])

    feature_map = (
        deduplicated.groupby("annotation_value")["contig_id"].agg(list).to_dict()
    )

    return feature_map, annotation_counts, contig_annotation_counts


def _merge_maps(maps: list[dict]) -> dict:
    """Merges a list of feature maps into a single dictionary."""
    merged = defaultdict(set)
    for d in maps:
        for k, v in d.items():
            merged[k].update(v)
    return {k: list(v) for k, v in merged.items()}


def extract_annotations(
    ortholog_annotations: OrthologAnnotationDirFmt,
    annotation: str,
    max_evalue: float = 1.0,
    min_score: float = 0.0,
) -> (pd.DataFrame, dict, pd.DataFrame):
    extract_method = extraction_methods.get(annotation)
    if not extract_method:
        raise NotImplementedError(f"Annotation '{annotation}' not supported.")
    else:
        col, func = extract_method

    annotations, feature_maps, contig_counts = [], [], []
    for _id, fp in ortholog_annotations.annotation_dict().items():
        annot_df = pd.read_csv(
            fp, sep="\t", skiprows=4, index_col=0
        )  # skip the first 4 rows as they contain comments
        annot_df = annot_df.iloc[:-3, :]  # remove the last 3 comment rows
        annot_df = _filter(annot_df, max_evalue, min_score)

        # get feature map, annotation counts, and contig annotation counts
        feature_map, annot_df, contig_annot_counts = _extract_generic(
            annot_df, col, func
        )
        annot_df.name = _id

        feature_maps.append(feature_map)
        annotations.append(annot_df)
        contig_counts.append(contig_annot_counts)

    result = pd.concat(annotations, axis=1).fillna(0).T
    result.index.name = "id"

    merged_maps = _merge_maps(feature_maps)
    contig_result = pd.concat(contig_counts, axis=0).fillna(0)

    return result, dict(merged_maps), contig_result
