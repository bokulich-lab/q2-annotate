# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import os
import subprocess
from typing import List

import pandas as pd
import qiime2.util
from qiime2.sdk import Context

from q2_types.feature_data import FeatureData
from q2_types.feature_data_mag import MAG, MAGSequencesDirFmt
from q2_types.genome_data import (
    SeedOrthologDirFmt,
    OrthologFileFmt,
    LociDirectoryFormat,
)
from q2_types.per_sample_sequences import (
    Contigs,
    MAGs,
    ContigSequencesDirFmt,
    MultiMAGSequencesDirFmt,
)
from q2_types.sample_data import SampleData


def _create_symlinks(source_dirs: list, target_dir: str):
    """
    Create symbolic links for files from source_dirs in the target directory.

    Parameters:
    - source_dirs: A list of source directories containing files to be linked.
    - target_dir: The target directory where symbolic links will be created.
    """
    for src in source_dirs:
        for filename in os.listdir(str(src)):
            os.symlink(
                os.path.join(str(src), filename), os.path.join(target_dir, filename)
            )


def _run_eggnog_search_pipeline(
    ctx: Context,
    sequences: qiime2.Artifact,
    db: list,
    num_cpus: int,
    db_in_memory: bool,
    num_partitions: int,
    search_action: str,
):
    """
    Run the eggNOG search pipeline on the given sequences.

    This function determines the appropriate partitioning method based on
    the type of input sequences, partitions the sequences, runs the
    specified search action on each partition, and then collates the results.

    Parameters:
    - ctx: The QIIME 2 context object.
    - sequences: The input sequences to be analyzed.
    - db: A list containing database-related inputs
        (e.g.: pressed_hmm_db, idmap, seed_alignments).
    - num_cpus: Number of CPUs to use for the search.
    - db_in_memory: Boolean indicating whether the database should
        be loaded into memory.
    - num_partitions: Number of partitions to divide the input sequences into.
    - search_action: The specific search action (function name as
        a string) to be executed.

    Returns:
    - collated_hits: The results_all ortholog hits.
    - collated_tables: The results_all feature tables.
    """
    if sequences.type <= FeatureData[MAG]:
        plugin, action_name = "types", "partition_feature_data_mags"
    elif sequences.type <= SampleData[Contigs]:
        plugin, action_name = "assembly", "partition_contigs"
    elif sequences.type <= SampleData[MAGs]:
        plugin, action_name = "types", "partition_sample_data_mags"
    else:
        raise NotImplementedError()

    partition_method = ctx.get_action(plugin, action_name)
    _eggnog_search = ctx.get_action("annotate", search_action)
    collate_hits = ctx.get_action("types", "collate_orthologs")
    _eggnog_feature_table = ctx.get_action("annotate", "_eggnog_feature_table")
    collate_loci = ctx.get_action("types", "collate_loci")
    (partitioned_sequences,) = partition_method(sequences, num_partitions)

    hits, loci = [], []
    for seq in partitioned_sequences.values():
        (hit, _, locus) = _eggnog_search(seq, *db, num_cpus, db_in_memory)
        hits.append(hit)
        loci.append(locus)

    (collated_hits,) = collate_hits(hits)
    (collated_tables,) = _eggnog_feature_table(collated_hits)
    (collated_loci,) = collate_loci(loci)
    return collated_hits, collated_tables, collated_loci


def _search_runner(
    input_path,
    sample_label: str,
    output_loc: str,
    num_cpus: int,
    db_in_memory: bool,
    runner_args: List[str],
):
    """
    Execute the eggNOG-mapper command with specified arguments for
    a single sample.

    Parameters:
    - input_path: Path to the input file containing sequences for the sample.
    - sample_label: Label for the sample, used in naming output files.
    - output_loc: Directory where the output files will be stored.
    - num_cpus: Number of CPUs to use for the search.
    - db_in_memory: Boolean indicating whether the database should be loaded
        into memory.
    - runner_args: Additional arguments to pass to the eggNOG-mapper command.
    """
    cmd = [
        "emapper.py",
        "-i",
        str(input_path),
        "-o",
        sample_label,
        "-m",
        *runner_args,
        "--genepred",
        "prodigal",
        "--itype",
        "metagenome",
        "--output_dir",
        output_loc,
        "--cpu",
        str(num_cpus),
        "--no_annot",
    ]
    if db_in_memory:
        cmd.append("--dbmem")

    try:
        subprocess.run(cmd, check=True, cwd=output_loc)
    except subprocess.CalledProcessError as e:
        raise Exception("Error running eggNOG-mapper. " f"The exception was: {e}")


def _eggnog_search(
    sequences, search_runner, output_loc
) -> (SeedOrthologDirFmt, pd.DataFrame):
    # run analysis
    if isinstance(sequences, ContigSequencesDirFmt):
        for sample_id, contigs_fp in sequences.sample_dict().items():
            search_runner(input_path=contigs_fp, sample_label=sample_id)
    elif isinstance(sequences, MAGSequencesDirFmt):
        for mag_id, mag_fp in sequences.feature_dict().items():
            search_runner(input_path=mag_fp, sample_label=mag_id)
    elif isinstance(sequences, MultiMAGSequencesDirFmt):
        for sample_id, mags in sequences.sample_dict().items():
            for mag_id, mag_fp in mags.items():
                search_runner(input_path=mag_fp, sample_label=mag_id)

    # iterate over the gff files and move them to the correct location
    loci = LociDirectoryFormat()
    gff_fp = [
        os.path.basename(x) for x in glob.glob(f"{output_loc}/*.emapper.genepred.gff")
    ]
    for fn in gff_fp:
        new_fn = fn.replace(".emapper.genepred.gff", ".gff")
        qiime2.util.duplicate(
            os.path.join(output_loc, fn), os.path.join(loci.path, new_fn)
        )

    orthologs = SeedOrthologDirFmt()
    ortholog_fps = [
        os.path.basename(x) for x in glob.glob(f"{output_loc}/*.seed_orthologs")
    ]
    for item in ortholog_fps:
        qiime2.util.duplicate(
            os.path.join(output_loc, item), os.path.join(orthologs.path, item)
        )

    ft = _eggnog_feature_table(orthologs)
    return orthologs, ft, loci


def _eggnog_feature_table(seed_orthologs: SeedOrthologDirFmt) -> pd.DataFrame:
    per_sample_counts = []

    for sample_path, obj in seed_orthologs.seed_orthologs.iter_views(OrthologFileFmt):
        # TODO: put filename to sample name logic on OrthologFileFmt object
        sample_name = str(sample_path).replace(".emapper.seed_orthologs", "")
        sample_df = obj.view(pd.DataFrame)
        sample_feature_counts = sample_df.value_counts("sseqid")
        sample_feature_counts.name = str(sample_name)
        per_sample_counts.append(sample_feature_counts)
    df = pd.DataFrame(per_sample_counts)
    df.fillna(0, inplace=True)
    df.columns = df.columns.astype("str")

    return df
