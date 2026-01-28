# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import subprocess
import hashlib
import warnings
from typing import List, Union

import pandas as pd
import skbio
import biom
import numpy as np
from q2_types.feature_data import DNASequencesDirectoryFormat

from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_types.feature_table import FeatureTable, PresenceAbsence, RelativeFrequency
from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt, ContigSequencesDirFmt

EXTERNAL_CMD_WARNING = (
    "Running external command line application(s). "
    "This may print messages to stdout and/or stderr.\n"
    "The command(s) being run are below. These commands "
    "cannot be manually re-run as they will depend on "
    "temporary files that no longer exist."
)


def run_command(cmd, env=None, verbose=True, pipe=False, **kwargs):
    if verbose:
        print(EXTERNAL_CMD_WARNING)
        print("\nCommand:", end=" ")
        print(" ".join(cmd), end="\n\n")

    if pipe:
        result = subprocess.run(
            cmd, env=env, check=True, capture_output=True, text=True
        )
        return result

    if env:
        subprocess.run(cmd, env=env, check=True, **kwargs)
    else:
        subprocess.run(cmd, check=True, **kwargs)


def _construct_param(arg_name):
    """Converts argument name into a command line parameter."""
    return f'--{arg_name.replace("_", "-")}'


def _process_common_input_params(processing_func, params: dict) -> List[str]:
    """Converts provided arguments and their values.

    Conversion is entirely dependent on the passed 'processing_func'
    that processes individual arguments. The output is a list of
    parameters with their values that can be directly passed to the
    respective command.

    Arguments without any value are skipped.
    Any other argument is processed using the 'processing_func' and
    appended to the final list.

    Args:
        processing_func: Function to be used for formatting a single argument.
        params (dict): Dictionary of parameter: value pairs to be processed.

    Returns:
        processed_args (list): List of processed arguments and their values.

    """
    processed_args = []
    for arg_key, arg_val in params.items():
        # This if condition excludes arguments which are falsy
        # (False, None, "", []), except for integers and floats.
        if type(arg_val) == int or type(arg_val) == float or arg_val:  # noqa: E721
            processed_args.extend(processing_func(arg_key, arg_val))

    return processed_args


def colorify(string: str):
    return "%s%s%s" % ("\033[1;32m", string, "\033[0m")


def _calculate_md5_from_file(file_path: str) -> str:
    md5_hash = hashlib.md5()
    with open(file_path, "rb") as f:
        # Read the file in chunks to handle large files
        for chunk in iter(lambda: f.read(4096), b""):
            md5_hash.update(chunk)
    return md5_hash.hexdigest()


def get_feature_lengths(
    features: Union[
        MAGSequencesDirFmt,
        MultiMAGSequencesDirFmt,
        ContigSequencesDirFmt,
        DNASequencesDirectoryFormat,
    ],
) -> pd.DataFrame:
    """Calculate lengths of features in a feature data object."""
    ids, lengths = [], []
    per_sequence = False

    # for FeatureData[Sequence]
    if isinstance(features, DNASequencesDirectoryFormat):
        sample_dict = {"": {"": os.path.join(str(features), "dna-sequences.fasta")}}
        per_sequence = True

        # For SampleData[Contigs]
    elif isinstance(features, ContigSequencesDirFmt):
        sample_dict = {"": features.sample_dict()}
        per_sequence = True

    # For SampleData[MAGs]
    elif isinstance(features, MultiMAGSequencesDirFmt):
        sample_dict = features.sample_dict()

    # For FeatureData[MAG]
    else:
        sample_dict = {"": features.feature_dict()}

    for _, file_dict in sample_dict.items():
        for _id, fp in file_dict.items():
            sequences = skbio.io.read(fp, format="fasta", verify=False)

            if per_sequence:
                for seq in sequences:
                    ids.append(seq.metadata["id"])
                    lengths.append(len(seq))
            else:
                ids.append(_id)
                lengths.append(sum(len(seq) for seq in sequences))

    df = pd.DataFrame({"id": ids, "length": lengths})
    df.set_index("id", inplace=True)
    return df


def _multiply(table1: biom.Table, table2: biom.Table) -> biom.Table:
    """Calculate dot product of two biom tables."""
    # Subset table1 to only include observations present in table2's samples
    table2_sample_ids = set(table2.ids(axis="sample"))
    table1_obs_to_keep = [
        obs_id
        for obs_id in table1.ids(axis="observation")
        if obs_id in table2_sample_ids
    ]

    if not table1_obs_to_keep:
        raise ValueError(
            "No overlapping features found between table1 observations and table2 samples."
        )

    if len(table1_obs_to_keep) < len(table1.ids(axis="observation")):
        warnings.warn(
            f"Removed {len(table1.ids(axis='observation')) - len(table1_obs_to_keep)} features from table1."
        )

    table1 = table1.filter(table1_obs_to_keep, axis="observation")

    # Reorder table2 samples to match table1 observations
    table2 = table2.sort_order(table1.ids(axis="observation"), axis="sample")

    # Perform sparse matrix multiplication
    # table1: samples x observations, table2: samples x observations
    # result: table1.samples x table2.observations
    result_matrix = table1.matrix_data.T.dot(table2.matrix_data.T)

    result_table = biom.Table(
        result_matrix.T,
        observation_ids=table2.ids(axis="observation"),
        sample_ids=table1.ids(axis="sample"),
    )

    return result_table


def _multiply_tables(table1: biom.Table, table2: biom.Table) -> biom.Table:
    """Calculate dot product of two feature tables."""
    result = _multiply(table1, table2)
    return result


def _multiply_tables_relative(table1: biom.Table, table2: biom.Table) -> biom.Table:
    """Calculate dot product of two feature tables and convert to
    a relative frequency table."""
    result = _multiply(table1, table2)
    result.norm(axis="sample", inplace=True)
    return result


def _multiply_tables_pa(table1: biom.Table, table2: biom.Table) -> biom.Table:
    """Calculate dot product of two feature tables and convert to
    a presence-absence table."""
    result = _multiply(table1, table2)
    # Convert to presence-absence (1 if non-zero, 0 otherwise)
    result_data = result.matrix_data.copy()
    result_data.data = np.ones_like(result_data.data)
    result = biom.Table(
        result_data,
        observation_ids=result.ids(axis="observation"),
        sample_ids=result.ids(axis="sample"),
    )
    return result


def multiply_tables(ctx, table1, table2):
    """Calculate dot product of two feature tables."""
    if (
        table1.type <= FeatureTable[PresenceAbsence]
        or table2.type <= FeatureTable[PresenceAbsence]
    ):
        multiply = ctx.get_action("annotate", "_multiply_tables_pa")
    elif (
        table1.type <= FeatureTable[RelativeFrequency]
        or table2.type <= FeatureTable[RelativeFrequency]
    ):
        multiply = ctx.get_action("annotate", "_multiply_tables_relative")
    else:
        multiply = ctx.get_action("annotate", "_multiply_tables")
    (result,) = multiply(table1, table2)
    return result
