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
from typing import List, Union

import pandas as pd
import skbio
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


def _multiply(table1: pd.DataFrame, table2: pd.DataFrame) -> pd.DataFrame:
    """Calculate dot product of two tables."""
    if not table1.columns.equals(table2.index):
        raise ValueError(
            "Columns of the first table do not match the index of the second."
        )
    return table1.dot(table2)


def _multiply_tables(table1: pd.DataFrame, table2: pd.DataFrame) -> pd.DataFrame:
    """Calculate dot product of two feature tables."""
    result = _multiply(table1, table2)
    return result


def _multiply_tables_relative(
    table1: pd.DataFrame, table2: pd.DataFrame
) -> pd.DataFrame:
    """Calculate dot product of two feature tables and convert to
    a relative frequency table."""
    result = _multiply(table1, table2)
    sum_per_sample = result.sum(axis=1)
    result = result.div(sum_per_sample, axis=0)
    return result


def _multiply_tables_pa(table1: pd.DataFrame, table2: pd.DataFrame) -> pd.DataFrame:
    """Calculate dot product of two feature tables and convert to
    a presence-absence table."""
    result = _multiply(table1, table2)
    result = result.applymap(lambda x: 1 if x != 0 else 0)
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
