# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


def _process_semibin2_arg(arg_key, arg_val):
    """Creates a list with argument and its value to be consumed by SemiBin2.

    Argument names will be converted to command line parameters by
    appending a '--' prefix and replacing underscores with hyphens,
    e.g.: 'min_len' -> '--min-len'.

    Args:
        arg_key (str): Argument name.
        arg_val: Argument value.

    Returns:
        [converted_arg, arg_value]: List containing a prepared command line
            parameter and, optionally, its value.
    """
    # Map parameter names to SemiBin2 CLI arguments
    param_mapping = {
        "threads": "--threads",
        "min_len": "--min-len",
        "batch_size": "--batch-size", 
        "epochs": "--epochs",
        "random_seed": "--random-seed",
        "sequencing_type": "--data-type",
    }
    
    arg_key_flag = param_mapping.get(arg_key, f"--{arg_key.replace('_', '-')}")

    if isinstance(arg_val, bool) and arg_val:
        return [arg_key_flag]
    else:
        return [arg_key_flag, str(arg_val)]