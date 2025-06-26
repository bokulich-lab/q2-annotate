# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


def _capitalize_and_join(key):
    split_key = key.split("_")
    caps = [x.capitalize() for x in split_key[1:]]
    return f'--{split_key[0]}{"".join(caps)}'


def _process_metabat2_arg(arg_key, arg_val):
    """Creates a list with argument and its value to be consumed by MetaBAT 2.

    Argument names will be converted to command line parameters by
    appending a '--' prefix and concatenating words separated by a '_',
    e.g.: 'some_parameter_x' -> '--someParameterX'.

    Args:
        arg_key (str): Argument name.
        arg_val: Argument value.

    Returns:
        [converted_arg, arg_value]: List containing a prepared command line
            parameter and, optionally, its value.
    """
    if arg_key == "p_tnf":
        arg_key_flag = "--pTNF"
    elif arg_key == "min_cv":
        arg_key_flag = "--minCV"
    elif arg_key == "min_cv_sum":
        arg_key_flag = "--minCVSum"
    else:
        arg_key_flag = _capitalize_and_join(arg_key)

    if isinstance(arg_val, bool) and arg_val:
        return [arg_key_flag]
    else:
        return [arg_key_flag, str(arg_val)]
