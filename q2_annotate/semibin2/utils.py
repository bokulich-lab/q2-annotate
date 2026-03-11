# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


def _format_arg(key):
    key_with_dashes = key.replace("_", "-")
    return "--{}".format(key_with_dashes)


def _process_semibin2_arg(arg_key, arg_val):
    """Creates a list with argument and its value to be consumed by SemiBin2.

    Argument names will be converted to command line parameters by
    appending a '--' prefix and concatenating words separated by a '_',
    e.g.: 'some_parameter_x' -> '--some-parameter-x'.

    Args:
        arg_key (str): Argument name.
        arg_val: Argument value.

    Returns:
        [converted_arg, arg_value]: List containing a prepared command line
            parameter and, optionally, its value.
    """
    if arg_key == "training_type":
        if arg_val == "semi":
            return ["--semi-supervised"]
        else:
            return ["--self-supervised"]

    arg_key_flag = _format_arg(arg_key)

    if isinstance(arg_val, bool) and arg_val:
        return [arg_key_flag]

    return [arg_key_flag, str(arg_val)]
