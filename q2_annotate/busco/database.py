# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import subprocess
from q2_annotate._utils import colorify, run_command
from q2_annotate.busco.types import BuscoDatabaseDirFmt


def fetch_busco_db(
    virus: bool = False,
    prok: bool = False,
    euk: bool = False
) -> BuscoDatabaseDirFmt:
    busco_db = BuscoDatabaseDirFmt(path=None, mode='w')

    # Parse kwargs
    if all([virus, prok, euk]):
        args = ["all"]
    else:
        variable_and_flag = [
            ('virus', virus),
            ('prokaryota', prok),
            ('eukaryota', euk)
        ]
        args = [name for name, flag in variable_and_flag if flag]

    # Download
    print(colorify("Downloading BUSCO database..."))
    try:
        run_command(cmd=["busco", "--download", *args], cwd=str(busco_db))
    except subprocess.CalledProcessError as e:
        raise Exception(
            f"Error during BUSCO database download: {e.returncode}"
        )

    # Let user know that the process is complete but it still needs
    # some time to copy files over.
    print(colorify(
        "Download completed. \n"
        "Copying files from temporary directory to final location..."
    ))

    return busco_db
