# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import subprocess
from typing import List

from q2_annotate._utils import colorify, run_command
from q2_annotate.busco.types import BuscoDatabaseDirFmt


def _process_lineages(lineages):
    """Process and validate BUSCO lineage selections.

    This function validates and normalizes lineage selections for BUSCO
    database downloads. It handles special cases like 'all' and domain-level
    lineages (prokaryota, eukaryota, virus).

    Parameters
    ----------
    lineages : list
        List of lineage names to process. Can include 'all' or domain names.

    Returns
    -------
    list
        Processed list of lineage names.

    Raises
    ------
    ValueError
        If no lineages are provided (empty list or None).
    """
    if not lineages:
        raise ValueError("No lineages provided.")

    if "all" in lineages:
        lineages = ["all"]
        return lineages

    domain_lineages = ["prokaryota", "eukaryota", "virus"]
    domain_matches = [lin for lin in lineages if lin in domain_lineages]
    if domain_matches:
        print(
            colorify(
                f"Domain lineages were provided ({', '.join(domain_matches)}) - "
                "other lineages, if any, will be ignored."
            )
        )
        lineages = domain_matches

    return lineages


def delete_symlinks(root_dir):
    """Delete all symbolic links in a directory tree.

    This function walks through a directory tree and removes any symbolic
    links found in the file system. It does not follow symbolic links while
    traversing the directory structure.

    Parameters
    ----------
    root_dir : str
        Path to the root directory to start searching for symbolic links.
    """
    for dirpath, dirnames, filenames in os.walk(root_dir, followlinks=False):
        for _file in filenames:
            if os.path.islink(os.path.join(dirpath, _file)):
                os.unlink(os.path.join(dirpath, _file))


def fetch_busco_db(lineages: List[str] = None) -> BuscoDatabaseDirFmt:
    busco_db = BuscoDatabaseDirFmt(path=None, mode="w")

    lineages = _process_lineages(lineages)

    print(colorify(f"Fetching lineages: {', '.join(lineages)}."))
    cmd = ["busco", "--download_path", str(busco_db), "--download", *lineages]
    try:
        run_command(cmd)
    except subprocess.CalledProcessError as e:
        raise Exception(f"Error during BUSCO database download: {e.returncode}")

    # There is a symlink in the BUSCO database that needs to be removed
    delete_symlinks(str(busco_db))

    print(
        colorify(
            "Download completed. \n"
            "Copying files from temporary directory to the final location..."
        )
    )

    return busco_db
