# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
from q2_annotate.busco.database import fetch_busco_db, delete_symlinks
from unittest.mock import patch
from qiime2.plugin.testing import TestPluginBase


class TestFetchBUSCO(TestPluginBase):
    package = "q2_annotate.busco.tests"

    @patch("subprocess.run")
    def test_fetch_busco_db_virus(self, subp_run):
        busco_db = fetch_busco_db(virus=True, prok=False, euk=False)

        # Check that command was called in the expected way
        cmd = ["busco", "--download", "virus"]
        subp_run.assert_called_once_with(cmd, check=True, cwd=str(busco_db))

    @patch("subprocess.run")
    def test_fetch_busco_db_prok_euk(self, subp_run):
        busco_db = fetch_busco_db(virus=False, prok=True, euk=True)

        # Check that command was called in the expected way
        cmd = ["busco", "--download", "prokaryota", "eukaryota"]
        subp_run.assert_called_once_with(cmd, check=True, cwd=str(busco_db))

    @patch("subprocess.run")
    def test_fetch_busco_db_all(self, subp_run):
        busco_db = fetch_busco_db(virus=True, prok=True, euk=True)

        # Check that command was called in the expected way
        cmd = ["busco", "--download", "all"]
        subp_run.assert_called_once_with(cmd, check=True, cwd=str(busco_db))

    def test_delete_symlinks(self):
        temp_dir = self.temp_dir.name
        sub_dir = os.path.join(temp_dir, "subdir")
        os.makedirs(sub_dir)

        # Create symlinks
        symlink1 = os.path.join(temp_dir, "symlink1")
        symlink2 = os.path.join(sub_dir, "symlink2")

        # Run delete function
        delete_symlinks(temp_dir)

        # Check that symlinks were deleted
        self.assertFalse(os.path.exists(symlink1))
        self.assertFalse(os.path.exists(symlink2))
