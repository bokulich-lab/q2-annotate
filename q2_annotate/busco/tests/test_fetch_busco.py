# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
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
        busco_db = fetch_busco_db(lineages=["virus"])
        cmd = ["busco", "--download_path", str(busco_db), "--download", "virus"]
        subp_run.assert_called_once_with(cmd, check=True)

    @patch("subprocess.run")
    def test_fetch_busco_db_prok_euk(self, subp_run):
        busco_db = fetch_busco_db(lineages=["prokaryota", "eukaryota"])
        cmd = [
            "busco",
            "--download_path",
            str(busco_db),
            "--download",
            "prokaryota",
            "eukaryota",
        ]
        subp_run.assert_called_once_with(cmd, check=True)

    @patch("subprocess.run")
    def test_fetch_busco_db_all(self, subp_run):
        busco_db = fetch_busco_db(lineages=["all"])
        cmd = ["busco", "--download_path", str(busco_db), "--download", "all"]
        subp_run.assert_called_once_with(cmd, check=True)

    @patch("subprocess.run")
    def test_fetch_busco_db_two_lineages(self, subp_run):
        busco_db = fetch_busco_db(lineages=["lineage1", "lineage2"])
        cmd = [
            "busco",
            "--download_path",
            str(busco_db),
            "--download",
            "lineage1",
            "lineage2",
        ]
        subp_run.assert_called_once_with(cmd, check=True)

    @patch("subprocess.run")
    def test_fetch_busco_db_lineages_with_all(self, subp_run):
        busco_db = fetch_busco_db(lineages=["lineage1", "lineage2", "all"])
        cmd = ["busco", "--download_path", str(busco_db), "--download", "all"]
        subp_run.assert_called_once_with(cmd, check=True)

    @patch("subprocess.run")
    def test_fetch_busco_db_lineages_with_domains(self, subp_run):
        busco_db = fetch_busco_db(
            lineages=["lineage1", "lineage2", "prokaryota", "virus"]
        )
        cmd = [
            "busco",
            "--download_path",
            str(busco_db),
            "--download",
            "prokaryota",
            "virus",
        ]
        subp_run.assert_called_once_with(cmd, check=True)

    @patch("subprocess.run")
    def test_fetch_busco_db_no_lineage(self, subp_run):
        with self.assertRaisesRegex(ValueError, "No lineages provided."):
            fetch_busco_db()

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
