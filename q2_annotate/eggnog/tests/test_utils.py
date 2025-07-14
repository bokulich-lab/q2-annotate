# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
from unittest.mock import patch, call, MagicMock
from qiime2.plugin.testing import TestPluginBase
from qiime2.core.exceptions import ValidationError
from ..utils import (
    _download_and_build_hmm_db,
    _download_fastas_into_hmmer_db,
    _merge_hmms_and_write_idmap,
    COMMON_URL,
    _validate_eggnog_hmmer_taxon_id,
)
from q2_annotate.eggnog.types import EggnogHmmerIdmapDirectoryFmt
from q2_types.genome_data import ProteinsDirectoryFormat
from q2_types.profile_hmms import (
    PressedProfileHmmsDirectoryFmt,
    ProteinMultipleProfileHmmDirectoryFmt,
)


class TestEggnogUtils(TestPluginBase):
    package = "q2_annotate.eggnog.tests"

    @patch("urllib.request.urlopen")
    def test_validate_eggnog_hmmer_taxon_id_valid(self, mock_urlopen):
        # Mock HTML content
        mock_html = """
        <html>
            <body>
                <a href="1/">1/</a>
                <a href="2/">2/</a>
                <a href="3/">3/</a>
            </body>
        </html>
        """
        mock_response = MagicMock()
        mock_response.read.return_value = mock_html.encode("utf-8")
        mock_urlopen.return_value.__enter__.return_value = mock_response

        _validate_eggnog_hmmer_taxon_id(2)

    @patch("urllib.request.urlopen")
    def test_validate_eggnog_hmmer_taxon_id_invalid(self, mock_urlopen):
        # Mock HTML content
        mock_html = """
        <html>
            <body>
                <a href="1/">1/</a>
                <a href="2/">2/</a>
                <a href="3/">3/</a>
            </body>
        </html>
        """
        mock_response = MagicMock()
        mock_response.read.return_value = mock_html.encode("utf-8")
        mock_urlopen.return_value.__enter__.return_value = mock_response

        with self.assertRaisesRegex(ValidationError, "4 is not a valid taxon ID."):
            _validate_eggnog_hmmer_taxon_id(4)

    @patch("subprocess.run")
    @patch("tempfile.TemporaryDirectory")
    def test_download_and_build_hmm_db(self, tmpdir, mock_run):
        tmp = self.get_data_path("hmmer/hmms")
        taxon_id = 1
        tmpdir.return_value.__enter__.return_value = tmp
        idmap, hmm_db, pressed_hmm_db = _download_and_build_hmm_db(taxon_id)

        self.assertIsInstance(idmap, EggnogHmmerIdmapDirectoryFmt)
        idmap.validate()

        self.assertIsInstance(hmm_db, ProteinMultipleProfileHmmDirectoryFmt)
        hmm_db.validate()

        self.assertIsInstance(pressed_hmm_db, PressedProfileHmmsDirectoryFmt)
        # No validation because files where not created because of patch

        mock_run.assert_has_calls(
            [
                call(
                    [
                        "wget",
                        "-O",
                        f"{tmp}/{taxon_id}_hmms.tar.gz",
                        f"{COMMON_URL}/{taxon_id}/{taxon_id}_hmms.tar.gz",
                    ],
                    check=True,
                ),
                call(["tar", "zxf", f"{taxon_id}_hmms.tar.gz"], check=True, cwd=tmp),
                call(["hmmpress", f"{str(pressed_hmm_db)}/{taxon_id}.hmm"], check=True),
            ]
        )

    @patch("glob.glob")
    @patch("subprocess.run")
    @patch("tempfile.TemporaryDirectory")
    def test_download_fastas_into_hmmer_db(self, tmpdir, mock_run, mock_glob):
        tmp = "tmp"
        taxon_id = 1
        directory_path = self.get_data_path("hmmer/fastas/1")
        tmpdir.return_value.__enter__.return_value = tmp
        mock_glob.return_value = [
            os.path.join(directory_path, f) for f in os.listdir(directory_path)
        ]

        fastas = _download_fastas_into_hmmer_db(1)
        self.assertIsInstance(fastas, ProteinsDirectoryFormat)
        fastas.validate()

        mock_run.assert_has_calls(
            [
                call(
                    [
                        "wget",
                        "-O",
                        f"{tmp}/{taxon_id}_raw_algs.tar",
                        f"{COMMON_URL}/{taxon_id}/{taxon_id}_raw_algs.tar",
                    ],
                    check=True,
                ),
                call(["tar", "xf", f"{taxon_id}_raw_algs.tar"], check=True, cwd=tmp),
            ]
        )

    def test_merge_hmms_and_write_idmap(self):
        hmms = self.get_data_path("hmmer/hmms")
        taxon_id = 1

        merged_hmms_obj = ProteinMultipleProfileHmmDirectoryFmt()
        idmap_obj = EggnogHmmerIdmapDirectoryFmt()

        hmms_merged_p = f"{str(merged_hmms_obj)}/{taxon_id}.hmm"
        idmap_p = f"{str(idmap_obj)}/{taxon_id}.hmm.idmap"

        _merge_hmms_and_write_idmap(hmms_merged_p, idmap_p, taxon_id, hmms)

        merged_hmms_obj.validate()
        idmap_obj.validate()
