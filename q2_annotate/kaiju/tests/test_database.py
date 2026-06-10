# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import tempfile
import unittest
from unittest.mock import patch, Mock

import pandas as pd
from qiime2.plugin.testing import TestPluginBase

from q2_annotate.kaiju.database import (
    _fetch_and_extract_db,
    _find_latest_db_url,
    fetch_kaiju_db,
    build_kaiju_db,
    _fetch_ncbi_taxonomy_files,
    _create_protein_fasta,
    CHUNK_SIZE,
    ERR_MSG,
    KAIJU_SERVER_URL,
)
from requests.exceptions import ConnectionError, RequestException

from q2_types.kaiju import KaijuDBDirectoryFormat
from q2_types.genome_data import ProteinsDirectoryFormat
from qiime2 import Metadata


class TestDatabaseFunctions(TestPluginBase):
    package = "q2_annotate.kaiju.tests"

    def setUp(self):
        super().setUp()
        self.html_content = b""" # noqa: E501
            <html><body>
            <table>
              <thead>
                <tr>
                  <th>Database</th>
                  <th>Date</th>
                  <th style="text-align: right">Archive size (GB)</th>
                  <th style="text-align: right">RAM needed (GB)</th>
                  <th>HTTPS URL</th>
                </tr>
              </thead>
              <tbody>
                <tr><td><strong>nr</strong></td><td colspan="4">Subset of NCBI BLAST <a href="https://ftp.ncbi.nlm.nih.gov/blast/db/v5/v5/FASTA/">nr</a> database containing Archaea, bacteria and viruses.</td></tr>
                <tr>
                  <td></td>
                  <td>2023-05-10</td>
                  <td style="text-align: right">67</td>
                  <td style="text-align: right">177</td>
                  <td><a href="https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_nr_2023-05-10.tgz">.tgz</a></td>
                </tr>
                <tr><td><strong>nr_euk</strong></td><td colspan="4">Like nr, but additionally including fungi and microbial eukaryotes, see <a href="https://github.com/bioinformatics-centre/kaiju/blob/master/util/kaiju-taxonlistEuk.tsv">taxon list</a></td></tr>
                <tr>
                  <td></td>
                  <td>2023-05-10</td>
                  <td style="text-align: right">82</td>
                  <td style="text-align: right">204</td>
                  <td><a href="https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_nr_euk_2023-05-10.tgz">.tgz</a></td>
                </tr>
                <tr><td><strong>refseq</strong></td><td colspan="4">Protein sequences from genome assemblies of Archaea and bacteria with assembly level "Complete Genome", as well as viral protein sequences from NCBI RefSeq.</td></tr>
                <tr>
                  <td></td>
                  <td>2023-05-10</td>
                  <td style="text-align: right">30</td>
                  <td style="text-align: right">87</td>
                  <td><a href="https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_refseq_2023-05-23.tgz">.tgz</a></td>
                </tr>
              </tbody>
            </table>
            </body></html>
            """

    @patch("requests.get")
    @patch("q2_annotate.kaiju.database.tqdm")
    @patch("tarfile.open")
    @patch("os.remove")
    def test_fetch_and_extract_db(
        self, mock_remove, mock_tarfile_open, mock_progress, mock_requests
    ):
        response = mock_requests.return_value
        response.headers = {"content-length": 1024}
        response.iter_content.return_value = [b"test"] * 1024
        mock_tar = Mock()
        mock_tarfile_open.return_value.__enter__.return_value = mock_tar

        with tempfile.TemporaryDirectory() as tmpdir:
            _fetch_and_extract_db("http://a/b/db.tar.gz", tmpdir)
            db_path = os.path.join(tmpdir, "db.tar.gz")

            mock_progress.assert_called_with(
                desc='Downloading the "db.tar.gz" database',
                total=1024,
                unit="B",
                unit_scale=True,
                unit_divisor=1024,
            )
            response.iter_content.assert_called_with(chunk_size=CHUNK_SIZE)
            mock_tarfile_open.assert_called_with(db_path, "r:gz")
            mock_tar.extractall.assert_called_with(path=tmpdir)
            mock_remove.assert_called_with(db_path)
            mock_requests.assert_called_with("http://a/b/db.tar.gz", stream=True)

    @patch("requests.get", side_effect=ConnectionError("some error"))
    def test_fetch_and_extract_db_exception(self, mock_requests):
        exp_error = ERR_MSG.format("some error")
        with self.assertRaisesRegex(Exception, exp_error):
            with tempfile.TemporaryDirectory() as tmpdir:
                _fetch_and_extract_db("http://a/b/db.tar.gz", tmpdir)

                mock_requests.assert_called_with("http://a/b/db.tar.gz", stream=True)

    def test_find_latest_db_url_ok(self):
        url = _find_latest_db_url(self.html_content, "nr")
        self.assertEqual(
            url,
            "https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/"
            "kaiju_db_nr_2023-05-10.tgz",
        )

    def test_find_latest_db_url_not_found(self):
        with self.assertRaises(ValueError) as context:
            _find_latest_db_url(self.html_content, "non_existing_db")
        self.assertIn(
            "URL for database type 'non_existing_db' not found.", str(context.exception)
        )

    @patch("requests.get")
    @patch("q2_annotate.kaiju.database._fetch_and_extract_db")
    def test_fetch_kaiju_db(self, mock_fetch, mock_requests):
        mock_requests.return_value = Mock(content=self.html_content)
        obs_db = fetch_kaiju_db("nr_euk")
        self.assertIsInstance(obs_db, KaijuDBDirectoryFormat)
        mock_requests.assert_called_with(KAIJU_SERVER_URL)
        mock_fetch.assert_called_with(
            "https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/"
            "kaiju_db_nr_euk_2023-05-10.tgz",
            str(obs_db.path),
        )

    @patch("requests.get", side_effect=RequestException("some error"))
    def test_fetch_kaiju_db_exception(self, mock_requests):
        with self.assertRaisesRegex(Exception, ERR_MSG.format("some error")):
            fetch_kaiju_db("nr_euk")

        mock_requests.assert_called_with(KAIJU_SERVER_URL)

    @patch("tarfile.open")
    @patch("requests.get")
    def test_fetch_ncbi_taxonomy_files(self, mock_requests, mock_tarfile):
        # Mock the response
        mock_response = Mock()
        mock_response.headers = {"content-length": "1024"}
        mock_response.iter_content.return_value = [b"test"] * 1024
        mock_requests.return_value = mock_response
        
        # Mock the tarfile
        mock_tar = Mock()
        mock_member1 = Mock()
        mock_member1.name = "nodes.dmp"
        mock_member2 = Mock()
        mock_member2.name = "names.dmp"
        mock_member3 = Mock()
        mock_member3.name = "other_file.txt"
        mock_tar.getmembers.return_value = [mock_member1, mock_member2, mock_member3]
        mock_tarfile.return_value.__enter__.return_value = mock_tar
        
        with tempfile.TemporaryDirectory() as tmp_dir:
            _fetch_ncbi_taxonomy_files(tmp_dir)
            
            # Verify only the required files were extracted
            mock_tar.extract.assert_any_call(mock_member1, tmp_dir)
            mock_tar.extract.assert_any_call(mock_member2, tmp_dir)
            # Should not extract other_file.txt
            self.assertEqual(mock_tar.extract.call_count, 2)

    def test_create_protein_fasta(self):
        # Create a temporary proteins directory
        with tempfile.TemporaryDirectory() as tmp_dir:
            proteins_dir = ProteinsDirectoryFormat(tmp_dir, mode='w')
            
            # Create some test protein files
            genome1_path = os.path.join(proteins_dir.path, "genome1.fasta")
            genome2_path = os.path.join(proteins_dir.path, "genome2.fasta")
            
            with open(genome1_path, 'w') as f:
                f.write(">protein1\nMKLLILLFFFLILSSPGGLGKNTKTKIIT\n")
                f.write(">protein2\nMLKIFASDFASDFASDFASDFASD\n")
            
            with open(genome2_path, 'w') as f:
                f.write(">protein3\nMQKLIFASDFASDFASDFASDFASD\n")
            
            # Create metadata
            metadata_dict = {
                'genome1': {'taxon_id': 12345},
                'genome2': {'taxon_id': 67890}
            }
            metadata_df = pd.DataFrame.from_dict(metadata_dict, orient='index')
            metadata = Metadata(metadata_df)
            
            # Test the function
            output_path = os.path.join(tmp_dir, "combined.fasta")
            _create_protein_fasta(proteins_dir, metadata, output_path)
            
            # Verify the output file exists and has correct content
            self.assertTrue(os.path.exists(output_path))
            with open(output_path, 'r') as f:
                content = f.read()
                self.assertIn("_taxid_12345", content)
                self.assertIn("_taxid_67890", content)
                self.assertIn("MKLLILLFFFLILSSPGGLGKNTKTKKIIT", content)

    @patch("q2_annotate.kaiju.database.run_command")
    @patch("q2_annotate.kaiju.database._fetch_ncbi_taxonomy_files")
    @patch("q2_annotate.kaiju.database._create_protein_fasta")
    @patch("shutil.copy2")
    def test_build_kaiju_db(self, mock_copy, mock_create_fasta, mock_fetch_tax, mock_run_command):
        # Create mock inputs
        with tempfile.TemporaryDirectory() as tmp_dir:
            proteins_dir = ProteinsDirectoryFormat(tmp_dir, mode='w')
            
            # Create simple metadata
            metadata_dict = {'genome1': {'taxon_id': 12345}}
            metadata_df = pd.DataFrame.from_dict(metadata_dict, orient='index')
            metadata = Metadata(metadata_df)
            
            # Call the function
            result = build_kaiju_db(proteins_dir, metadata)
            
            # Verify it returns a KaijuDBDirectoryFormat
            self.assertIsInstance(result, KaijuDBDirectoryFormat)
            
            # Verify the required functions were called
            mock_fetch_tax.assert_called_once()
            mock_create_fasta.assert_called_once()
            
            # Verify kaiju commands were run
            self.assertEqual(mock_run_command.call_count, 2)
            
            # Verify files were copied
            self.assertEqual(mock_copy.call_count, 3)  # fmi, nodes.dmp, names.dmp


if __name__ == "__main__":
    unittest.main()
