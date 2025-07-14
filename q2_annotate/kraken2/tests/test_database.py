# ----------------------------------------------------------------------------
# Copyright (c) 2023-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import io
import os
import shutil
import tarfile
import tempfile
import unittest
from copy import deepcopy
from requests.exceptions import ConnectionError
from subprocess import CalledProcessError
from tempfile import TemporaryDirectory
from unittest.mock import patch, ANY, call, Mock, MagicMock
from parameterized import parameterized

import pandas as pd

from q2_types.feature_data import DNAFASTAFormat
from qiime2 import Artifact
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import annotate

from q2_annotate.kraken2.database import (
    _fetch_taxonomy,
    _fetch_libraries,
    _add_seqs_to_library,
    _build_kraken2_database,
    _move_db_files,
    _build_bracken_database,
    _find_latest_db,
    _fetch_db_collection,
    S3_COLLECTIONS_URL,
    _build_dbs_from_seqs,
    _fetch_prebuilt_dbs,
    inspect_kraken2_db,
    _move_files_one_level_up,
)
from q2_types.kraken2 import (
    Kraken2DBDirectoryFormat,
    BrackenDBDirectoryFormat,
    Kraken2DBReportFormat,
    Kraken2DBReportDirectoryFormat,
)


class MockTempDir(tempfile.TemporaryDirectory):
    pass


class TestKraken2Database(TestPluginBase):
    package = "q2_annotate.kraken2.tests"

    def setUp(self):
        super().setUp()
        self.kraken2_db_dir = "fake/db/dir"
        self.kwargs = {
            "threads": 2,
            "fast_build": True,
            "kmer_len": 31,
            "use_ftp": False,
            "max_db_size": 1000,
            "load_factor": 0.5,
        }
        self.s3_response = b"""
            <ListBucketResult>
                <Contents>
                    <Key>kraken/k2_viral_20201202.tar.gz</Key>
                    <LastModified>2020-12-09T01:38:22.000Z</LastModified>
                </Contents>
                <Contents>
                    <Key>kraken/k2_viral_20230314.tar.gz</Key>
                    <LastModified>2023-03-22T01:29:11.000Z</LastModified>
                </Contents>
            </ListBucketResult>
        """

        self.s3_response_16S = b"""
            <ListBucketResult>
                <Contents>
                    <Key>kraken/16S_Greengenes13.5_20200326.tgz</Key>
                    <LastModified>2020-03-26T01:38:22.000Z</LastModified>
                </Contents>
                <Contents>
                    <Key>kraken/16S_Greengenes13.5_20200326.tgz</Key>
                    <LastModified>2024-02-12T01:29:11.000Z</LastModified>
                </Contents>
            </ListBucketResult>
        """

        self.temp_dir = tempfile.mkdtemp()
        self.temp_tar = os.path.join(self.temp_dir, "temp.tar.gz")

        with tarfile.open(self.temp_tar, "w:gz") as tar:
            data = io.BytesIO(b"sample data")
            tarinfo = tarfile.TarInfo(name="sample.txt")
            tarinfo.size = len(data.getbuffer())
            tar.addfile(tarinfo, data)

        with open(self.temp_tar, "rb") as f:
            self.tar_chunks = [chunk for chunk in iter(lambda: f.read(8192), b"")]

    def create_expected_calls(self, origin_dir, moving_dir, seq) -> list:
        expected_calls = []
        expected_files = [
            "database100mers.kmer_distrib",
            "database150mers.kmer_distrib",
            "database200mers.kmer_distrib",
            "database250mers.kmer_distrib",
            "database50mers.kmer_distrib",
            "database75mers.kmer_distrib",
            "database100mers.kmer_distrib",
            "hash.k2d",
            "opts.k2d",
            "taxo.k2d",
            "README.md",
        ]
        if seq == "silva132" or seq == "rdp":
            expected_files.append("seqid2taxid.map")

        for filename in expected_files:
            expected_calls.append(call(os.path.join(origin_dir, filename), moving_dir))
        return expected_calls

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    @patch("q2_annotate.kraken2.database.run_command")
    def test_fetch_taxonomy(self, p1):
        _fetch_taxonomy(self.kraken2_db_dir, threads=3, use_ftp=True)

        exp_cmd = [
            "kraken2-build",
            "--download-taxonomy",
            "--threads",
            "3",
            "--db",
            self.kraken2_db_dir,
            "--use-ftp",
        ]
        p1.assert_called_once_with(cmd=exp_cmd, verbose=True)

    @patch(
        "q2_annotate.kraken2.database.run_command",
        side_effect=CalledProcessError(123, "cmd"),
    )
    def test_fetch_taxonomy_exception(self, p1):
        with self.assertRaisesRegex(
            Exception,
            "An error was encountered .* downloading taxonomy, "
            r"\(return code 123\), please inspect .*",
        ):
            _fetch_taxonomy(self.kraken2_db_dir, threads=3, use_ftp=True)

    @patch("q2_annotate.kraken2.database.run_command")
    def test_fetch_libraries_skip(self, p1):
        all_kwargs = deepcopy(self.kwargs)
        all_kwargs["library_exists"] = "skip"
        libraries = ["plasmid", "human"]

        with TemporaryDirectory() as tmp_dir:
            for lib in libraries:
                os.makedirs(os.path.join(tmp_dir, "library", lib))
            open(
                os.path.join(tmp_dir, "library", libraries[0], "library.fna"), "w"
            ).close()

            _fetch_libraries(tmp_dir, libraries=libraries, all_kwargs=all_kwargs)

        exp_cmd = [
            "kraken2-build",
            "--download-library",
            libraries[1],
            "--threads",
            "2",
            "--db",
            tmp_dir,
        ]
        p1.assert_called_once_with(cmd=exp_cmd, verbose=True)

    @patch("q2_annotate.kraken2.database.run_command")
    def test_fetch_libraries_refetch(self, p1):
        all_kwargs = deepcopy(self.kwargs)
        all_kwargs["library_exists"] = "refetch"
        libraries = ["plasmid", "human"]

        _fetch_libraries(
            self.kraken2_db_dir, libraries=libraries, all_kwargs=all_kwargs
        )

        base_cmd = ["kraken2-build", "--download-library"]
        exp_common_args = ["--threads", "2", "--db", self.kraken2_db_dir]
        exp_cmds = [
            [*base_cmd, libraries[0], *exp_common_args],
            [*base_cmd, libraries[1], *exp_common_args],
        ]
        p1.assert_has_calls(
            [call(cmd=exp_cmds[0], verbose=True), call(cmd=exp_cmds[1], verbose=True)]
        )

    @patch(
        "q2_annotate.kraken2.database.run_command",
        side_effect=CalledProcessError(123, "cmd"),
    )
    def test_fetch_libraries_exception(self, p1):
        with self.assertRaisesRegex(
            Exception,
            "An error was encountered .* downloading the 'human' "
            r"library, \(return code 123\), please inspect .*",
        ):
            _fetch_libraries(
                self.kraken2_db_dir, libraries=["human"], all_kwargs=self.kwargs
            )

    @patch("q2_annotate.kraken2.database.run_command")
    def test_add_seqs_to_library(self, p1):
        seqs = DNAFASTAFormat(
            self.get_data_path(
                "mags/sample1/3b72d1a7-ddb0-4dc7-ac36-080ceda04aaa.fasta"
            ),
            "r",
        )

        _add_seqs_to_library(self.kraken2_db_dir, seqs=seqs, no_masking=True)

        exp_cmd = [
            "kraken2-build",
            "--add-to-library",
            str(seqs.path),
            "--db",
            self.kraken2_db_dir,
            "--no-mask",
        ]
        p1.assert_called_once_with(cmd=exp_cmd, verbose=True)

    @patch(
        "q2_annotate.kraken2.database.run_command",
        side_effect=CalledProcessError(123, "cmd"),
    )
    def test_add_seqs_to_library_exception(self, p1):
        seqs = DNAFASTAFormat(
            self.get_data_path(
                "mags/sample1/3b72d1a7-ddb0-4dc7-ac36-080ceda04aaa.fasta"
            ),
            "r",
        )

        with self.assertRaisesRegex(
            Exception,
            "An error was encountered .* adding sequences to the "
            r"library, \(return code 123\), please inspect .*",
        ):
            _add_seqs_to_library(self.kraken2_db_dir, seqs=seqs, no_masking=True)

    @patch("q2_annotate.kraken2.database.run_command")
    def test_build_kraken2_database(self, p1):
        _build_kraken2_database(self.kraken2_db_dir, all_kwargs=self.kwargs)

        exp_cmd = [
            "kraken2-build",
            "--build",
            "--db",
            self.kraken2_db_dir,
            "--threads",
            "2",
            "--fast-build",
            "--kmer-len",
            "31",
            "--load-factor",
            "0.5",
            "--max-db-size",
            "1000",
        ]
        p1.assert_called_once_with(cmd=exp_cmd, verbose=True)

    @patch("q2_annotate.kraken2.database.run_command")
    def test_build_kraken2_database_no_max_db(self, p1):
        all_kwargs = deepcopy(self.kwargs)
        all_kwargs["max_db_size"] = 0

        _build_kraken2_database(self.kraken2_db_dir, all_kwargs=all_kwargs)

        exp_cmd = [
            "kraken2-build",
            "--build",
            "--db",
            self.kraken2_db_dir,
            "--threads",
            "2",
            "--fast-build",
            "--kmer-len",
            "31",
            "--load-factor",
            "0.5",
        ]
        p1.assert_called_once_with(cmd=exp_cmd, verbose=True)

    @patch(
        "q2_annotate.kraken2.database.run_command",
        side_effect=CalledProcessError(123, "cmd"),
    )
    def test_build_kraken2_database_exception(self, p1):
        with self.assertRaisesRegex(
            Exception,
            "An error was encountered .* building the database, "
            r"\(return code 123\), please inspect .*",
        ):
            _build_kraken2_database(self.kraken2_db_dir, all_kwargs=self.kwargs)

    @patch("q2_annotate.kraken2.database.run_command")
    def test_build_bracken_database(self, p1):
        _build_bracken_database(
            kraken2_db_dir=self.kraken2_db_dir, threads=2, kmer_len=31, read_len=150
        )

        exp_cmd = [
            "bracken-build",
            "-d",
            self.kraken2_db_dir,
            "-t",
            "2",
            "-k",
            "31",
            "-l",
            "150",
        ]
        p1.assert_called_once_with(cmd=exp_cmd, verbose=True)

    @patch(
        "q2_annotate.kraken2.database.run_command",
        side_effect=CalledProcessError(123, "cmd"),
    )
    def test_build_bracken_database_exception(self, p1):
        with self.assertRaisesRegex(
            Exception,
            "An error was encountered while building the Bracken "
            r"database, \(return code 123\), please inspect .+",
        ):
            _build_bracken_database(
                kraken2_db_dir=self.kraken2_db_dir, threads=2, kmer_len=31, read_len=150
            )

    def test_find_latest_db(self):
        response = Mock(content=self.s3_response)

        obs_db = _find_latest_db("viral", response)
        exp_db = "kraken/k2_viral_20230314.tar.gz"
        self.assertEqual(obs_db, exp_db)

    def test_find_latest_16S_db(self):
        response = Mock(content=self.s3_response_16S)

        obs_db = _find_latest_db("greengenes", response)
        exp_db = "kraken/16S_Greengenes13.5_20200326.tgz"
        self.assertEqual(obs_db, exp_db)

    def test_find_latest_db_empty(self):
        response = Mock(content=b"""<ListBucketResult></ListBucketResult>""")

        with self.assertRaisesRegex(ValueError, r"No databases were found.+"):
            _find_latest_db("viral", response)

    @patch("requests.get")
    @patch("tarfile.open")
    @patch(
        "q2_annotate.kraken2.database._find_latest_db",
        return_value="kraken/k2_viral.tar.gz",
    )
    @patch("q2_annotate.kraken2.database.tqdm")
    def test_fetch_db_collection_success(
        self, mock_tqdm, mock_find, mock_tarfile_open, mock_requests_get
    ):
        mock_requests_get.side_effect = [
            MagicMock(status_code=200),
            MagicMock(
                status_code=200,
                iter_content=lambda chunk_size: self.tar_chunks,
                headers={},
            ),
        ]
        mock_tarfile_open.return_value.__enter__.return_value = MagicMock()

        _fetch_db_collection("viral", "/tmp")

        mock_requests_get.has_calls(
            [
                call(S3_COLLECTIONS_URL),
                call(f"{S3_COLLECTIONS_URL}/kraken/k2_viral.tar.gz", stream=True),
            ]
        )
        mock_tarfile_open.assert_called_once_with("/tmp/k2_viral.tar.gz", "r:gz")
        mock_find.assert_called_once_with("viral", ANY)
        mock_tqdm.assert_not_called()

    @parameterized.expand(
        [
            (
                "kraken/tmp/16S_Greengenes13.5.tgz",
                "16S_Greengenes13.5.tgz",
                "16S_Greengenes13.5",
                "/tmp/16S_Greengenes13.5/",
                "greengenes",
            ),
            (
                "kraken/tmp/16S_RDP11.5_20200326.tgz",
                "16S_RDP11.5_20200326.tgz",
                "16S_RDP_k2db",
                "/tmp/16S_RDP_k2db/",
                "rdp",
            ),
            (
                "kraken/tmp/16S_Silva132_20200326.tgz",
                "16S_Silva132_20200326.tgz",
                "16S_SILVA132_k2db",
                "/tmp/16S_SILVA132_k2db/",
                "silva132",
            ),
            (
                "kraken/tmp/16S_Silva138_20200326.tgz",
                "16S_Silva138_20200326.tgz",
                "16S_SILVA138_k2db",
                "/tmp/16S_SILVA138_k2db/",
                "silva138",
            ),
        ]
    )
    @patch("requests.get")
    @patch("shutil.move")
    @patch("os.listdir")
    @patch("tarfile.open")
    @patch("q2_annotate.kraken2.database.tqdm")
    def test_fetch_db_collection_16S_success(
        self,
        latest_db,
        zipped_folder,
        unzipped_folder,
        folder_name,
        collection,
        mock_tqdm,
        mock_tarfile_open,
        mock_os_listdir,
        mock_shutil_move,
        mock_requests_get,
    ):
        with patch(
            "q2_annotate.kraken2.database._find_latest_db", return_value=latest_db
        ):
            mock_requests_get.side_effect = [
                MagicMock(status_code=200),
                MagicMock(
                    status_code=200,
                    iter_content=lambda chunk_size: self.tar_chunks,
                    headers={},
                ),
            ]

            mock_os_listdir.return_value = [zipped_folder, unzipped_folder]
            expected_calls = self.create_expected_calls(
                folder_name, "/tmp", seq=collection
            )

            _fetch_db_collection(collection, "/tmp")

            mock_shutil_move.has_calls(expected_calls)

            mock_requests_get.has_calls(
                [
                    call(S3_COLLECTIONS_URL),
                    call(f"{S3_COLLECTIONS_URL}/kraken/{zipped_folder}", stream=True),
                ]
            )
            mock_tqdm.assert_not_called()

    @patch("requests.get")
    @patch("tarfile.open")
    @patch(
        "q2_annotate.kraken2.database._find_latest_db",
        return_value="kraken/k2_viral.tar.gz",
    )
    @patch("q2_annotate.kraken2.database.tqdm")
    def test_fetch_db_collection_success_with_tqdm(
        self, mock_tqdm, mock_find, mock_tarfile_open, mock_requests_get
    ):
        mock_requests_get.side_effect = [
            MagicMock(status_code=200),
            MagicMock(
                status_code=200,
                iter_content=lambda chunk_size: self.tar_chunks,
                headers={"content-length": "1000"},
            ),
        ]
        mock_tarfile_open.return_value.__enter__.return_value = MagicMock()

        _fetch_db_collection("viral", "/tmp")

        mock_requests_get.has_calls(
            [
                call(S3_COLLECTIONS_URL),
                call(f"{S3_COLLECTIONS_URL}/kraken/k2_viral.tar.gz", stream=True),
            ]
        )
        mock_tarfile_open.assert_called_once_with("/tmp/k2_viral.tar.gz", "r:gz")
        mock_find.assert_called_once_with("viral", ANY)
        mock_tqdm.assert_called_once_with(
            desc='Downloading the "kraken/k2_viral.tar.gz" database',
            total=1000,
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
        )

    @parameterized.expand(
        [
            (
                "kraken/tmp/16S_Greengenes13.5.tgz",
                "16S_Greengenes13.5.tgz",
                "16S_Greengenes13.5",
                "/tmp/16S_Greengenes13.5/",
                "greengenes",
            ),
            (
                "kraken/tmp/16S_RDP11.5_20200326.tgz",
                "16S_RDP11.5_20200326.tgz",
                "16S_RDP_k2db",
                "/tmp/16S_RDP_k2db/",
                "rdp",
            ),
            (
                "kraken/tmp/16S_Silva132_20200326.tgz",
                "16S_Silva132_20200326.tgz",
                "16S_SILVA132_k2db",
                "/tmp/16S_SILVA132_k2db/",
                "silva132",
            ),
            (
                "kraken/tmp/16S_Silva138_20200326.tgz",
                "16S_Silva138_20200326.tgz",
                "16S_SILVA138_k2db",
                "/tmp/16S_SILVA138_k2db/",
                "silva138",
            ),
        ]
    )
    @patch("requests.get")
    @patch("shutil.move")
    @patch("os.listdir")
    @patch("tarfile.open")
    @patch("q2_annotate.kraken2.database.tqdm")
    def test_fetch_db_collection_16S_tqdm_success(
        self,
        latest_db,
        zipped_folder,
        unzipped_folder,
        folder_name,
        collection,
        mock_tqdm,
        mock_tarfile_open,
        mock_os_listdir,
        mock_shutil_move,
        mock_requests_get,
    ):
        with patch(
            "q2_annotate.kraken2.database._find_latest_db", return_value=latest_db
        ):
            mock_requests_get.side_effect = [
                MagicMock(status_code=200),
                MagicMock(
                    status_code=200,
                    iter_content=lambda chunk_size: self.tar_chunks,
                    headers={"content-length": "1000"},
                ),
            ]

            mock_os_listdir.return_value = [zipped_folder, unzipped_folder]
            expected_calls = self.create_expected_calls(
                folder_name, "/tmp", seq=collection
            )

            _fetch_db_collection(collection, "/tmp")

            mock_shutil_move.has_calls(expected_calls)

            mock_requests_get.has_calls(
                [
                    call(S3_COLLECTIONS_URL),
                    call(f"{S3_COLLECTIONS_URL}/kraken/{zipped_folder}", stream=True),
                ]
            )
            mock_tqdm.assert_called_once_with(
                desc=f'Downloading the "kraken/tmp/{zipped_folder}" database',
                total=1000,
                unit="B",
                unit_scale=True,
                unit_divisor=1024,
            )

    @patch("requests.get")
    def test_fetch_db_collection_connection_error(self, mock_get):
        mock_get.side_effect = ConnectionError("Some error.")
        with self.assertRaisesRegex(ValueError, r".+The error was\: Some error\."):
            _fetch_db_collection("my_collection", "/tmp")

    @patch("requests.get")
    def test_fetch_db_collection_status_non200(self, mock_get):
        mock_get.return_value = Mock(status_code=404)
        with self.assertRaisesRegex(ValueError, r".+Status code was\: 404"):
            _fetch_db_collection("my_collection", "/tmp")

    def test_file_move_one_level_up(self):
        with TemporaryDirectory() as tmp_dir:
            fake_folder = os.path.join(tmp_dir, "test_move")
            fake_file = os.path.join(fake_folder, "folder_1.tgz")
            fake_subfolder = os.path.join(fake_folder, "folder_1")
            fake_subfiles = ["file_1.txt", "file_2.txt"]

            os.makedirs(fake_folder)
            os.makedirs(fake_subfolder)
            open(os.path.join(fake_folder, fake_file), "w").close()

            for f in fake_subfiles:
                open(os.path.join(fake_subfolder, f), "w").close()

            # in the beginning the directory has only the zipped
            # and unzipped folder
            assert len(os.listdir(fake_folder)) == 2

            _move_files_one_level_up(fake_folder)

            # after the moving there must be 4 files/folders,
            # the contents of the unzipped folder, the unzipped folder
            # and the zipped folder
            assert len(os.listdir(fake_folder)) == 4

    def test_move_db_files(self):
        with TemporaryDirectory() as tmp_dir:
            fake_src = os.path.join(tmp_dir, "fake_src")
            fake_dest = os.path.join(tmp_dir, "fake_dest")
            fake_files = ["fake_db.k2d", "fake_file.k2d", "other.file"]

            os.makedirs(fake_src)
            os.makedirs(fake_dest)

            for f in fake_files:
                open(os.path.join(tmp_dir, "fake_src", f), "w").close()

            _move_db_files(fake_src, fake_dest)

            for f in fake_files[:2]:
                self.assertTrue(os.path.exists(os.path.join(fake_dest, f)))

    @patch("q2_annotate.kraken2.database._fetch_taxonomy")
    @patch("q2_annotate.kraken2.database._add_seqs_to_library")
    @patch("q2_annotate.kraken2.database._build_kraken2_database")
    @patch("q2_annotate.kraken2.database._build_bracken_database")
    @patch("q2_annotate.kraken2.database._move_db_files")
    def test_build_dbs_from_seqs(
        self, mock_move, mock_bracken, mock_kraken, mock_add_seqs, mock_fetch_tax
    ):
        bracken_db, kraken2_db = MagicMock(), MagicMock()
        seqs, tmp_dir = ["seq1", "seq2"], "/tmp"
        common_args = {
            "threads": 1,
            "use_ftp": False,
            "no_masking": False,
            "read_len": [100, 150],
            "kmer_len": 35,
        }

        _build_dbs_from_seqs(bracken_db, kraken2_db, seqs, tmp_dir, common_args)

        mock_fetch_tax.assert_called_once_with(db_dir=tmp_dir, threads=1, use_ftp=False)
        mock_add_seqs.assert_has_calls(
            [
                call(db_dir=tmp_dir, seqs="seq1", no_masking=False),
                call(db_dir=tmp_dir, seqs="seq2", no_masking=False),
            ]
        )
        mock_kraken.assert_called_once_with(db_dir=tmp_dir, all_kwargs=common_args)
        mock_bracken.assert_has_calls(
            [
                call(kraken2_db_dir=tmp_dir, threads=1, kmer_len=35, read_len=100),
                call(kraken2_db_dir=tmp_dir, threads=1, kmer_len=35, read_len=150),
            ]
        )
        mock_move.has_calls(
            [
                call(tmp_dir, str(kraken2_db.path), extension="k2d"),
                call(tmp_dir, str(bracken_db.path), extension="kmer_distrib"),
            ]
        )

    @patch("q2_annotate.kraken2.database._fetch_db_collection")
    @patch("q2_annotate.kraken2.database._move_db_files")
    def test_fetch_prebuilt_dbs(self, mock_move, mock_fetch):
        bracken_db = MagicMock(path="/path/to/bracken_db")
        kraken2_db = MagicMock(path="/path/to/kraken2_db")

        _fetch_prebuilt_dbs(bracken_db, kraken2_db, "some_collection", "/tmp")

        mock_fetch.assert_called_once_with(collection="some_collection", tmp_dir="/tmp")
        mock_move.assert_has_calls(
            [
                call("/tmp", str(kraken2_db.path), extension="k2d"),
                call("/tmp", str(bracken_db.path), extension="kmer_distrib"),
            ]
        )

    @patch("tempfile.TemporaryDirectory", return_value=MockTempDir())
    @patch("q2_annotate.kraken2.database.Kraken2DBDirectoryFormat")
    @patch("q2_annotate.kraken2.database.BrackenDBDirectoryFormat")
    @patch("q2_annotate.kraken2.database._fetch_prebuilt_dbs")
    def test_build_kraken_db_action_with_prebuilt(
        self, mock_fetch, mock_bracken, mock_kraken, mock_tmp
    ):
        fake_kraken_dir_fmt = Kraken2DBDirectoryFormat(self.get_data_path("db"), "r")
        mock_kraken.return_value = fake_kraken_dir_fmt
        fake_bracken_dir_fmt = BrackenDBDirectoryFormat(
            self.get_data_path("bracken-db"), "r"
        )
        mock_bracken.return_value = fake_bracken_dir_fmt

        annotate.actions.build_kraken_db(collection="viral")

        mock_fetch.assert_called_once_with(
            fake_bracken_dir_fmt,
            fake_kraken_dir_fmt,
            "viral",
            str(mock_tmp.return_value.name),
        )

    @patch("tempfile.TemporaryDirectory", return_value=MockTempDir())
    @patch("q2_annotate.kraken2.database.Kraken2DBDirectoryFormat")
    @patch("q2_annotate.kraken2.database.BrackenDBDirectoryFormat")
    @patch("q2_annotate.kraken2.database._build_dbs_from_seqs")
    def test_build_kraken_db_action_with_seqs(
        self, mock_build, mock_bracken, mock_kraken, mock_tmp
    ):
        seqs = Artifact.import_data("FeatureData[Sequence]", self.get_data_path("seqs"))
        fake_kraken_dir_fmt = Kraken2DBDirectoryFormat(self.get_data_path("db"), "r")
        mock_kraken.return_value = fake_kraken_dir_fmt
        fake_bracken_dir_fmt = BrackenDBDirectoryFormat(
            self.get_data_path("bracken-db"), "r"
        )
        mock_bracken.return_value = fake_bracken_dir_fmt

        annotate.actions.build_kraken_db(seqs=[seqs], threads=2, fast_build=True)

        exp_common_args = {
            "threads": 2,
            "kmer_len": 35,
            "minimizer_len": 31,
            "minimizer_spaces": 7,
            "no_masking": False,
            "max_db_size": 0,
            "use_ftp": False,
            "load_factor": 0.7,
            "fast_build": True,
            "read_len": [50, 75, 100, 150, 200, 250, 300],
            "kraken2_db": fake_kraken_dir_fmt,
            "bracken_db": fake_bracken_dir_fmt,
            "tmp": str(mock_tmp.return_value.name),
        }
        mock_build.assert_called_once_with(
            fake_bracken_dir_fmt,
            fake_kraken_dir_fmt,
            [ANY],
            str(mock_tmp.return_value.name),
            exp_common_args,
        )

    @patch("tempfile.TemporaryDirectory", return_value=MockTempDir())
    def test_build_kraken_db_action_with_error(self, mock_tmp):
        with self.assertRaisesRegex(
            ValueError, r"You need to either provide a list .+"
        ):
            annotate.actions.build_kraken_db()


class TestInspectKraken2Database(unittest.TestCase):
    package = "q2_annotate.kraken2.tests"

    @classmethod
    def setUpClass(cls):
        datadir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
        db_path = os.path.join(datadir, "simulated-sequences", "kraken2-db")
        db = Kraken2DBDirectoryFormat(db_path, "r")
        cls.report = inspect_kraken2_db(db)
        cls.report_view = cls.report.file.view(pd.DataFrame)

        cls.species_to_ncbi_id = {
            "Bacillus anthracis": 1392,
            "Mus musculus": 10090,
            "Staphylococcus aureus": 1280,
            "Staphylococcus epidermidis": 1282,
        }

    def test_format(self):
        self.assertIsInstance(self.report, Kraken2DBReportDirectoryFormat)
        self.report.validate()

    def test_report(self):
        self.assertGreater(len(self.report_view), 0)
        self.assertEqual(
            list(self.report_view.columns), list(Kraken2DBReportFormat.COLUMNS.keys())
        )

        self.report_view["name"] = self.report_view["name"].apply(
            lambda cell: cell.strip()
        )
        for species, ncbi_id in self.species_to_ncbi_id.items():
            self.assertIn(species, self.report_view["name"].values)
            self.assertIn(ncbi_id, self.report_view["taxon_id"].values)

        # now rows skipped when skipping header
        self.assertEqual("root", self.report_view.loc[0, "name"])


if __name__ == "__main__":
    unittest.main()
