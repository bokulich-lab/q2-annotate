# ----------------------------------------------------------------------------
# Copyright (c) 2023-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import platform
import shutil
import tempfile
import unittest
from subprocess import CalledProcessError
from unittest.mock import patch

import pandas as pd
from pandas._testing import assert_frame_equal
from qiime2.plugin.testing import TestPluginBase

from q2_annotate.kraken2.bracken import (
    _assert_read_lens_available, _run_bracken_one_sample, _estimate_bracken,
    estimate_bracken
)
from q2_types.kraken2 import (BrackenDBDirectoryFormat,
                              Kraken2ReportDirectoryFormat)


class MockTempDir(tempfile.TemporaryDirectory):
    pass


class TestBracken(TestPluginBase):
    package = "q2_annotate.kraken2.tests"

    def setUp(self):
        super().setUp()
        self.bracken_db_dir = 'fake/db/dir'
        self.kwargs = {
            'threshold': 5,
            'read_len': 150,
            'level': 'S',
            'include_unclassified': True
        }

        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_assert_read_lens_avail(self):
        db = BrackenDBDirectoryFormat(self.get_data_path('bracken-db'), 'r')

        obs = _assert_read_lens_available(
            bracken_db=db, read_len=150
        )
        self.assertIsNone(obs)

    def test_assert_read_lens_not_avail(self):
        db = BrackenDBDirectoryFormat(self.get_data_path('bracken-db'), 'r')

        with self.assertRaisesRegex(
                ValueError,
                r'read length \(250\) is not available.*'
                r'The available values are\: 100\, 150\.'
        ):
            _assert_read_lens_available(bracken_db=db, read_len=250)

    @patch('q2_annotate.kraken2.bracken.run_command')
    def test_run_bracken_one_sample(self, p1):
        shutil.copyfile(
            self.get_data_path('bracken-report/sample1.bracken.output.txt'),
            os.path.join(self.temp_dir, 'sample1.bracken.output.txt')
        )

        kraken2_report_fp = self.get_data_path(
            'kraken2-reports-select/samples/sample1.report.txt'
        )
        bracken_report_dir = self.get_data_path('bracken-report')
        obs_table = _run_bracken_one_sample(
            bracken_db=self.bracken_db_dir,
            kraken2_report_fp=kraken2_report_fp,
            bracken_report_dir=bracken_report_dir,
            tmp_dir=self.temp_dir,
            threshold=self.kwargs['threshold'],
            read_len=self.kwargs['read_len'],
            level=self.kwargs['level'],
        )
        exp_table = pd.read_csv(
            self.get_data_path('bracken-report/sample1.table.csv'),
            index_col=0
        )
        exp_table["taxonomy_id"] = exp_table["taxonomy_id"].astype(str)

        assert_frame_equal(obs_table, exp_table)
        p1.assert_called_once_with(
            cmd=[
                "bracken", "-d", self.bracken_db_dir,
                "-i",
                kraken2_report_fp,
                "-o",
                os.path.join(self.temp_dir, "sample1.bracken.output.txt"),
                "-w",
                os.path.join(
                    bracken_report_dir, "sample1.report.txt"
                ),
                "-t", str(self.kwargs['threshold']),
                "-r", str(self.kwargs['read_len']),
                "-l", self.kwargs['level']
            ], verbose=True
        )

    @patch(
        'q2_annotate.kraken2.bracken.run_command',
        side_effect=CalledProcessError(returncode=123, cmd='bracken')
    )
    def test_run_bracken_one_sample_error(self, p1):
        kraken2_report_fp = tempfile.NamedTemporaryFile()
        bracken_report_dir = self.get_data_path('bracken-report')

        with self.assertRaisesRegex(Exception, 'return code 123'):
            _run_bracken_one_sample(
                bracken_db=self.bracken_db_dir,
                kraken2_report_fp=kraken2_report_fp.name,
                bracken_report_dir=bracken_report_dir,
                tmp_dir=self.temp_dir,
                threshold=self.kwargs['threshold'],
                read_len=self.kwargs['read_len'],
                level=self.kwargs['level'],
            )

    @patch('q2_annotate.kraken2.bracken._run_bracken_one_sample')
    def test_estimate_bracken(self, p1):
        kraken_reports = Kraken2ReportDirectoryFormat(
            self.get_data_path('reports-mags'), 'r'
        )
        bracken_db = BrackenDBDirectoryFormat()

        tables = [
            pd.read_csv(
                self.get_data_path('bracken-report/sample1.table.csv'),
                index_col=0
            ),
            pd.read_csv(
                self.get_data_path('bracken-report/sample2.table.csv'),
                index_col=0
            )
        ]
        p1.side_effect = tables

        obs_table, obs_reports = _estimate_bracken(
            kraken2_reports=kraken_reports,
            db=bracken_db,
            threshold=self.kwargs['threshold'],
            read_len=self.kwargs['read_len'],
            level=self.kwargs['level'],
        )
        exp_table = pd.read_csv(
            self.get_data_path('bracken-report/samples-merged.csv'),
            index_col=0
        )
        exp_table.columns = [int(x) for x in exp_table.columns]
        exp_table.columns.name = "taxonomy_id"
        exp_table.index = pd.Index(exp_table.index, name='sample_id')

        assert_frame_equal(obs_table, exp_table)
        self.assertIsInstance(obs_reports, Kraken2ReportDirectoryFormat)

    @patch('platform.system', return_value='Linux')
    @patch('q2_annotate.kraken2.bracken._assert_read_lens_available')
    @patch('q2_annotate.kraken2.bracken._estimate_bracken')
    def test_estimate_bracken_with_unclassified(self, p1, p2, p3):
        kraken_reports = Kraken2ReportDirectoryFormat(
            self.get_data_path('bracken-report-with-unclassified/'
                               'kraken-reports'), 'r'
        )
        bracken_db = BrackenDBDirectoryFormat()

        table = pd.read_csv(self.get_data_path(
            'bracken-report-with-unclassified/samples-merged.csv'
        ), index_col=0)
        p1.return_value = (table, kraken_reports)

        obs_reports, obs_taxonomy, obs_table = estimate_bracken(
            kraken2_reports=kraken_reports,
            db=bracken_db,
            threshold=self.kwargs['threshold'],
            read_len=self.kwargs['read_len'],
            level=self.kwargs['level'],
            include_unclassified=self.kwargs['include_unclassified']
        )
        exp_table = pd.read_csv(
            self.get_data_path('bracken-report-with-unclassified/'
                               'samples-merged-corrected.csv'),
            index_col=0
        )
        exp_table.index = pd.Index(exp_table.index, name='sample_id')
        exp_table['0'] = exp_table['0'].astype(int)
        exp_taxonomy = pd.read_csv(
            self.get_data_path('bracken-report-with-unclassified/'
                               'taxonomy.csv'),
            index_col=0
        )
        exp_taxonomy.index = exp_taxonomy.index.astype(str)

        assert_frame_equal(obs_table, exp_table)
        assert_frame_equal(obs_taxonomy, exp_taxonomy)
        self.assertIsInstance(obs_reports, Kraken2ReportDirectoryFormat)

    def test_estimate_bracken_linux_vs_osx(self):
        from qiime2.plugins import annotate
        try:
            annotate.methods.estimate_bracken
        except AttributeError:
            if platform.system() == "Darwin":
                pass
            else:
                raise


if __name__ == "__main__":
    unittest.main()
