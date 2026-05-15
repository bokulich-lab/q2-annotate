# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob

import shutil

import os
import tempfile
import unittest
import uuid
from q2_types.per_sample_sequences import ContigSequencesDirFmt, BAMDirFmt
from q2_types.per_sample_sequences import MultiFASTADirectoryFormat
from unittest.mock import MagicMock, patch, call

from qiime2.plugin.testing import TestPluginBase

from q2_annotate.semibin2.semibin2 import (
    bin_contigs_semibin2,
    _bin_contigs_semibin2,
    _run_semibin2,
    _process_sample,
)


class TestSemibin2(TestPluginBase):
    package = "q2_annotate.semibin2.tests"

    @patch("subprocess.run")
    def test_run_semibin2_ok(self, p1):
        fake_props = {"map": "/some/where/map.bam", "contigs": "/a/b/co.fa"}
        fake_args = [
            "--epochs",
            "16",
            "--min-len",
            "1000",
        ]

        with tempfile.TemporaryDirectory() as fake_loc:
            obs_fp = _run_semibin2(
                "samp1", fake_props, fake_loc, "single_easy_bin", fake_args
            )
            exp_fp = os.path.join(fake_loc, "samp1")

            self.assertEqual(exp_fp, obs_fp)

            exp_cmd = [
                "SemiBin2",
                "single_easy_bin",
                "--input-fasta",
                fake_props["contigs"],
                "--input-bam",
                fake_props["map"],
                "--output",
                os.path.join(fake_loc, "samp1", "bin"),
                "--compression",
                "none",
                "--verbose",
            ]
            exp_cmd.extend(fake_args)
            p1.assert_called_once_with(exp_cmd, check=True)

    @patch("tempfile.TemporaryDirectory")
    @patch("q2_annotate.semibin2.semibin2.uuid4")
    @patch("q2_annotate.semibin2.semibin2._run_semibin2")
    def test_process_sample(self, p1, p2, p3):
        fake_props = {
            "map": "some/where/samp1_alignment.bam",
            "contigs": "some/where/samp1_contigs.fasta",
        }
        fake_args = [
            "--compression",
            "none",
            "--verbose",
            "--epochs",
            "16",
            "--min-len",
            "1000",
        ]

        fake_temp_dir = tempfile.mkdtemp()
        p2.side_effect = [
            uuid.UUID("93fa393d-504a-4257-8a38-931af0cc4b1f"),
            uuid.UUID("a69ab2c1-b23a-400f-be45-7dfc35bf343c"),
        ]
        p3.return_value.__enter__.return_value = fake_temp_dir

        with tempfile.TemporaryDirectory() as fake_loc:
            p1.return_value = os.path.join(fake_loc, "bins", "samp1")

            # copy two expected bins to the new location
            samp1_bins_fp = self.get_data_path("bins-no-uuid/samp1")
            shutil.copytree(
                samp1_bins_fp,
                os.path.join(fake_loc, "bins", "samp1", "bin", "output_bins"),
                dirs_exist_ok=True,
            )

            _process_sample("samp1", fake_props, "single_easy_bin", fake_args, fake_loc)

            # find the newly formed bins
            obs_bins = set(
                [
                    x.split("/")[-1]
                    for x in glob.glob(os.path.join(fake_loc, "samp1", "*.fa"))
                ]
            )
            exp_bins = {
                "93fa393d-504a-4257-8a38-931af0cc4b1f.fa",
                "a69ab2c1-b23a-400f-be45-7dfc35bf343c.fa",
            }
            self.assertSetEqual(exp_bins, obs_bins)

    @patch("q2_annotate.semibin2.semibin2.ContigSequencesDirFmt")
    @patch("q2_annotate.semibin2.semibin2.MultiFASTADirectoryFormat")
    @patch("q2_annotate.semibin2.semibin2._process_sample")
    def test_bin_contigs_semibin2(self, p1, p2, _):
        input_contigs = self.get_data_path("contigs")
        input_maps = self.get_data_path("maps")
        contigs = ContigSequencesDirFmt(input_contigs, mode="r")
        maps = BAMDirFmt(input_maps, mode="r")

        args = [
            "--compression",
            "none",
            "--verbose",
            "--epochs",
            "16",
            "--min-len",
            "1000",
        ]

        mock_bins = MultiFASTADirectoryFormat(self.get_data_path("bins"), "r")
        p2.return_value = mock_bins

        obs_bins, obs_map = _bin_contigs_semibin2(
            contigs, maps, "single_easy_bin", args
        )

        self.assertIsInstance(obs_bins, MultiFASTADirectoryFormat)
        p1.assert_has_calls(
            [
                call(
                    "samp1",
                    {
                        "contigs": self.get_data_path("/contigs/samp1_contigs.fa"),
                        "map": self.get_data_path("/maps/samp1_alignment.bam"),
                    },
                    "single_easy_bin",
                    args,
                    str(mock_bins),
                ),
                call(
                    "samp2",
                    {
                        "contigs": self.get_data_path("/contigs/samp2_contigs.fa"),
                        "map": self.get_data_path("/maps/samp2_alignment.bam"),
                    },
                    "single_easy_bin",
                    args,
                    str(mock_bins),
                ),
            ]
        )

        # find the newly formed bins
        obs_bins = sorted(
            [
                sorted(
                    [
                        "/".join(x.split("/")[-2:])
                        for x in glob.glob(
                            os.path.join(str(obs_bins), f"samp{y}", "*.fa")
                        )
                    ]
                )
                for y in (1, 2)
            ]
        )
        exp_bins = [
            [
                "samp1/522775d4-b1c6-4ee3-8b47-cd990f17eb8b.fa",
                "samp1/684db670-6304-4f33-a0ea-7f570532e178.fa",
            ],
            [
                "samp2/37356c23-b8db-4bbe-b4c9-d35e1cef615b.fa",
                "samp2/51c19113-31f0-4e4c-bbb3-b9df26b949f3.fa",
            ],
        ]
        self.assertListEqual(exp_bins, obs_bins)

    @patch("q2_annotate.semibin2.semibin2.MultiFASTADirectoryFormat")
    @patch("q2_annotate.semibin2.semibin2._process_sample")
    def test_bin_contigs_metabat_no_mags(self, p1, p2):
        input_contigs = self.get_data_path("contigs")
        input_maps = self.get_data_path("maps")
        contigs = ContigSequencesDirFmt(input_contigs, mode="r")
        maps = BAMDirFmt(input_maps, mode="r")

        args = [
            "--compression",
            "none",
            "--verbose",
            "--epochs",
            "16",
            "--min-len",
            "1000",
        ]

        mock_bins = MultiFASTADirectoryFormat()
        p2.return_value = mock_bins

        with self.assertRaisesRegex(ValueError, "No MAGs were formed"):
            _bin_contigs_semibin2(contigs, maps, "single_easy_bin", args)

    @patch("q2_annotate.semibin2.semibin2._bin_contigs_semibin2")
    @patch("q2_annotate.semibin2.semibin2._process_common_input_params")
    def test_bin_contigs_semibin2_wrapper(self, p1, p2):
        p1.return_value = ["--epochs", "16"]
        p2.return_value = ("bins", {"contigA": "bin1"})

        contigs = MagicMock()
        bams = MagicMock()

        result = bin_contigs_semibin2(
            contigs=contigs,
            alignment_maps=bams,
            epochs=10,
            verbose=True,
        )

        p1.assert_called_once()

        p2.assert_called_once_with(
            contigs=contigs,
            alignment_maps=bams,
            mode="single_easy_bin",
            common_args=["--epochs", "16"],
        )

        assert result == ("bins", {"contigA": "bin1"})


if __name__ == "__main__":
    unittest.main()
