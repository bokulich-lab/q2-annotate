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
from pathlib import Path
import tempfile
import unittest
from q2_types.per_sample_sequences import ContigSequencesDirFmt, BAMDirFmt
from q2_types.per_sample_sequences import MultiFASTADirectoryFormat
from unittest.mock import patch, ANY, call

from qiime2.plugin.testing import TestPluginBase

from q2_annotate.metabat2.metabat2 import (
    _assert_samples,
    _sort_bams,
    _estimate_depth,
    _run_metabat2,
    _process_sample,
    _bin_contigs_metabat,
    _generate_contig_map,
)


class TestMetabat2(TestPluginBase):
    package = "q2_annotate.metabat2.tests"

    def test_assert_samples_ok(self):
        contigs_path = self.get_data_path("contigs")
        maps_path = self.get_data_path("maps")

        contigs = ContigSequencesDirFmt(contigs_path, mode="r")
        maps = BAMDirFmt(maps_path, mode="r")

        obs_samples = _assert_samples(contigs, maps)

        exp_samples = {
            "samp1": {
                "contigs": str(Path(contigs_path) / "samp1_contigs.fa"),
                "map": str(Path(maps_path) / "samp1_alignment.bam"),
            },
            "samp2": {
                "contigs": str(Path(contigs_path) / "samp2_contigs.fa"),
                "map": str(Path(maps_path) / "samp2_alignment.bam"),
            },
        }
        self.assertDictEqual(exp_samples, obs_samples)

    def test_assert_samples_uneven(self):
        contigs_path = self.get_data_path("contigs")
        with tempfile.TemporaryDirectory() as maps_path:
            map_path = Path(self.get_data_path("maps")) / "samp1_alignment.bam"
            shutil.copy(map_path, maps_path)

            contigs = ContigSequencesDirFmt(contigs_path, mode="r")
            maps = BAMDirFmt(maps_path, mode="r")

            with self.assertRaisesRegex(
                Exception,
                "Contigs and alignment maps should belong to the same sample"
                " set. You provided contigs for samples: samp1,samp2 but maps "
                "for samples: samp1. Please check your inputs and try again.",
            ):
                _assert_samples(contigs, maps)

    def test_assert_samples_non_matching(self):
        with tempfile.TemporaryDirectory() as tempdir:
            contigs_path = Path(tempdir) / "contigs-path"
            maps_path = Path(tempdir) / "maps-path"
            os.makedirs(contigs_path)
            os.makedirs(maps_path)

            contig_path = Path(self.get_data_path("contigs")) / "samp1_contigs.fa"
            map_path = Path(self.get_data_path("maps")) / "samp2_alignment.bam"

            shutil.copy(contig_path, contigs_path)
            shutil.copy(map_path, maps_path)

            contigs = ContigSequencesDirFmt(contigs_path, mode="r")
            maps = BAMDirFmt(maps_path, mode="r")

            with self.assertRaisesRegex(
                Exception,
                "Contigs and alignment maps should belong to the same sample"
                " set. You provided contigs for samples: samp1 but maps "
                "for samples: samp2. Please check your inputs and try again.",
            ):
                _assert_samples(contigs, maps)

    @patch("subprocess.run")
    def test_sort_bams_ok(self, p1):
        fake_props = {"map": "/some/where/map.bam", "fake_key": "abc"}

        obs_props = _sort_bams("samp1", fake_props, "/new/location")
        exp_props = {
            "map": "/new/location/samp1_alignment_sorted.bam",
            "fake_key": "abc",
        }

        self.assertDictEqual(exp_props, obs_props)
        p1.assert_called_once_with(
            [
                "samtools",
                "sort",
                fake_props["map"],
                "-o",
                "/new/location/samp1_alignment_sorted.bam",
            ],
            check=True,
        )

    @patch("subprocess.run")
    def test_estimate_depth_ok(self, p1):
        fake_props = {"map": "/some/where/map.bam", "fake_key": "abc"}

        obs_fp = _estimate_depth("samp1", fake_props, "/new/location")
        exp_fp = "/new/location/samp1_depth.txt"

        self.assertEqual(exp_fp, obs_fp)
        p1.assert_called_once_with(
            [
                "jgi_summarize_bam_contig_depths",
                "--outputDepth",
                exp_fp,
                fake_props["map"],
            ],
            check=True,
        )

    @patch("subprocess.run")
    def test_run_metabat2_ok(self, p1):
        fake_props = {"map": "/some/where/map.bam", "contigs": "/a/b/co.fa"}
        fake_args = ["--verbose", "--pTNF", "85", "--minContig", "2500"]

        with tempfile.TemporaryDirectory() as fake_loc:
            obs_fp = _run_metabat2(
                "samp1", fake_props, fake_loc, "/some/depth/file.txt", fake_args
            )
            exp_fp = os.path.join(fake_loc, "samp1")

            self.assertEqual(exp_fp, obs_fp)

            exp_cmd = [
                "metabat2",
                "-i",
                fake_props["contigs"],
                "-a",
                "/some/depth/file.txt",
                "-o",
                os.path.join(fake_loc, "samp1", "bin"),
                "--unbinned",
            ]
            exp_cmd.extend(fake_args)
            p1.assert_called_once_with(exp_cmd, check=True)

    @patch("tempfile.TemporaryDirectory")
    @patch("q2_annotate.metabat2.uuid4")
    @patch("q2_annotate.metabat2._sort_bams")
    @patch("q2_annotate.metabat2._estimate_depth")
    @patch("q2_annotate.metabat2._run_metabat2")
    def test_process_sample(self, p1, p2, p3, p4, p5):
        fake_props = {
            "map": "some/where/samp1_alignment.bam",
            "contigs": "some/where/samp1_contigs.fasta",
        }
        fake_props_mod = {
            "map": "some/where/samp1_alignment_sorted.bam",
            "contigs": "some/where/samp1_contigs.fasta",
        }
        fake_args = ["--verbose", "--minContig", "1500", "--minClsSize", "10000"]
        fake_unbinned = ContigSequencesDirFmt()

        p3.return_value = fake_props_mod
        p2.return_value = "some/where/samp1_depth.txt"
        p4.side_effect = [
            "522775d4-b1c6-4ee3-8b47-cd990f17eb8b",
            "684db670-6304-4f33-a0ea-7f570532e178",
            "37356c23-b8db-4bbe-b4c9-d35e1cef615b",
            "51c19113-31f0-4e4c-bbb3-b9df26b949f3",
        ]
        fake_temp_dir = tempfile.mkdtemp()
        p5.return_value.__enter__.return_value = fake_temp_dir

        with tempfile.TemporaryDirectory() as fake_loc:
            p1.return_value = os.path.join(fake_loc, "bins", "samp1")

            # copy two expected bins to the new location
            samp1_bins_fp = self.get_data_path("bins-no-uuid/samp1")
            shutil.copytree(
                samp1_bins_fp,
                os.path.join(fake_loc, "bins", "samp1"),
                dirs_exist_ok=True,
            )

            # copy expected unbinned contigs to the new location
            samp1_unbinned_fp = self.get_data_path("contigs/samp1_contigs.fa")
            os.makedirs(os.path.join(fake_loc, "samp1"))
            shutil.copy(
                samp1_unbinned_fp,
                os.path.join(fake_temp_dir, "bins", "samp1", "bin.unbinned.fa"),
            )

            _process_sample("samp1", fake_props, fake_args, fake_loc, fake_unbinned)

            # find the newly formed bins
            obs_bins = set(
                [
                    x.split("/")[-1]
                    for x in glob.glob(os.path.join(fake_loc, "samp1", "*.fa"))
                ]
            )
            exp_bins = {
                "522775d4-b1c6-4ee3-8b47-cd990f17eb8b.fa",
                "684db670-6304-4f33-a0ea-7f570532e178.fa",
            }
            self.assertSetEqual(exp_bins, obs_bins)

            # find the unbinned contigs
            obs_unbinned = set(
                [
                    x.split("/")[-1]
                    for x in glob.glob(os.path.join(fake_unbinned.path, "*.fa"))
                ]
            )
            exp_unbinned = {
                "samp1_contigs.fa",
            }
            self.assertSetEqual(exp_unbinned, obs_unbinned)

            p3.assert_called_once_with("samp1", fake_props, ANY)
            p2.assert_called_once_with("samp1", fake_props_mod, ANY)
            p1.assert_called_once_with(
                "samp1", fake_props_mod, ANY, "some/where/samp1_depth.txt", fake_args
            )

    @patch("q2_annotate.metabat2.ContigSequencesDirFmt")
    @patch("q2_annotate.metabat2.MultiFASTADirectoryFormat")
    @patch("q2_annotate.metabat2._process_sample")
    def test_bin_contigs_metabat(self, p1, p2, p3):
        input_contigs = self.get_data_path("contigs")
        input_maps = self.get_data_path("maps")
        contigs = ContigSequencesDirFmt(input_contigs, mode="r")
        maps = BAMDirFmt(input_maps, mode="r")

        args = ["--verbose", "--minContig", "1500", "--minClsSize", "10000"]

        mock_bins = MultiFASTADirectoryFormat(self.get_data_path("bins"), "r")
        p2.return_value = mock_bins

        mock_unbinned = ContigSequencesDirFmt(
            self.get_data_path("contigs/samp1_contigs.fa"), "r"
        )
        p3.return_value = mock_unbinned

        obs_bins, obs_map, obs_unbinned = _bin_contigs_metabat(contigs, maps, args)

        self.assertIsInstance(obs_bins, MultiFASTADirectoryFormat)
        p1.assert_has_calls(
            [
                call(
                    "samp1",
                    {
                        "contigs": self.get_data_path("/contigs/samp1_contigs.fa"),
                        "map": self.get_data_path("/maps/samp1_alignment.bam"),
                    },
                    args,
                    str(mock_bins),
                    str(mock_unbinned),
                ),
                call(
                    "samp2",
                    {
                        "contigs": self.get_data_path("/contigs/samp2_contigs.fa"),
                        "map": self.get_data_path("/maps/samp2_alignment.bam"),
                    },
                    args,
                    str(mock_bins),
                    str(mock_unbinned),
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

    @patch("q2_annotate.metabat2.MultiFASTADirectoryFormat")
    @patch("q2_annotate.metabat2._process_sample")
    def test_bin_contigs_metabat_no_mags(self, p1, p2):
        input_contigs = self.get_data_path("contigs")
        input_maps = self.get_data_path("maps")
        contigs = ContigSequencesDirFmt(input_contigs, mode="r")
        maps = BAMDirFmt(input_maps, mode="r")

        args = ["--verbose", "--minContig", "1500", "--minClsSize", "10000"]

        mock_bins = MultiFASTADirectoryFormat()
        p2.return_value = mock_bins

        with self.assertRaisesRegex(ValueError, "No MAGs were formed"):
            _bin_contigs_metabat(contigs, maps, args)

    def test_generate_contig_map(self):
        contigs = MultiFASTADirectoryFormat(self.get_data_path("bins-small"), "r")
        obs = _generate_contig_map(contigs)
        exp = {
            "684db670-6304-4f33-a0ea-7f570532e178": ["NODE_2", "NODE_6", "NODE_7"],
            "522775d4-b1c6-4ee3-8b47-cd990f17eb8b": ["NODE_8", "NODE_11"],
            "51c19113-31f0-4e4c-bbb3-b9df26b949f3": ["NODE_12", "NODE_13", "NODE_14"],
            "37356c23-b8db-4bbe-b4c9-d35e1cef615b": [
                "NODE_2",
                "NODE_6",
                "NODE_7",
                "NODE_15",
            ],
        }
        self.assertDictEqual(exp, obs)


if __name__ == "__main__":
    unittest.main()
