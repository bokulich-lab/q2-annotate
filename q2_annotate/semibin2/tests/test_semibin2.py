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

from q2_annotate.semibin2.semibin2 import (
    _assert_samples,
    _sort_bams,
    _run_semibin2_single,
    _run_semibin2_multi,
    _concatenate_contigs_with_semibin2,
    _process_semibin2_outputs,
    _process_sample,
    _bin_contigs_semibin2,
    _generate_contig_map,
)


class TestSemiBin2(TestPluginBase):
    package = "q2_annotate.semibin2.tests"

    def test_assert_samples_ok(self):
        """Test that sample validation works correctly."""
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
        """Test that uneven sample sets raise appropriate errors."""
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
        """Test that non-matching sample names raise appropriate errors."""
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
        """Test BAM sorting functionality."""
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
    def test_run_semibin2_single_ok(self, p1):
        """Test SemiBin2 single-sample execution."""
        fake_props = {"map": "/some/where/map.bam", "contigs": "/a/b/co.fa"}
        fake_args = ["--threads", "4", "--min-len", "1000"]

        with tempfile.TemporaryDirectory() as fake_loc:
            obs_fp = _run_semibin2_single(
                "samp1", fake_props, fake_loc, fake_args, multi_sample=False
            )
            exp_fp = os.path.join(fake_loc, "samp1")

            self.assertEqual(exp_fp, obs_fp)

            exp_cmd = [
                "SemiBin2",
                "single_easy_bin",
                "-i", fake_props["contigs"],
                "-b", fake_props["map"],
                "-o", exp_fp,
                "--verbose",
            ]
            exp_cmd.extend(fake_args)
            p1.assert_called_once_with(exp_cmd, check=True)

    @patch("subprocess.run")
    def test_concatenate_contigs_with_semibin2(self, p1):
        """Test SemiBin2 concatenate_fasta helper function."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create fake contig files and sample set within temp directory
            fake_sample_set = {
                "samp1": {
                    "contigs": os.path.join(temp_dir, "samp1.fa"), 
                    "map": os.path.join(temp_dir, "samp1.bam")
                },
                "samp2": {
                    "contigs": os.path.join(temp_dir, "samp2.fa"), 
                    "map": os.path.join(temp_dir, "samp2.bam")
                },
            }

            # Create fake contig files for testing
            for samp in fake_sample_set:
                contig_file = fake_sample_set[samp]["contigs"]
                with open(contig_file, "w") as f:
                    f.write(f">contig1_{samp}\nATCG\n>contig2_{samp}\nGCTA\n")

            obs_combined = _concatenate_contigs_with_semibin2(fake_sample_set, temp_dir)
            exp_combined = os.path.join(temp_dir, "concatenated_contigs", "concatenated.fa")

            self.assertEqual(exp_combined, obs_combined)

            # Verify SemiBin2 concatenate_fasta was called correctly
            expected_call_args = p1.call_args[0][0]
            self.assertEqual(expected_call_args[0], "SemiBin2")
            self.assertEqual(expected_call_args[1], "concatenate_fasta")
            self.assertIn("-i", expected_call_args)
            self.assertIn("-o", expected_call_args)
            
            # Verify that the -o parameter is the directory, not the file
            o_index = expected_call_args.index("-o")
            output_dir = expected_call_args[o_index + 1]
            self.assertEqual(output_dir, os.path.join(temp_dir, "concatenated_contigs"))
            
            # Verify that the individual contig files are passed directly
            i_index = expected_call_args.index("-i")
            contig_files = expected_call_args[i_index + 1:o_index]
            expected_contig_files = [
                os.path.join(temp_dir, "samp1.fa"),
                os.path.join(temp_dir, "samp2.fa")
            ]
            self.assertEqual(sorted(contig_files), sorted(expected_contig_files))

    @patch("q2_annotate.semibin2.semibin2._concatenate_contigs_with_semibin2")
    @patch("subprocess.run")
    def test_run_semibin2_multi_ok(self, p1, p_concat):
        """Test SemiBin2 multi-sample execution."""
        fake_args = ["--threads", "4", "--min-len", "1000"]

        with tempfile.TemporaryDirectory() as temp_dir:
            fake_sample_set = {
                "samp1": {
                    "contigs": os.path.join(temp_dir, "samp1.fa"), 
                    "map": os.path.join(temp_dir, "samp1.bam")
                },
                "samp2": {
                    "contigs": os.path.join(temp_dir, "samp2.fa"), 
                    "map": os.path.join(temp_dir, "samp2.bam")
                },
            }

            combined_contigs = os.path.join(temp_dir, "concatenated_contigs", "concatenated.fa")
            p_concat.return_value = combined_contigs

            obs_fp = _run_semibin2_multi(fake_sample_set, temp_dir, fake_args)
            exp_fp = os.path.join(temp_dir, "multi_sample_bins")

            self.assertEqual(exp_fp, obs_fp)

            # Check that concatenate_contigs_with_semibin2 was called
            p_concat.assert_called_once_with(fake_sample_set, temp_dir)

            # Check that the command was called with the expected arguments
            expected_call_args = p1.call_args[0][0]
            self.assertEqual(expected_call_args[0], "SemiBin2")
            self.assertEqual(expected_call_args[1], "multi_easy_bin")
            self.assertIn("-i", expected_call_args)
            self.assertIn(combined_contigs, expected_call_args)
            self.assertIn("-b", expected_call_args)
            self.assertIn("-o", expected_call_args)
            self.assertIn("--verbose", expected_call_args)

    def test_process_semibin2_outputs(self):
        """Test processing of SemiBin2 output files."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            # Create fake SemiBin2 output structure
            bins_dp = os.path.join(tmp_dir, "semibin_output")
            recluster_dir = os.path.join(bins_dp, "output_recluster_bins")
            os.makedirs(recluster_dir)

            # Create fake bin files
            bin1_path = os.path.join(recluster_dir, "bin.1.fa")
            bin2_path = os.path.join(recluster_dir, "bin.2.fa")
            
            with open(bin1_path, "w") as f:
                f.write(">contig1\nATCG\n")
            with open(bin2_path, "w") as f:
                f.write(">contig2\nGCTA\n")

            result_loc = os.path.join(tmp_dir, "results")
            unbinned_loc = os.path.join(tmp_dir, "unbinned")
            os.makedirs(result_loc)
            os.makedirs(unbinned_loc)

            _process_semibin2_outputs(bins_dp, "samp1", result_loc, unbinned_loc)

            # Check that bins were moved and renamed
            sample_bins_dir = os.path.join(result_loc, "samp1")
            self.assertTrue(os.path.exists(sample_bins_dir))
            
            bin_files = glob.glob(os.path.join(sample_bins_dir, "*.fa"))
            self.assertEqual(len(bin_files), 2)
            
            # Check that original bin files no longer exist
            self.assertFalse(os.path.exists(bin1_path))
            self.assertFalse(os.path.exists(bin2_path))

    @patch("tempfile.TemporaryDirectory")
    @patch("q2_annotate.semibin2.uuid4")
    @patch("q2_annotate.semibin2._sort_bams")
    @patch("q2_annotate.semibin2._run_semibin2_single")
    @patch("q2_annotate.semibin2._process_semibin2_outputs")
    def test_process_sample(self, p_outputs, p_run, p_sort, p_uuid, p_temp):
        """Test single sample processing pipeline."""
        fake_props = {
            "map": "some/where/samp1_alignment.bam",
            "contigs": "some/where/samp1_contigs.fasta",
        }
        fake_props_mod = {
            "map": "some/where/samp1_alignment_sorted.bam",
            "contigs": "some/where/samp1_contigs.fasta",
        }
        fake_args = ["--threads", "4", "--min-len", "1000"]

        p_sort.return_value = fake_props_mod
        p_run.return_value = "/tmp/bins_output"
        fake_temp_dir = tempfile.mkdtemp()
        p_temp.return_value.__enter__.return_value = fake_temp_dir

        _process_sample("samp1", fake_props, fake_args, "/result/loc", "/unbinned/loc", False)

        p_sort.assert_called_once_with("samp1", fake_props, ANY)
        p_run.assert_called_once_with("samp1", fake_props_mod, ANY, fake_args, False)
        p_outputs.assert_called_once_with("/tmp/bins_output", "samp1", "/result/loc", "/unbinned/loc")

    @patch("q2_annotate.semibin2.ContigSequencesDirFmt")
    @patch("q2_annotate.semibin2.MultiFASTADirectoryFormat")
    @patch("q2_annotate.semibin2._process_sample")
    def test_bin_contigs_semibin2_single_sample(self, p_process, p_multifasta, p_contigs):
        """Test the main binning function in single-sample mode."""
        input_contigs = self.get_data_path("contigs")
        input_maps = self.get_data_path("maps")
        contigs = ContigSequencesDirFmt(input_contigs, mode="r")
        maps = BAMDirFmt(input_maps, mode="r")

        args = ["--threads", "4", "--min-len", "1000"]

        mock_bins = MultiFASTADirectoryFormat(self.get_data_path("bins"), "r")
        p_multifasta.return_value = mock_bins

        mock_unbinned = ContigSequencesDirFmt(
            self.get_data_path("contigs/samp1_contigs.fa"), "r"
        )
        p_contigs.return_value = mock_unbinned

        obs_bins, obs_map, obs_unbinned = _bin_contigs_semibin2(
            contigs, maps, args, multi_sample=False
        )

        self.assertIsInstance(obs_bins, MultiFASTADirectoryFormat)
        p_process.assert_has_calls(
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
                    False,
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
                    False,
                ),
            ]
        )

    @patch("q2_annotate.semibin2.MultiFASTADirectoryFormat")
    @patch("q2_annotate.semibin2._process_sample")
    def test_bin_contigs_semibin2_no_mags(self, p_process, p_multifasta):
        """Test that an error is raised when no MAGs are formed."""
        input_contigs = self.get_data_path("contigs")
        input_maps = self.get_data_path("maps")
        contigs = ContigSequencesDirFmt(input_contigs, mode="r")
        maps = BAMDirFmt(input_maps, mode="r")

        args = ["--threads", "4", "--min-len", "1000"]

        mock_bins = MultiFASTADirectoryFormat()
        p_multifasta.return_value = mock_bins

        with self.assertRaisesRegex(ValueError, "No MAGs were formed"):
            _bin_contigs_semibin2(contigs, maps, args, multi_sample=False)

    def test_generate_contig_map(self):
        """Test contig mapping generation."""
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