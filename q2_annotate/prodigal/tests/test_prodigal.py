# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
from q2_annotate.prodigal.prodigal import predict_genes_prodigal
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt, ContigSequencesDirFmt
from unittest.mock import patch, call
from q2_types.genome_data import (
    LociDirectoryFormat,
    GenesDirectoryFormat,
    ProteinsDirectoryFormat,
)


class TestBUSCO(TestPluginBase):
    package = "q2_annotate.prodigal.tests"

    @patch("subprocess.run")
    def test_run_prodigal_feature_data_1_mag(self, subp_run):
        # Run prodigal with dummy data
        p = self.get_data_path("mags/dir_with_1_mag")
        mags = MAGSequencesDirFmt(path=p, mode="r")
        loci, genes, proteins = predict_genes_prodigal(seqs=mags)

        # Check that output is correct type
        self.assertIsInstance(loci, LociDirectoryFormat)
        self.assertIsInstance(genes, GenesDirectoryFormat)
        self.assertIsInstance(proteins, ProteinsDirectoryFormat)

        # Get names of fasta files from test data dir
        fasta_file = [
            os.path.splitext(file)[0]
            for file in os.listdir(mags.path)
            if file.endswith(".fa") or file.endswith(".fasta")
        ]

        # Unpack list
        fasta_file = fasta_file[0]

        # Assert that patch was called once
        subp_run.assert_called_once_with(
            [
                "prodigal",
                "-g",
                "11",
                "-f",
                "gff",
                "-i",
                os.path.join(mags.path, f"{fasta_file}.fasta"),
                "-o",
                os.path.join(loci.path, f"{fasta_file}.gff"),
                "-a",
                os.path.join(proteins.path, f"{fasta_file}.fasta"),
                "-d",
                os.path.join(genes.path, f"{fasta_file}.fasta"),
            ],
            check=True,
        )

    @patch("subprocess.run")
    def test_run_prodigal_feature_data_3_mag(self, subp_run):
        # Run prodigal with dummy data
        p = self.get_data_path("mags/dir_with_3_mag")
        mags = MAGSequencesDirFmt(path=p, mode="r")
        loci, genes, proteins = predict_genes_prodigal(seqs=mags)

        # Check that output is correct type
        self.assertIsInstance(loci, LociDirectoryFormat)
        self.assertIsInstance(genes, GenesDirectoryFormat)
        self.assertIsInstance(proteins, ProteinsDirectoryFormat)

        # Get names of fasta files from test data dir
        fasta_files = [
            os.path.splitext(file)[0]
            for file in os.listdir(mags.path)
            if file.endswith(".fa") or file.endswith(".fasta")
        ]

        # Define calls
        three_calls = [
            call(
                [
                    "prodigal",
                    "-g",
                    "11",
                    "-f",
                    "gff",
                    "-i",
                    os.path.join(mags.path, f"{fasta_file}.fasta"),
                    "-o",
                    os.path.join(loci.path, f"{fasta_file}.gff"),
                    "-a",
                    os.path.join(proteins.path, f"{fasta_file}.fasta"),
                    "-d",
                    os.path.join(genes.path, f"{fasta_file}.fasta"),
                ],
                check=True,
            )
            for fasta_file in fasta_files
        ]

        # Assert that patch was called 3 times
        subp_run.assert_has_calls(three_calls, any_order=True)

    @patch("subprocess.run")
    def test_run_prodigal_sample_data(self, subp_run):
        p = self.get_data_path("mags")
        mags = MultiMAGSequencesDirFmt(path=p, mode="r")
        loci, genes, prot = predict_genes_prodigal(seqs=mags)

        # Check that output is correct type
        self.assertIsInstance(loci, LociDirectoryFormat)
        self.assertIsInstance(genes, GenesDirectoryFormat)
        self.assertIsInstance(prot, ProteinsDirectoryFormat)

        # Get names of fasta files from test data dir
        calls = []
        for sample in os.listdir(mags.path):
            for fasta_file in os.listdir(f"{mags.path}/{sample}"):
                file_id = os.path.splitext(fasta_file)[0]
                # Define calls
                calls.append(
                    call(
                        [
                            "prodigal",
                            "-g",
                            "11",
                            "-f",
                            "gff",
                            "-i",
                            os.path.join(mags.path, sample, f"{file_id}.fasta"),
                            "-o",
                            os.path.join(loci.path, f"{sample}/{file_id}.gff"),
                            "-a",
                            os.path.join(prot.path, f"{sample}/{file_id}.fasta"),
                            "-d",
                            os.path.join(genes.path, f"{sample}/{file_id}.fasta"),
                        ],
                        check=True,
                    )
                )

        # Assert that patch was called 3 times
        subp_run.assert_has_calls(calls, any_order=True)

    @patch("subprocess.run")
    def test_run_prodigal_contigs(self, subp_run):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), mode="r")
        loci, genes, prot = predict_genes_prodigal(seqs=contigs)

        subp_run.assert_called_once_with(
            [
                "prodigal",
                "-g",
                "11",
                "-f",
                "gff",
                "-i",
                os.path.join(contigs.path, "sample1_contigs.fasta"),
                "-o",
                os.path.join(loci.path, "sample1.gff"),
                "-a",
                os.path.join(prot.path, "sample1.fasta"),
                "-d",
                os.path.join(genes.path, "sample1.fasta"),
            ],
            check=True,
        )
