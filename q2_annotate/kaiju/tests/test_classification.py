# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json
import os
import unittest
from subprocess import CalledProcessError
from unittest.mock import patch, Mock, ANY, MagicMock, call

import numpy as np
import pandas as pd
import qiime2
from pandas._testing import assert_frame_equal
from q2_types.per_sample_sequences import (
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
    CasavaOneEightSingleLanePerSampleDirFmt,
    ContigSequencesDirFmt,
)

from q2_annotate.kaiju.classification import (
    _construct_feature_table,
    _rename_taxon,
    _encode_unclassified_ids,
    _fix_id_types,
    _process_kaiju_reports,
    _classify_kaiju,
    classify_kaiju,
)
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import annotate


class TestKaijuClassification(TestPluginBase):
    package = "q2_annotate.kaiju.tests"

    def setUp(self):
        super().setUp()
        self.ctx = MagicMock()

        self.mock_classify_kaiju = MagicMock(
            side_effect=[("table1", "taxonomy1"), ("table2", "taxonomy2")]
        )
        self.ctx.get_action.side_effect = lambda domain, action_name: {
            ("annotate", "_classify_kaiju"): self.mock_classify_kaiju,
            ("feature_table", "merge"): self.mock_merge,
            ("feature_table", "merge_taxa"): self.mock_merge_taxa,
            ("demux", "partition_samples_single"): self.mock_partition_single,
            ("demux", "partition_samples_paired"): self.mock_partition_paired,
            ("types", "partition_contigs"): self.mock_partition_contigs,
            ("types", "partition_sample_data_mags"): self.mock_partition_mags,
        }[domain, action_name]

        # Additional action mocks
        self.mock_merge = MagicMock(return_value=("merged_table",))
        self.mock_merge_taxa = MagicMock(return_value=("merged_taxa",))
        self.mock_partition_single = MagicMock()
        self.mock_partition_paired = MagicMock()
        self.mock_partition_contigs = MagicMock()
        self.mock_partition_mags = MagicMock()

        self.mock_partition_single.return_value = {"part1": Mock()}
        self.mock_partition_paired.return_value = {"part1": Mock()}
        self.mock_partition_contigs.return_value = {"part1": Mock()}
        self.mock_partition_mags.return_value = {"part1": Mock()}
        self.classify_kaiju = annotate.pipelines.classify_kaiju
        with open(self.get_data_path("taxa-map.json")) as f:
            self.taxa_map = json.load(f)

    def test_construct_feature_table(self):
        obs_table, obs_taxonomy = _construct_feature_table(
            table_fp=self.get_data_path("kaiju-table.tsv")
        )

        exp_table = pd.read_csv(self.get_data_path("kaiju-ft-ok.csv"), index_col=0)
        exp_table.columns.name = "taxon_id"
        exp_taxonomy = pd.read_csv(
            self.get_data_path("kaiju-taxonomy-ok.csv"), index_col=0
        )

        assert_frame_equal(exp_table, obs_table)
        assert_frame_equal(exp_taxonomy, obs_taxonomy)

    def test_rename_taxon_unspecified(self):
        obs_taxon = _rename_taxon("1236", self.taxa_map)
        exp_taxon = (
            "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;"
            "o__Unspecified;f__Unspecified;g__Unspecified;s__Unspecified"
        )

        self.assertEqual(exp_taxon, obs_taxon)

    def test_rename_taxon_full(self):
        obs_taxon = _rename_taxon("1773", self.taxa_map)
        exp_taxon = (
            "d__Bacteria;p__Actinobacteria;c__Actinomycetia;"
            "o__Corynebacteriales;f__Mycobacteriaceae;"
            "g__Mycobacterium;s__Mycobacterium tuberculosis"
        )

        self.assertEqual(exp_taxon, obs_taxon)

    def test_rename_taxon_unclassified(self):
        obs_taxon = _rename_taxon("dW5jbGFz", self.taxa_map)
        exp_taxon = "d__unclassified"

        self.assertEqual(exp_taxon, obs_taxon)

    def test_encode_unclassified(self):
        df_in = pd.DataFrame(
            {
                "taxon_name": ["Alpha", "Beta", "Alpha", "Gamma"],
                "taxon_id": [1, 2, 3, 4],
            }
        )

        obs_table = _encode_unclassified_ids(df_in, "Alph")
        exp_table = pd.DataFrame(
            {
                "taxon_name": ["Alpha", "Beta", "Alpha", "Gamma"],
                "taxon_id": ["QWxwaGE=", 2, "QWxwaGE=", 4],
            }
        )

        assert_frame_equal(exp_table, obs_table)

    def test_fix_id_types(self):
        df_in = pd.DataFrame(
            {
                "taxon_name": [
                    "cannot be assigned",
                    "Beta",
                    "unclassified",
                    "unclassified",
                ],
                "taxon_id": [np.nan, 21.0, np.nan, np.nan],
            }
        )

        obs_table = _fix_id_types(df_in)
        exp_table = pd.DataFrame(
            {
                "taxon_name": [
                    "cannot be assigned",
                    "Beta",
                    "unclassified",
                    "unclassified",
                ],
                "taxon_id": ["Y2Fubm90", 21, "dW5jbGFz", "dW5jbGFz"],
            }
        )

        assert_frame_equal(exp_table, obs_table, check_dtype=False)

    @patch("subprocess.run")
    @patch("q2_annotate.kaiju.classification._construct_feature_table")
    def test_process_kaiju_reports_c_float(self, p1, p2):
        open(os.path.join(str(self.temp_dir.name), "sample1.out"), "w").close()
        open(os.path.join(str(self.temp_dir.name), "sample2.out"), "w").close()
        args = {
            "r": "species",
            "db": Mock(path=self.temp_dir.name),
            "exp": True,
            "u": True,
            "c": 0.6,
        }

        _process_kaiju_reports(self.temp_dir.name, args)

        exp_cmd = [
            "kaiju2table",
            "-v",
            "-o",
            f"{self.temp_dir.name}/results.tsv",
            "-r",
            "species",
            "-t",
            f"{self.temp_dir.name}/nodes.dmp",
            "-n",
            f"{self.temp_dir.name}/names.dmp",
            "-l",
            "superkingdom,phylum,class,order,family,genus,species",
            "-e",
            "-u",
            "-m",
            "0.6",
            f"{self.temp_dir.name}/sample1.out",
            f"{self.temp_dir.name}/sample2.out",
        ]

        p1.assert_called_once_with(f"{self.temp_dir.name}/results.tsv")
        p2.assert_called_once_with(exp_cmd, check=True)

    @patch("subprocess.run")
    @patch("q2_annotate.kaiju.classification._construct_feature_table")
    def test_process_kaiju_reports_c_int(self, p1, p2):
        open(os.path.join(str(self.temp_dir.name), "sample1.out"), "w").close()
        open(os.path.join(str(self.temp_dir.name), "sample2.out"), "w").close()
        args = {
            "r": "species",
            "db": Mock(path=self.temp_dir.name),
            "exp": True,
            "u": True,
            "c": 2,
        }

        _process_kaiju_reports(self.temp_dir.name, args)

        exp_cmd = [
            "kaiju2table",
            "-v",
            "-o",
            f"{self.temp_dir.name}/results.tsv",
            "-r",
            "species",
            "-t",
            f"{self.temp_dir.name}/nodes.dmp",
            "-n",
            f"{self.temp_dir.name}/names.dmp",
            "-l",
            "superkingdom,phylum,class,order,family,genus,species",
            "-e",
            "-u",
            "-c",
            "2",
            f"{self.temp_dir.name}/sample1.out",
            f"{self.temp_dir.name}/sample2.out",
        ]

        p1.assert_called_once_with(f"{self.temp_dir.name}/results.tsv")
        p2.assert_called_once_with(exp_cmd, check=True)

    @patch("subprocess.run")
    @patch("q2_annotate.kaiju.classification._process_kaiju_reports")
    def test_classify_kaiju_single(self, p1, p2):
        seqs = SingleLanePerSampleSingleEndFastqDirFmt(
            self.get_data_path("single-end"), "r"
        )
        r_out = r".*sample1.out,.*sample2.out"
        r_in = r".*reads1_R1.fastq.gz,.*reads2_R1.fastq.gz"
        self.classify_kaiju_test_helper(p1, p2, seqs, r_in, r_out)

    @patch("subprocess.run")
    @patch("q2_annotate.kaiju.classification._process_kaiju_reports")
    def test_classify_kaiju_paired(self, p1, p2):
        seqs = SingleLanePerSamplePairedEndFastqDirFmt(
            self.get_data_path("paired-end"), "r"
        )
        r_out = r".*sample1.out,.*sample2.out"
        r_in = r".*reads1_R2.fastq.gz,.*reads2_R2.fastq.gz"
        r_paired = r".*reads1_R1.fastq.gz,.*reads2_R1.fastq.gz"

        self.classify_kaiju_test_helper(p1, p2, seqs, r_in, r_out, r_paired)

    @patch("subprocess.run")
    @patch("q2_annotate.kaiju.classification._process_kaiju_reports")
    def test_classify_kaiju_contigs(self, p1, p2):
        seqs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        r_out = r".*sample1.out,.*sample2.out"
        r_in = r".*sample1_contigs.fasta,.*sample2_contigs.fasta"
        self.classify_kaiju_test_helper(p1, p2, seqs, r_in, r_out)

    def classify_kaiju_test_helper(self, p1, p2, seqs, r_in, r_out, r_paired=None):
        db_path = self.temp_dir.name
        open(os.path.join(db_path, "kaiju_123.fmi"), "w").close()
        p1.return_value = [pd.DataFrame(), pd.DataFrame()]

        with patch("tempfile.TemporaryDirectory"):
            _classify_kaiju(
                seqs=seqs,
                db=Mock(path=self.temp_dir.name),
                z=3,
                a="greedy",
                e=2,
                m=10,
                s=66,
                evalue=0,
                x=False,
                r="class",
                c=0.1,
                exp=True,
                u=True,
            )

        exp_cmd = [
            "kaiju-multi",
            "-v",
            "-z",
            "3",
            "-a",
            "greedy",
            "-e",
            "2",
            "-m",
            "10",
            "-s",
            "66",
            "-E",
            "0",
            "-X",
            "-t",
            os.path.join(db_path, "nodes.dmp"),
            "-f",
            os.path.join(db_path, "kaiju_123.fmi"),
            "-i",
            ANY,
            "-o",
            ANY,
        ]
        obs_cmd = p2.call_args.args[0]

        if r_paired:
            self.assertRegex(obs_cmd[-5], r_paired)
            exp_cmd[21:2] = ["-j", ANY]

        self.assertRegex(obs_cmd[-1], r_out)
        self.assertRegex(obs_cmd[-3], r_in)
        p2.assert_called_once_with(exp_cmd, check=True)

    @patch("subprocess.run", side_effect=CalledProcessError(1, "hello"))
    def test_classify_kaiju_exception(self, p1):
        seqs = SingleLanePerSampleSingleEndFastqDirFmt(
            self.get_data_path("single-end"), "r"
        )
        open(os.path.join(self.temp_dir.name, "kaiju_123.fmi"), "w").close()

        with self.assertRaisesRegex(Exception, r"\(return code 1\), please inspect"):
            _classify_kaiju(
                seqs=seqs,
                db=Mock(path=self.temp_dir.name),
                z=3,
                a="greedy",
                e=2,
                m=10,
                s=66,
                evalue=0,
                x=False,
                r="class",
                c=0.1,
                exp=True,
                u=True,
            )

    def test_classify_kaiju_single_partition_single_end(self):
        fake_seqs = qiime2.Artifact.import_data(
            "SampleData[SequencesWithQuality]",
            self.get_data_path("single-end-casava"),
            CasavaOneEightSingleLanePerSampleDirFmt,
        )
        fake_db = Mock()

        self.mock_partition_single.side_effect = [({0: "part1"},)]

        out_table, out_taxonomy = classify_kaiju(
            self.ctx, fake_seqs, fake_db, num_partitions=1
        )

        self.ctx.get_action.assert_any_call("annotate", "_classify_kaiju")
        self.ctx.get_action.assert_any_call("feature_table", "merge")
        self.ctx.get_action.assert_any_call("feature_table", "merge_taxa")
        self.ctx.get_action.assert_any_call("demux", "partition_samples_single")

        self.mock_partition_single.assert_called_once_with(fake_seqs, 1)
        self.mock_partition_paired.assert_not_called()
        self.mock_classify_kaiju.assert_called_once_with(
            "part1",
            fake_db,
            z=1,
            a="greedy",
            e=3,
            m=11,
            s=65,
            evalue=0.01,
            x=True,
            r="species",
            c=0.0,
            exp=False,
            u=False,
        )
        self.mock_merge.assert_called_once_with(["table1"])
        self.mock_merge_taxa.assert_called_once_with(["taxonomy1"])
        self.assertEqual("merged_table", out_table)
        self.assertEqual("merged_taxa", out_taxonomy)

    def test_classify_kaiju_single_partition_paired_end(self):
        fake_seqs = qiime2.Artifact.import_data(
            "SampleData[PairedEndSequencesWithQuality]",
            self.get_data_path("paired-end-casava"),
            CasavaOneEightSingleLanePerSampleDirFmt,
        )
        fake_db = Mock()

        self.mock_partition_paired.side_effect = [({0: "part1"},)]

        out_table, out_taxonomy = classify_kaiju(
            self.ctx, fake_seqs, fake_db, num_partitions=1
        )

        self.ctx.get_action.assert_any_call("annotate", "_classify_kaiju")
        self.ctx.get_action.assert_any_call("feature_table", "merge")
        self.ctx.get_action.assert_any_call("feature_table", "merge_taxa")
        self.ctx.get_action.assert_any_call("demux", "partition_samples_paired")

        self.mock_partition_single.assert_not_called()
        self.mock_partition_paired.assert_called_once_with(fake_seqs, 1)
        self.mock_classify_kaiju.assert_called_once_with(
            "part1",
            fake_db,
            z=1,
            a="greedy",
            e=3,
            m=11,
            s=65,
            evalue=0.01,
            x=True,
            r="species",
            c=0.0,
            exp=False,
            u=False,
        )
        self.mock_merge.assert_called_once_with(["table1"])
        self.mock_merge_taxa.assert_called_once_with(["taxonomy1"])
        self.assertEqual("merged_table", out_table)
        self.assertEqual("merged_taxa", out_taxonomy)

    def test_classify_kaiju_multiple_partitions_single_end(self):
        fake_seqs = qiime2.Artifact.import_data(
            "SampleData[SequencesWithQuality]",
            self.get_data_path("single-end-casava"),
            CasavaOneEightSingleLanePerSampleDirFmt,
        )
        fake_db = Mock()

        self.mock_partition_single.side_effect = [({0: "part1", 1: "part2"},)]

        out_table, out_taxonomy = classify_kaiju(
            self.ctx, fake_seqs, fake_db, num_partitions=2
        )

        self.ctx.get_action.assert_any_call("annotate", "_classify_kaiju")
        self.ctx.get_action.assert_any_call("feature_table", "merge")
        self.ctx.get_action.assert_any_call("feature_table", "merge_taxa")
        self.ctx.get_action.assert_any_call("demux", "partition_samples_single")

        self.mock_partition_single.assert_called_once_with(fake_seqs, 2)
        self.mock_partition_paired.assert_not_called()
        self.mock_classify_kaiju.assert_has_calls(
            [
                call(
                    "part1",
                    fake_db,
                    z=1,
                    a="greedy",
                    e=3,
                    m=11,
                    s=65,
                    evalue=0.01,
                    x=True,
                    r="species",
                    c=0.0,
                    exp=False,
                    u=False,
                ),
                call(
                    "part2",
                    fake_db,
                    z=1,
                    a="greedy",
                    e=3,
                    m=11,
                    s=65,
                    evalue=0.01,
                    x=True,
                    r="species",
                    c=0.0,
                    exp=False,
                    u=False,
                ),
            ]
        )
        self.mock_merge.assert_called_once_with(["table1", "table2"])
        self.mock_merge_taxa.assert_called_once_with(["taxonomy1", "taxonomy2"])
        self.assertEqual("merged_table", out_table)
        self.assertEqual("merged_taxa", out_taxonomy)

    def test_classify_kaiju_multiple_partitions_paired_end(self):
        fake_seqs = qiime2.Artifact.import_data(
            "SampleData[PairedEndSequencesWithQuality]",
            self.get_data_path("paired-end-casava"),
            CasavaOneEightSingleLanePerSampleDirFmt,
        )
        fake_db = Mock()

        self.mock_partition_paired.side_effect = [({0: "part1", 1: "part2"},)]

        out_table, out_taxonomy = classify_kaiju(
            self.ctx, fake_seqs, fake_db, num_partitions=2
        )

        self.ctx.get_action.assert_any_call("annotate", "_classify_kaiju")
        self.ctx.get_action.assert_any_call("feature_table", "merge")
        self.ctx.get_action.assert_any_call("feature_table", "merge_taxa")
        self.ctx.get_action.assert_any_call("demux", "partition_samples_paired")

        self.mock_partition_single.assert_not_called()
        self.mock_partition_paired.assert_called_once_with(fake_seqs, 2)
        self.mock_classify_kaiju.assert_has_calls(
            [
                call(
                    "part1",
                    fake_db,
                    z=1,
                    a="greedy",
                    e=3,
                    m=11,
                    s=65,
                    evalue=0.01,
                    x=True,
                    r="species",
                    c=0.0,
                    exp=False,
                    u=False,
                ),
                call(
                    "part2",
                    fake_db,
                    z=1,
                    a="greedy",
                    e=3,
                    m=11,
                    s=65,
                    evalue=0.01,
                    x=True,
                    r="species",
                    c=0.0,
                    exp=False,
                    u=False,
                ),
            ]
        )
        self.mock_merge.assert_called_once_with(["table1", "table2"])
        self.mock_merge_taxa.assert_called_once_with(["taxonomy1", "taxonomy2"])
        self.assertEqual("merged_table", out_table)
        self.assertEqual("merged_taxa", out_taxonomy)

    def test_classify_kaiju_single_partition_contigs(self):
        fake_seqs = qiime2.Artifact.import_data(
            "SampleData[Contigs]", self.get_data_path("contigs"), ContigSequencesDirFmt
        )
        fake_db = Mock()

        self.mock_partition_contigs.side_effect = [({0: "part1"},)]

        out_table, out_taxonomy = classify_kaiju(
            self.ctx, fake_seqs, fake_db, num_partitions=1
        )

        self.ctx.get_action.assert_any_call("annotate", "_classify_kaiju")
        self.ctx.get_action.assert_any_call("feature_table", "merge")
        self.ctx.get_action.assert_any_call("feature_table", "merge_taxa")
        self.ctx.get_action.assert_any_call("types", "partition_contigs")

        self.mock_partition_contigs.assert_called_once_with(fake_seqs, 1)
        self.mock_classify_kaiju.assert_called_once_with(
            "part1",
            fake_db,
            z=1,
            a="greedy",
            e=3,
            m=11,
            s=65,
            evalue=0.01,
            x=True,
            r="species",
            c=0.0,
            exp=False,
            u=False,
        )
        self.mock_merge.assert_called_once_with(["table1"])
        self.mock_merge_taxa.assert_called_once_with(["taxonomy1"])
        self.assertEqual("merged_table", out_table)
        self.assertEqual("merged_taxa", out_taxonomy)


if __name__ == "__main__":
    unittest.main()
