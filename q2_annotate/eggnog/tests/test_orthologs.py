# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import shutil
import tempfile
from filecmp import dircmp
from unittest.mock import MagicMock, patch, ANY, call

import pandas as pd
import pandas.testing as pdt
import qiime2
from qiime2.plugin.testing import TestPluginBase

from q2_annotate.eggnog import (
    _eggnog_diamond_search,
    search_orthologs_hmmer,
    _eggnog_hmmer_search,
)
from q2_annotate.eggnog.orthologs.common import (
    _eggnog_search,
    _search_runner,
    _create_symlinks,
)
from q2_annotate.eggnog.types import EggnogHmmerIdmapDirectoryFmt
from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_types.genome_data import ProteinsDirectoryFormat
from q2_types.per_sample_sequences import ContigSequencesDirFmt, MultiMAGSequencesDirFmt
from q2_types.profile_hmms import PressedProfileHmmsDirectoryFmt
from q2_types.reference_db import DiamondDatabaseDirFmt


class TestHMMER(TestPluginBase):
    package = "q2_annotate.eggnog.tests"

    def setUp(self):
        super().setUp()
        self.idmap_artifact = qiime2.Artifact.import_data(
            "EggnogHmmerIdmap", self.get_data_path("idmap")
        )
        self.idmap = self.idmap_artifact.view(EggnogHmmerIdmapDirectoryFmt)

        self.pressed_hmm_artifact = qiime2.Artifact.import_data(
            "ProfileHMM[PressedProtein]", self.get_data_path("pressed_hmm")
        )
        self.pressed_hmm = self.pressed_hmm_artifact.view(
            PressedProfileHmmsDirectoryFmt
        )

        self.fastas_artifact = qiime2.Artifact.import_data(
            "GenomeData[Proteins]", self.get_data_path("fastas")
        )
        self.fastas = self.fastas_artifact.view(ProteinsDirectoryFormat)

        self.mags_artifact = qiime2.Artifact.import_data(
            "SampleData[MAGs]", self.get_data_path("mag-sequences-per-sample")
        )
        self.mags = self.mags_artifact.view(MultiMAGSequencesDirFmt)

        self.eggnog_hmmer_search = self.plugin.pipelines["search_orthologs_hmmer"]
        self._eggnog_hmmer_search = self.plugin.methods["_eggnog_hmmer_search"]
        self._eggnog_annotate = self.plugin.methods["_eggnog_annotate"]

    def test_eggnog_hmmer_search_pipeline(self):
        mock_action = MagicMock(
            side_effect=[
                lambda sequences, num_partitions: ({"mag1": {}, "mag2": {}},),
                lambda seq, pressed, idmap, fastas, num_cpus, db_in_memory: (0, 0, 0),
                lambda hits: ("collated_hits",),
                lambda collated_hits: ("collated_tables",),
                lambda collated_loci: ("collated_loci",),
            ]
        )
        mock_ctx = MagicMock(get_action=mock_action)
        obs = search_orthologs_hmmer(
            ctx=mock_ctx,
            seqs=self.mags_artifact,
            pressed_hmm_db=self.pressed_hmm_artifact,
            idmap=self.idmap_artifact,
            seed_alignments=self.fastas_artifact,
        )
        exp = ("collated_hits", "collated_tables", "collated_loci")
        self.assertTupleEqual(obs, exp)

    def test_symlink_files_to_target_dir(self):
        with tempfile.TemporaryDirectory() as tmp1:
            for dir in ["idmap", "fastas", "pressed_hmm"]:
                shutil.copytree(self.get_data_path(dir), tmp1, dirs_exist_ok=True)
            with tempfile.TemporaryDirectory() as tmp2:
                _create_symlinks(
                    [
                        self.get_data_path("pressed_hmm"),
                        self.get_data_path("idmap"),
                        self.get_data_path("fastas"),
                    ],
                    tmp2,
                )
                comp = dircmp(tmp1, tmp2)
                self.assertFalse(comp.diff_files or comp.left_only or comp.right_only)

    @patch("os.makedirs")
    @patch("tempfile.TemporaryDirectory")
    @patch("q2_annotate.eggnog.orthologs.hmmer._create_symlinks")
    @patch("q2_annotate.eggnog.orthologs.hmmer._eggnog_search")
    def test_eggnog_hmmer_search(
        self, mock_eggnog_search, mock_symlink, mock_tmpdir, mock_makedirs
    ):
        mock_tmpdir.return_value.__enter__.return_value = "tmp"
        mock_eggnog_search.return_value = (0, 1, 2)
        result, ft, loci = _eggnog_hmmer_search(
            seqs=self.mags,
            idmap=self.idmap,
            pressed_hmm_db=self.pressed_hmm,
            seed_alignments=self.fastas,
        )
        mock_symlink.assert_called_once_with(
            [self.pressed_hmm, self.idmap, self.fastas], "tmp/hmmer/1100069"
        )
        mock_eggnog_search.assert_called_once_with(
            self.mags, ANY, "tmp"  # partial() method not patchable or comparable
        )
        self.assertTupleEqual((result, ft, loci), (0, 1, 2))

    def test_eggnog_search_mags(self):
        sequences = MultiMAGSequencesDirFmt(
            self.get_data_path("mag-sequences-per-sample"), "r"
        )
        output_loc = self.get_data_path("hits")
        search_runner = MagicMock()

        result, ft, _ = _eggnog_search(sequences, search_runner, output_loc)
        result.validate()
        self.assertIsInstance(ft, pd.DataFrame)

        search_runner.assert_has_calls(
            [
                call(input_path=mag_fp, sample_label=mag_id)
                for sample_id, mags in sequences.sample_dict().items()
                for mag_id, mag_fp in mags.items()
            ]
        )

    def test_eggnog_search_contigs(self):
        sequences = ContigSequencesDirFmt(self.get_data_path("contig-sequences-1"), "r")
        output_loc = self.get_data_path("hits")
        search_runner = MagicMock()

        result, ft, _ = _eggnog_search(sequences, search_runner, output_loc)
        result.validate()
        self.assertIsInstance(ft, pd.DataFrame)

        search_runner.assert_has_calls(
            [
                call(input_path=contigs_fp, sample_label=sample_id)
                for sample_id, contigs_fp in sequences.sample_dict().items()
            ]
        )

    def test_eggnog_search_mags_derep(self):
        sequences = MAGSequencesDirFmt(self.get_data_path("mag-sequences"), "r")
        output_loc = self.get_data_path("hits")
        search_runner = MagicMock()

        result, ft, _ = _eggnog_search(sequences, search_runner, output_loc)
        result.validate()
        self.assertIsInstance(ft, pd.DataFrame)

        search_runner.assert_has_calls(
            [
                call(input_path=mag_fp, sample_label=mag_id)
                for mag_id, mag_fp in sequences.feature_dict().items()
            ]
        )

    @patch("subprocess.run")
    def test_search_runner(self, mock_run):
        _search_runner("a", "b", "c", "d", True, ["f", "g"])
        mock_run.called_once_with(
            [
                "emapper.py",
                "-i",
                "a",
                "-o",
                "b",
                "-m",
                "f",
                "g",
                "--itype",
                "metagenome",
                "--output_dir",
                "c",
                "--cpu",
                "d",
                "--no_annot",
                "--dbmem",
            ],
            check=True,
        )


class TestDiamond(TestPluginBase):
    package = "q2_annotate.eggnog.tests"

    def setUp(self):
        super().setUp()
        self.diamond_db_artifact = qiime2.Artifact.import_data(
            "ReferenceDB[Diamond]", self.get_data_path("random-db-1")
        )
        self.diamond_db = self.diamond_db_artifact.view(DiamondDatabaseDirFmt)

        self.eggnog_diamond_search = self.plugin.pipelines["search_orthologs_diamond"]
        self._eggnog_diamond_search = self.plugin.methods["_eggnog_diamond_search"]
        self._eggnog_annotate = self.plugin.methods["_eggnog_annotate"]

    def test_good_small_search_contigs(self):
        contigs = qiime2.Artifact.import_data(
            "SampleData[Contigs]", self.get_data_path("contig-sequences-1")
        ).view(ContigSequencesDirFmt)

        _, obs, _ = _eggnog_diamond_search(seqs=contigs, db=self.diamond_db)
        exp = pd.DataFrame(
            {"0": [1.0, 0.0], "2": [1.0, 0.0], "8": [0.0, 3.0]}, index=["s1", "s2"]
        )
        exp.columns.name = "sseqid"

        pdt.assert_frame_equal(obs, exp)

    def test_good_small_search_mags_derep(self):
        mags = qiime2.Artifact.import_data(
            "FeatureData[MAG]", self.get_data_path("mag-sequences")
        ).view(MAGSequencesDirFmt)

        _, obs, _ = _eggnog_diamond_search(seqs=mags, db=self.diamond_db)
        exp = pd.DataFrame(
            {"8": [3.0, 0.0], "0": [0.0, 1.0], "2": [0.0, 1.0]},
            index=[
                "194b1aca-9373-4298-ba5c-05b382d1f553",
                "e1ba1d20-b466-4fef-ae19-6bf9c5c63d6f",
            ],
        )
        exp.columns.name = "sseqid"

        pdt.assert_frame_equal(obs, exp)

    def test_good_small_search_mags(self):
        mags = qiime2.Artifact.import_data(
            "SampleData[MAGs]", self.get_data_path("mag-sequences-per-sample")
        ).view(MultiMAGSequencesDirFmt)

        _, obs, _ = _eggnog_diamond_search(seqs=mags, db=self.diamond_db)
        exp = pd.DataFrame(
            {
                "8": [3.0, 3.0, 0.0, 0.0],
                "0": [0.0, 0.0, 1.0, 1.0],
                "2": [0.0, 0.0, 1.0, 1.0],
            },
            index=[
                "194b1aca-9373-4298-ba5c-05b382d1f553",
                "8f40a3bc-14f0-426b-b2ba-7fb3b1cd0c02",
                "c6bd5123-35cf-473a-bebf-85cf544bcd48",
                "e1ba1d20-b466-4fef-ae19-6bf9c5c63d6f",
            ],
        )
        exp.columns.name = "sseqid"

        pdt.assert_frame_equal(obs, exp)

    def test_eggnog_search_parallel_contigs(self):
        contigs = qiime2.Artifact.import_data(
            "SampleData[Contigs]", self.get_data_path("contig-sequences-1")
        )

        with self.test_config:
            _, parallel, _ = self.eggnog_diamond_search.parallel(
                contigs, self.diamond_db_artifact
            )._result()

        _, single, _ = self._eggnog_diamond_search(
            seqs=contigs, db=self.diamond_db_artifact
        )

        parallel = parallel.view(pd.DataFrame)
        single = single.view(pd.DataFrame)

        pdt.assert_frame_equal(parallel, single)

    def test_eggnog_search_parallel_mags_derep(self):
        mags = qiime2.Artifact.import_data(
            "FeatureData[MAG]", self.get_data_path("mag-sequences")
        )

        with self.test_config:
            _, parallel, _ = self.eggnog_diamond_search.parallel(
                mags, self.diamond_db_artifact
            )._result()

        _, single, _ = self._eggnog_diamond_search(
            seqs=mags, db=self.diamond_db_artifact
        )

        parallel = parallel.view(pd.DataFrame)
        single = single.view(pd.DataFrame)

        pdt.assert_frame_equal(parallel, single)

    def test_eggnog_search_parallel_mags(self):
        mags = qiime2.Artifact.import_data(
            "SampleData[MAGs]", self.get_data_path("mag-sequences-per-sample")
        )

        with self.test_config:
            _, parallel, _ = self.eggnog_diamond_search.parallel(
                mags, self.diamond_db_artifact
            )._result()

        _, single, _ = self._eggnog_diamond_search(
            seqs=mags, db=self.diamond_db_artifact
        )

        parallel = parallel.view(pd.DataFrame)
        single = single.view(pd.DataFrame)

        pdt.assert_frame_equal(parallel, single)
