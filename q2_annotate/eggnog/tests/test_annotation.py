# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import filecmp

import pandas as pd
import pandas.testing as pdt
import qiime2
from qiime2.plugin.testing import TestPluginBase

from q2_annotate.eggnog import _eggnog_annotate, extract_annotations
from q2_annotate.eggnog.annotation import _extract_generic, _filter, extraction_methods
from q2_types.genome_data import (
    OrthologAnnotationDirFmt,
    SeedOrthologDirFmt,
    OrthologFileFmt,
)
from q2_types.reference_db import EggnogRefDirFmt


class TestAnnotate(TestPluginBase):
    package = "q2_annotate.eggnog.tests"

    def setUp(self):
        super().setUp()
        self.eggnog_db = EggnogRefDirFmt(self.get_data_path("eggnog_db/"), mode="r")
        self.eggnog_db_artifact = qiime2.Artifact.import_data(
            "ReferenceDB[Eggnog]", self.get_data_path("eggnog_db/")
        )
        self.eggnog_annotate = self.plugin.pipelines["map_eggnog"]
        self._eggnog_annotate = self.plugin.methods["_eggnog_annotate"]

    def test_small_good_hits(self):
        seed_orthologs = SeedOrthologDirFmt(self.get_data_path("good_hits/"), mode="r")

        obs_obj = _eggnog_annotate(eggnog_hits=seed_orthologs, db=self.eggnog_db)

        exp_fp = self.get_data_path("expected/test_output.emapper.annotations")
        exp = OrthologFileFmt(exp_fp, mode="r").view(pd.DataFrame)

        objs = list(obs_obj.annotations.iter_views(OrthologFileFmt))
        self.assertEqual(len(objs), 1)
        df = objs[0][1].view(pd.DataFrame)
        pdt.assert_frame_equal(df, exp)

    def test_eggnog_annotate_parallel(self):
        orthologs = qiime2.Artifact.import_data(
            "SampleData[Orthologs]", self.get_data_path("good_hits/")
        )

        with self.test_config:
            (parallel,) = self.eggnog_annotate.parallel(
                orthologs, self.eggnog_db_artifact
            )._result()

        (single,) = self._eggnog_annotate(
            eggnog_hits=orthologs, db=self.eggnog_db_artifact
        )

        parallel = parallel.view(OrthologAnnotationDirFmt)
        single = single.view(OrthologAnnotationDirFmt)

        compare_dir = filecmp.dircmp(parallel.path, single.path)
        self.assertEqual(len(compare_dir.common), 1)

        # TODO: add exact file comparison


class TestAnnotationExtraction(TestPluginBase):
    package = "q2_annotate.eggnog.tests"

    def setUp(self):
        super().setUp()
        self.eggnog_ftf = pd.DataFrame(
            data={
                "ortholog1": [2, 0, 1, 10],
                "ortholog2": [0, 0, 1, 7],
                "ortholog3": [5, 2, 1, 0],
                "ortholog4": [1, 0, 1, 0],
                "ortholog5": [1, 3, 0, 9],
            },
            index=[
                "b9b4ab71-8e5f-48d7-bb23-df2726df1393",
                "62e07985-2556-435c-9e02-e7f94b8df07d",
                "1e9ffc02-0847-4f2c-b1e2-3965a4a78b15",
                "ab4f5ff0-45a1-41c9-9711-620765d5e92c",
            ],
        )
        self.mags_tpm = pd.DataFrame(
            data={
                "b9b4ab71-8e5f-48d7-bb23-df2726df1393": [0, 150.0, 10.0],
                "62e07985-2556-435c-9e02-e7f94b8df07d": [20.5, 15.0, 0.01],
                "1e9ffc02-0847-4f2c-b1e2-3965a4a78b15": [10.0, 1.5, 20.0],
                "ab4f5ff0-45a1-41c9-9711-620765d5e92c": [0.5, 0, 15.0],
            },
            index=["sample1", "sample2", "sample3"],
        )
        self.annotations = OrthologAnnotationDirFmt(
            self.get_data_path("annotations/"), mode="r"
        )
        self.annotation_df = pd.read_csv(
            self.get_data_path(
                "annotations/1e9ffc02-0847-4f2c-b1e2-3965a4a78b15."
                "emapper.annotations"
            ),
            sep="\t",
            skiprows=4,
            index_col=0,
        ).iloc[:-3, :]
        self.df = pd.DataFrame(
            {
                "evalue": [0, 0.1, 0.2, 0.3],
                "score": [400.0, 300.0, 200.0, 100.0],
                "column": ["val1", "val2", "val3", "val4"],
            },
            index=["s1", "s2", "s3", "s4"],
        )

    def test_extract_annotations(self):
        obs_ft = extract_annotations(
            ortholog_annotations=self.annotations, annotation="cog"
        )
        exp_ft = pd.DataFrame(
            data={
                "L": [2.0, 0.0, 10.0, 3.0],
                "F": [1.0, 2.0, 0.0, 5.0],
                "A": [1.0, 0.0, 7.0, 0.0],
            },
            index=pd.Index(
                [
                    "1e9ffc02-0847-4f2c-b1e2-3965a4a78b15",
                    "62e07985-2556-435c-9e02-e7f94b8df07d",
                    "ab4f5ff0-45a1-41c9-9711-620765d5e92c",
                    "b9b4ab71-8e5f-48d7-bb23-df2726df1393",
                ],
                name="id",
            ),
        )
        pd.testing.assert_frame_equal(obs_ft, exp_ft)

    def test_extract_annotations_not_implemented(self):
        with self.assertRaisesRegex(
            NotImplementedError, "Annotation 'hello' not supported."
        ):
            extract_annotations(
                ortholog_annotations=self.annotations, annotation="hello"
            )

    def test_extract_generic(self):
        obs = _extract_generic(
            self.annotation_df, "EC", lambda x: pd.Series(x.split("."))
        )
        exp = pd.Series([3, 3, 2, 2, 2], ["6", "3", "5", "4", "12"], name="count")
        pd.testing.assert_series_equal(obs, exp)

    def test_extract_cog(self):
        col, func = extraction_methods["cog"]
        obs = _extract_generic(self.annotation_df, col, func)
        exp = pd.Series([2, 1, 1], ["L", "F", "A"], name="count")
        pd.testing.assert_series_equal(obs, exp)

    def test_extract_kegg_ko(self):
        col, func = extraction_methods["kegg_ko"]
        obs = _extract_generic(self.annotation_df, col, func)
        exp = pd.Series(
            [1, 1, 1, 1], ["K01955", "K02621", "K16898", "K03722"], name="count"
        )
        pd.testing.assert_series_equal(obs, exp)

    def test_extract_kegg_pathway(self):
        col, func = extraction_methods["kegg_pathway"]
        obs = _extract_generic(self.annotation_df, col, func)
        exp = pd.Series([1, 1, 1], ["map00240", "map00250", "map01100"], name="count")
        pd.testing.assert_series_equal(obs, exp)

    def test_extract_kegg_module(self):
        col, func = extraction_methods["kegg_module"]
        obs = _extract_generic(self.annotation_df, col, func)
        exp = pd.Series([1], ["M00051"], name="count")
        pd.testing.assert_series_equal(obs, exp)

    def test_extract_kegg_reaction(self):
        col, func = extraction_methods["kegg_reaction"]
        obs = _extract_generic(self.annotation_df, col, func)
        exp = pd.Series(
            [1, 1, 1, 1, 1],
            ["R00256", "R00575", "R01395", "R10948", "R10949"],
            name="count",
        )
        pd.testing.assert_series_equal(obs, exp)

    def test_extract_brite(self):
        col, func = extraction_methods["brite"]
        obs = _extract_generic(self.annotation_df, col, func)
        exp = pd.Series(
            [4, 4, 2, 1, 1, 1, 1, 1],
            [
                "ko00000",
                "ko01000",
                "ko03400",
                "ko00001",
                "ko00002",
                "ko02048",
                "ko03032",
                "ko03036",
            ],
            name="count",
        )
        pd.testing.assert_series_equal(obs, exp)

    def test_extract_caz(self):
        col, func = extraction_methods["caz"]
        obs = _extract_generic(self.annotation_df, col, func)
        exp = pd.Series(name="count")
        pd.testing.assert_series_equal(obs, exp, check_index=False, check_dtype=False)

    def test_extract_ec(self):
        col, func = extraction_methods["ec"]
        obs = _extract_generic(self.annotation_df, col, func)
        exp = pd.Series([2, 1], ["3.6.4.12", "6.3.5.5"], name="count")
        pd.testing.assert_series_equal(obs, exp)

    def test_filter(self):
        obs = _filter(self.df, 0.2, 300.0)
        exp = pd.DataFrame(
            {"evalue": [0.0, 0.1], "score": [400.0, 300.0], "column": ["val1", "val2"]},
            index=["s1", "s2"],
        )
        pd.testing.assert_frame_equal(obs, exp)

    def test_filter_empty(self):
        with self.assertRaisesRegex(ValueError, " resulted in an empty table"):
            _filter(self.df, 0.1, 500.0)
