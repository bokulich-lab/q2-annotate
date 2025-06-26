# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import csv
from qiime2.core.exceptions import ValidationError
from qiime2.plugin import model
from q2_types.feature_data import AlignedProteinFASTAFormat


class BUSCOResultsFormat(model.TextFileFormat):
    HEADER = [
        "mag_id",
        "sample_id",
        "input_file",
        "dataset",
        "complete",
        "single",
        "duplicated",
        "fragmented",
        "missing",
        "n_markers",
        "scaffold_n50",
        "contigs_n50",
        "percent_gaps",
        "scaffolds",
        "length",
        "completeness",
        "contamination",
    ]
    HEADER_2 = HEADER[:-2]

    def _validate(self, n_records=None):
        with self.open() as fh:
            reader = csv.reader(fh, delimiter="\t")
            headers = next(reader)

            if set(headers) == set(self.HEADER):
                expected = self.HEADER
            elif set(headers) == set(self.HEADER_2):
                expected = self.HEADER_2
            else:
                raise ValidationError(
                    f"Invalid header: {headers}, expected: {self.HEADER}, "
                    '"completness" and "contamination" columns are optional.'
                )

            for i, row in enumerate(reader, start=2):
                if len(row) != len(expected):
                    raise ValidationError(
                        f"Line {i} has {len(row)} columns, " f"expected {len(expected)}"
                    )

                if n_records is not None and i - 1 >= n_records:
                    break

    def _validate_(self, level):
        record_count_map = {"min": 100, "max": None}
        self._validate(record_count_map[level])


BUSCOResultsDirectoryFormat = model.SingleFileDirectoryFormat(
    "BUSCOResultsDirectoryFormat", "busco_results.tsv", BUSCOResultsFormat
)


class BuscoGenericTextFileFmt(model.TextFileFormat):
    def _validate_(self, level):
        pass


class BuscoGenericBinaryFileFmt(model.BinaryFileFormat):
    def _validate_(self, level):
        pass


class BuscoDatabaseDirFmt(model.DirectoryFormat):
    # File collections for text files.
    # Optional because some of those are not present in some lineages.
    (
        ancestral,
        ancestral_variants,
        dataset,
        hmms,
        lengths_cutoff,
        links_to_ODB,
        ogs_id,
        refseq_db_md5,
        scores_cutoff,
        species,
    ) = [
        model.FileCollection(
            rf"lineages\/.+\/{pattern}", format=BuscoGenericTextFileFmt, optional=True
        )
        for pattern in [
            r"ancestral$",
            r"ancestral_variants$",
            r"dataset\.cfg$",
            r"hmms\/.+\.hmm$",
            r"lengths_cutoff$",
            r"links_to_ODB.+\.txt$",
            r"info\/ogs\.id\.info$",
            r"refseq_db\.faa\.gz\.md5",
            r"scores_cutoff$",
            r"info\/species\.info$",
        ]
    ]

    # Placement files. Optional because they are not in virus DB
    (
        list_of_reference_markers,
        mapping_taxid_lineage,
        mapping_taxids_busco_dataset_name,
        tree,
        tree_metadata,
    ) = [
        model.FileCollection(
            rf"placement_files\/{pattern}",
            format=BuscoGenericTextFileFmt,
            optional=True,
        )
        for pattern in [
            r"list_of_reference_markers\..+\.txt$",
            r"mapping_taxid-lineage\..+\.txt$",
            r"mapping_taxids-busco_dataset_name\..+\.txt$",
            r"tree\..+\.nwk$",
            r"tree_metadata\..+\.txt$",
        ]
    ]

    # Others
    supermatrix_aln = model.FileCollection(
        r"placement_files\/supermatrix\.aln\..+\.faa$",
        format=AlignedProteinFASTAFormat,
        optional=True,
    )
    prfls = model.FileCollection(
        r"lineages\/.+\/prfl\/.+\.prfl$", format=BuscoGenericTextFileFmt, optional=True
    )
    version_file = model.File("file_versions.tsv", format=BuscoGenericTextFileFmt)
    refseq_db = model.FileCollection(
        r"lineages\/.+refseq_db\.faa(\.gz)?",
        format=BuscoGenericBinaryFileFmt,
        optional=True,
    )
    information = model.FileCollection(
        r"information\/.+\.txt$", format=BuscoGenericTextFileFmt, optional=True
    )
    missing_parasitic = model.File(
        r"lineages\/.+\/missing_in_parasitic\.txt$",
        format=BuscoGenericTextFileFmt,
        optional=True,
    )
    no_hits = model.File(
        r"lineages\/.+\/no_hits$", format=BuscoGenericTextFileFmt, optional=True
    )

    def _path_maker(self, name):
        return str(name)

    def __init__(self, path, mode):
        super().__init__(path, mode)

        # Overwrite path maker methods for all file collections
        for var_name, var_value in vars(self.__class__).items():
            if isinstance(var_value, model.FileCollection):
                var_value.set_path_maker(self._path_maker)

    def _validate_(self, level):
        pass
