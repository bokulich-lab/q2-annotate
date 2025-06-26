# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import re
from qiime2.plugin import model
from qiime2.core.exceptions import ValidationError


class EggnogHmmerIdmapFileFmt(model.TextFileFormat):
    def _validate_(self, level):
        with open(str(self), "r") as file:
            # Set the number of rows to be parsed
            max_lines = {"min": 100, "max": float("inf")}[level]
            lines = file.readlines()
            for i, line in enumerate(lines, 1):
                # Check number of lines parsed so far
                if i > max_lines:
                    break

                # Validate line
                if not re.match(r"^(\d+) ([A-Z0-9]+)$", line):
                    raise ValidationError(
                        f"Invalid line {i}.\n"
                        f"{line} \n"
                        "Expected index and an alphanumeric code separated "
                        "by a single space."
                    )

                # Check index is equal to line number
                idx, code = line.rstrip("\n").split(sep=" ")
                if not idx == str(i):
                    raise ValidationError(
                        f"Invalid line {i}.\n"
                        f"{line} \n"
                        f"Expected index {i} but got {idx} instead.\n"
                    )


class EggnogHmmerIdmapDirectoryFmt(model.DirectoryFormat):
    idmap = model.File(r".*\.hmm\.idmap", format=EggnogHmmerIdmapFileFmt)
