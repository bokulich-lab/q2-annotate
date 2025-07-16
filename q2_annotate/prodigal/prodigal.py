# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import copy as cp
from typing import Union, List
from .._utils import run_command
from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt, ContigSequencesDirFmt
from q2_types.genome_data import (
    LociDirectoryFormat,
    GenesDirectoryFormat,
    ProteinsDirectoryFormat,
)


def _run_prodigal(
    cmd: List[str],
    loci: LociDirectoryFormat,
    genes: GenesDirectoryFormat,
    proteins: ProteinsDirectoryFormat,
    path_to_input: str,
    _id: str,
    subdir: str = None,
):
    # If subdirectory is not None, append a "/" s.t. the command
    # below is defined correctly. Otw subdir = ""
    subdir = subdir + "/" if subdir else ""

    # Complete command and run
    cmd = cp.deepcopy(cmd)
    cmd.extend(
        [
            "-i",
            path_to_input,
            "-o",
            os.path.join(loci.path, f"{subdir}{_id}.gff"),
            "-a",
            os.path.join(proteins.path, f"{subdir}{_id}.fasta"),
            "-d",
            os.path.join(genes.path, f"{subdir}{_id}.fasta"),
        ]
    )
    run_command(cmd)


def predict_genes_prodigal(
    seqs: Union[MAGSequencesDirFmt, MultiMAGSequencesDirFmt, ContigSequencesDirFmt],
    translation_table_number: str = "11",
    mode: str = "single",
    closed: bool = False,
    no_shine_dalgarno: bool = False,
    mask: bool = False,
) -> (LociDirectoryFormat, GenesDirectoryFormat, ProteinsDirectoryFormat):

    # Instantiate output directories
    loci = LociDirectoryFormat()
    genes = GenesDirectoryFormat()
    proteins = ProteinsDirectoryFormat()

    base_cmd = ["prodigal", "-g", translation_table_number, "-f", "gff", "-p", mode]

    if closed:
        base_cmd.append("-c")
    if no_shine_dalgarno:
        base_cmd.append("-n")
    if mask:
        base_cmd.append("-m")

    if isinstance(seqs, MAGSequencesDirFmt):
        for mag_id, mag_fp in seqs.feature_dict().items():
            _run_prodigal(base_cmd, loci, genes, proteins, mag_fp, mag_id)
    elif isinstance(seqs, ContigSequencesDirFmt):
        for sample_id, contigs_fp in seqs.sample_dict().items():
            _run_prodigal(base_cmd, loci, genes, proteins, contigs_fp, sample_id)
    elif isinstance(seqs, MultiMAGSequencesDirFmt):
        for sample_id, mags_dict in seqs.sample_dict().items():
            # Make sample_id folders in output locations
            for output_object in [loci, genes, proteins]:
                os.makedirs(os.path.join(output_object.path, sample_id))

            # Run prodigal for each mag
            for mag_id, mag_fp in mags_dict.items():
                _run_prodigal(
                    base_cmd, loci, genes, proteins, mag_fp, mag_id, sample_id
                )

    # Return output directories
    return loci, genes, proteins
