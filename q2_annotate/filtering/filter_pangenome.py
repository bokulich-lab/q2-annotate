# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import os
import shutil
import subprocess
import tempfile

from q2_annotate._utils import run_command

EBI_SERVER_URL = (
    "ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ127/" "ERZ12792464/hprc-v1.0-pggb.gfa.gz"
)


def _fetch_and_extract_pangenome(uri: str, dest_dir: str):
    """
    Fetches and extracts the human pangenome GFA file.

    Args:
        uri (str): The URI of the genome to fetch. Should be in the form
                    ftp://host/path/to/file.
        dest_dir (str): The directory where the data will be saved.
    """
    filename = os.path.basename(uri)
    dest_fp = os.path.join(dest_dir, filename)

    try:
        print("Fetching the GFA file...")
        run_command(["wget", uri, "-q", "-O", dest_fp])
    except Exception as e:
        raise Exception(
            "Unable to connect to the server. Please try again later. "
            f"The error was: {e}"
        )

    print("Download finished. Extracting files...")
    run_command(["gunzip", dest_fp])


def _extract_fasta_from_gfa(gfa_fp: str, fasta_fp: str):
    """
    Extracts a FASTA file from a GFA file using gfatools.

    This function runs a subprocess calling 'gfatools' to convert a GFA file
    into a FASTA file. If the conversion is successful, the original GFA
    file is removed. Otherwise, an exception is raised.

    Args:
        gfa_fp (str): The file path to the input GFA file.
        fasta_fp (str): The file path where the output FASTA will be saved.
    """
    cmd = ["gfatools", "gfa2fa", gfa_fp]
    with open(fasta_fp, "w") as f_out:
        try:
            subprocess.run(cmd, stdout=f_out)
        except Exception as e:
            raise Exception(
                f"Failed to extract the fasta file from the GFA. " f"The error was: {e}"
            )
    os.remove(gfa_fp)


def _fetch_and_extract_grch38(get_ncbi_genomes: callable, dest_dir: str):
    """
    Fetches and extracts the GRCh38 human reference genome.

    This function uses the RESCRIPt method provided through the callable
    `get_ncbi_genomes` to fetch the GRCh38 human reference genome.
    The fetched 'dna-sequences.fasta' file is renamed to 'grch38.fasta'.

    Args:
        get_ncbi_genomes (callable): A function to fetch genomes from NCBI.
        dest_dir (str): The directory where the genome data will be saved.
    """
    results = get_ncbi_genomes(
        taxa=["Homo sapiens"],
        only_reference=True,
        assembly_levels=["chromosome"],
        assembly_source="refseq",
        only_genomic=False,
    )
    results.genome_assemblies.export_data(dest_dir)
    shutil.move(
        os.path.join(dest_dir, "dna-sequences.fasta"),
        os.path.join(dest_dir, "grch38.fasta"),
    )


def _combine_fasta_files(*fasta_in_fp, fasta_out_fp):
    """
    Combines multiple FASTA files into a single FASTA file.

    This function uses 'seqtk' to format and combine multiple FASTA files
    into a single file. Each input FASTA file is appended to the output file.
    After processing, the input files are removed.

    Args:
        *fasta_in_fp: Variable length argument list of paths to input
            FASTA files.
        fasta_out_fp (str): The file path where the combined output FASTA
            file should be saved.
    """
    with open(fasta_out_fp, "a") as f_out:
        for f_in in fasta_in_fp:
            try:
                subprocess.run(["seqtk", "seq", "-U", f_in], stdout=f_out)
            except Exception as e:
                raise Exception(
                    f"Failed to add the {f_in} to the reference FASTA file. "
                    f"The error was: {e}"
                )
            os.remove(f_in)


def construct_pangenome_index(ctx, threads=1):
    """Constructs the pangenome index.

    This action will fetch the human pangenome and GRCh38 reference genome,
    combine them into a single FASTA file, and generate a Bowtie 2 index.
    """
    get_ncbi_genomes = ctx.get_action("rescript", "get_ncbi_genomes")
    build_index = ctx.get_action("quality_control", "bowtie2_build")

    with tempfile.TemporaryDirectory() as tmp:
        print("Fetching the human pangenome GFA file...")
        _fetch_and_extract_pangenome(EBI_SERVER_URL, tmp)

        print("Fetching the human GRCh38 reference genome...")
        _fetch_and_extract_grch38(get_ncbi_genomes, tmp)

        print("Converting pangenome GFA to FASTA...")
        gfa_fp = glob.glob(os.path.join(tmp, "*.gfa"))[0]
        pan_fasta_fp = os.path.join(tmp, "pangenome.fasta")
        _extract_fasta_from_gfa(gfa_fp, pan_fasta_fp)

        print("Generating an index of the combined reference...")
        combined_fasta_fp = os.path.join(tmp, "combined.fasta")
        _combine_fasta_files(
            pan_fasta_fp,
            os.path.join(tmp, "grch38.fasta"),
            fasta_out_fp=combined_fasta_fp,
        )
        combined_reference = ctx.make_artifact(
            "FeatureData[Sequence]", combined_fasta_fp
        )
        (index,) = build_index(sequences=combined_reference, n_threads=threads)
    return index


def filter_reads_pangenome(
    ctx,
    reads,
    index=None,
    threads=1,
    mode="local",
    sensitivity="sensitive",
    ref_gap_open_penalty=5,
    ref_gap_ext_penalty=3,
):
    """
    Filters reads against a pangenome index, optionally generating the index
    if not provided.

    This function fetches and processes the human pangenome and GRCh38
    reference genome, combines them into a single FASTA file, and then
    generates a Bowtie 2 index if not already provided. It then filters
    reads against this index according to the specified sensitivity.
    """

    filter_reads = ctx.get_action("quality_control", "filter_reads")
    construct_index = ctx.get_action("annotate", "construct_pangenome_index")

    if index is None:
        print("Reference index was not provided - it will be generated.")
        (index,) = construct_index(threads)

    print("Filtering reads against the index...")
    filter_params = {
        k: v
        for k, v in locals().items()
        if k in ["mode", "ref_gap_open_penalty", "ref_gap_ext_penalty"]
    }
    (filtered_reads,) = filter_reads(
        demultiplexed_sequences=reads,
        database=index,
        exclude_seqs=True,
        n_threads=threads,
        **filter_params,
    )

    return filtered_reads, index
