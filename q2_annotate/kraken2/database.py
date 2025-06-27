# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import os
import re
import shutil
import subprocess
import tarfile
import tempfile
from copy import deepcopy
from typing import List

import requests
import xmltodict
from q2_types.feature_data import DNAFASTAFormat
from tqdm import tqdm

from q2_annotate._utils import _process_common_input_params, run_command
from q2_types.kraken2 import (
    Kraken2DBDirectoryFormat,
    BrackenDBDirectoryFormat,
    Kraken2DBReportDirectoryFormat,
)

from q2_annotate.kraken2.utils import _process_kraken2_arg

COLLECTIONS = {
    "standard": "standard",
    "viral": "viral",
    "minusb": "minusb",
    "standard8": "standard_08gb",
    "standard16": "standard_16gb",
    "pluspf": "pluspf",
    "pluspf8": "pluspf_08gb",
    "pluspf16": "pluspf_16gb",
    "pluspfp": "pluspfp",
    "pluspfp8": "pluspfp_08gb",
    "pluspfp16": "pluspfp_16gb",
    "eupathdb": "eupathdb48",
    "nt": "nt",
    "corent": "core_nt",
    "gtdb": "gtdb_genome_reps",
    "greengenes": "Greengenes13.5",
    "rdp": "RDP11.5",
    "silva132": "Silva132",
    "silva138": "Silva138",
}
S16_DBS = ["greengenes", "rdp", "silva132", "silva138"]
S3_COLLECTIONS_URL = "https://genome-idx.s3.amazonaws.com"
CHUNK_SIZE = 8192


def _fetch_taxonomy(db_dir: str, threads: int, use_ftp: bool):
    cmd = [
        "kraken2-build",
        "--download-taxonomy",
        "--threads",
        str(threads),
        "--db",
        str(db_dir),
    ]
    cmd.append("--use-ftp") if use_ftp else False
    try:
        run_command(cmd=cmd, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while downloading taxonomy, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )


def _fetch_libraries(db_dir: str, libraries: List[str], all_kwargs: dict):
    kwargs = {
        k: v for k, v in all_kwargs.items() if k in ["threads", "no_masking", "use_ftp"]
    }
    common_args = _process_common_input_params(
        processing_func=_process_kraken2_arg, params=kwargs
    )
    common_args.extend(["--db", db_dir])
    base_cmd = [
        "kraken2-build",
        "--download-library",
    ]
    fetch_behaviour = all_kwargs.get("library_exists", "refetch")
    for library in libraries:
        if fetch_behaviour == "skip":
            lib_path = os.path.join(db_dir, "library", library)
            if os.path.exists(lib_path) and len(os.listdir(lib_path)) > 0:
                print(
                    f"Skipping download of the '{library}' library, " f"already exists."
                )
                continue
        try:
            cmd = deepcopy(base_cmd)
            cmd.extend([library, *common_args])
            run_command(cmd=cmd, verbose=True)
        except subprocess.CalledProcessError as e:
            raise Exception(
                f"An error was encountered while downloading the "
                f"'{library}' library, (return code {e.returncode}), "
                "please inspect stdout and stderr to learn more."
            )


def _add_seqs_to_library(db_dir: str, seqs: DNAFASTAFormat, no_masking: bool):
    cmd = ["kraken2-build", "--add-to-library", str(seqs.path), "--db", db_dir]
    cmd.append("--no-mask") if no_masking else False
    try:
        run_command(cmd=cmd, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while adding sequences to the "
            f"library, (return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )


def _build_kraken2_database(db_dir: str, all_kwargs: dict):
    kwargs = {
        k: v
        for k, v in all_kwargs.items()
        if k
        in [
            "threads",
            "minimizer_len",
            "minimizer_spaces",
            "load_factor",
            "fast_build",
            "kmer_len",
        ]
    }
    common_args = _process_common_input_params(
        processing_func=_process_kraken2_arg, params=kwargs
    )
    if all_kwargs["max_db_size"] > 0:
        common_args.extend(["--max-db-size", str(all_kwargs["max_db_size"])])
    cmd = ["kraken2-build", "--build", "--db", db_dir, *common_args]
    try:
        run_command(cmd=cmd, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while building the database, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )


def _build_bracken_database(
    kraken2_db_dir: str, threads: int, kmer_len: int, read_len: int
):
    cmd = [
        "bracken-build",
        "-d",
        kraken2_db_dir,
        "-t",
        str(threads),
        "-k",
        str(kmer_len),
        "-l",
        str(read_len),
    ]
    try:
        run_command(cmd=cmd, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while building the Bracken "
            f"database, (return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )


def _find_latest_db(collection: str, response: requests.Response) -> str:
    collection_id = COLLECTIONS[collection]

    if collection in S16_DBS:
        pattern = rf"kraken\/16S_{collection_id}_\d{{8}}.tgz"
    else:
        pattern = rf"kraken\/k2_{collection_id}_\d{{8}}.tar.gz"

    s3_objects = xmltodict.parse(response.content)
    s3_objects = s3_objects.get("ListBucketResult")
    if not s3_objects:
        raise ValueError(
            "No databases were found in the response returned by S3. "
            "Please try again."
        )
    s3_objects = [
        obj for obj in s3_objects["Contents"] if re.match(pattern, obj["Key"])
    ]
    s3_objects = sorted(s3_objects, key=lambda x: x["LastModified"], reverse=True)

    latest_db = s3_objects[0]["Key"]
    return latest_db


def _move_files_one_level_up(tmp_dir):
    extr_db = [x for x in os.listdir(tmp_dir) if not x.endswith(".tgz")][0]
    src = os.path.join(tmp_dir, extr_db)
    for file in os.listdir(src):
        shutil.move(os.path.join(src, file), tmp_dir)


def _fetch_db_collection(collection: str, tmp_dir: str):
    err_msg = (
        "Could not connect to the server. Please check your internet "
        "connection and try again. The error was: {}."
    )
    try:
        response = requests.get(S3_COLLECTIONS_URL)
    except requests.exceptions.ConnectionError as e:
        raise ValueError(err_msg.format(e))

    if response.status_code == 200:
        latest_db = _find_latest_db(collection, response)
        print(f'Found the latest "{collection}" database: {latest_db}.')
    else:
        raise ValueError(
            "Could not fetch the list of available databases. "
            f"Status code was: {response.status_code}. "
            "Please try again later."
        )

    db_uri = f"{S3_COLLECTIONS_URL}/{latest_db}"
    try:
        response = requests.get(db_uri, stream=True)
        total_size = int(response.headers.get("content-length", 0))
        if total_size > 0:
            progress_bar = tqdm(
                desc=f'Downloading the "{latest_db}" database',
                total=total_size,
                unit="B",
                unit_scale=True,
                unit_divisor=1024,
            )
        db_path = os.path.join(tmp_dir, os.path.split(db_uri)[-1])
        with open(db_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=CHUNK_SIZE):
                f.write(chunk) if chunk else False
                if total_size > 0:
                    progress_bar.update(len(chunk))
            progress_bar.close() if total_size > 0 else False
    except requests.exceptions.ConnectionError as e:
        raise ValueError(err_msg.format(e))

    msg = "Download finished. Extracting database files..."
    print(f"{msg}", end="", flush=True)
    with tarfile.open(db_path, "r:gz") as tar:
        tar.extractall(path=tmp_dir)
    print(f"\r{msg} Done.", flush=True)

    if collection in S16_DBS:
        _move_files_one_level_up(tmp_dir)


def _move_db_files(source: str, destination: str, extension: str = "k2d"):
    files = glob.glob(f"{source}/*.{extension}")
    for file in files:
        new_file = os.path.join(destination, os.path.split(file)[-1])
        shutil.move(file, new_file)


def _fetch_prebuilt_dbs(bracken_db, kraken2_db, collection, tmp):
    # Find files with the latest version
    _fetch_db_collection(collection=collection, tmp_dir=tmp)
    # Move the Kraken2/Bracken database files to the final location
    _move_db_files(tmp, str(kraken2_db.path), extension="k2d")
    _move_db_files(tmp, str(bracken_db.path), extension="kmer_distrib")


def _build_dbs_from_seqs(bracken_db, kraken2_db, seqs, tmp_dir, common_args):
    # Fetch taxonomy (also needed for custom databases)
    _fetch_taxonomy(
        db_dir=tmp_dir, threads=common_args["threads"], use_ftp=common_args["use_ftp"]
    )
    for seq in seqs:
        _add_seqs_to_library(
            db_dir=tmp_dir, seqs=seq, no_masking=common_args["no_masking"]
        )
    # Build the Kraken2 database
    _build_kraken2_database(db_dir=tmp_dir, all_kwargs=common_args)
    # Build the Bracken database
    for rl in common_args["read_len"]:
        _build_bracken_database(
            kraken2_db_dir=tmp_dir,
            threads=common_args["threads"],
            kmer_len=common_args["kmer_len"],
            read_len=rl,
        )
    # Move the Kraken2/Bracken database files to the final location
    _move_db_files(tmp_dir, str(kraken2_db.path), extension="k2d")
    _move_db_files(tmp_dir, str(bracken_db.path), extension="kmer_distrib")


def build_kraken_db(
    seqs: DNAFASTAFormat = None,
    collection: str = None,
    threads: int = 1,
    kmer_len: int = 35,
    minimizer_len: int = 31,
    minimizer_spaces: int = 7,
    no_masking: bool = False,
    max_db_size: int = 0,
    use_ftp: bool = False,
    load_factor: float = 0.7,
    fast_build: bool = False,
    read_len: int = None,
) -> (Kraken2DBDirectoryFormat, BrackenDBDirectoryFormat):
    kraken2_db = Kraken2DBDirectoryFormat()
    bracken_db = BrackenDBDirectoryFormat()

    if not read_len:
        # use the same values as in the pre-built databases
        read_len = [50, 75, 100, 150, 200, 250, 300]

    with tempfile.TemporaryDirectory() as tmp:
        if seqs:
            # Construct the custom-made database
            common_args = {
                k: v for k, v in locals().items() if k not in ["seqs", "collection"]
            }

            # Fetch taxonomy (also needed for custom databases)
            _build_dbs_from_seqs(bracken_db, kraken2_db, seqs, tmp, common_args)
        elif collection:
            _fetch_prebuilt_dbs(bracken_db, kraken2_db, collection, tmp)
        else:
            raise ValueError(
                "You need to either provide a list of sequences to build the "
                "database from or a valid collection name to be fetched from "
                '"Kraken 2/Bracken Refseq indexes" resource.'
            )

    return kraken2_db, bracken_db


def inspect_kraken2_db(
    db: Kraken2DBDirectoryFormat, threads: int = 1
) -> Kraken2DBReportDirectoryFormat:
    cmd = ["kraken2-inspect", "--db", str(db.path), "--threads", str(threads)]
    try:
        result = run_command(cmd=cmd, pipe=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while building the database, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )

    report_fmt = Kraken2DBReportDirectoryFormat()
    report_fp = os.path.join(report_fmt.path, report_fmt.file.pathspec)
    with open(report_fp, "w") as fh:
        fh.write(result.stdout)

    return report_fmt
