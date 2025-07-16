# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import tempfile
import subprocess
import tarfile
from typing import List

import requests
from bs4 import BeautifulSoup
from tqdm import tqdm

from q2_types.kaiju import KaijuDBDirectoryFormat
from q2_types.feature_data import ProteinFASTAFormat
from q2_types.reference_db import ReferenceDB, NCBITaxonomy
from q2_annotate._utils import run_command

CHUNK_SIZE = 8192
KAIJU_SERVER_URL = "https://bioinformatics-centre.github.io/" "kaiju/downloads.html"
ERR_MSG = (
    "Unable to connect to the Kaiju server. Please try again later. "
    "The error was: {}"
)


def _fetch_and_extract_db(db_uri: str, db_dir: str):
    """
    Fetches and extracts the Kaiju database.

    Args:
        db_uri (str): The URI of the database to fetch.
        db_dir (str): Path to the final DB directory.
    """
    latest_db = os.path.basename(db_uri)
    db_path = os.path.join(db_dir, latest_db)
    try:
        response = requests.get(db_uri, stream=True)
        response.raise_for_status()
        total_size = int(response.headers.get("content-length", 0))
        if total_size > 0:
            progress_bar = tqdm(
                desc=f'Downloading the "{latest_db}" database',
                total=total_size,
                unit="B",
                unit_scale=True,
                unit_divisor=1024,
            )

        with open(db_path, "wb") as file:
            for chunk in response.iter_content(chunk_size=CHUNK_SIZE):
                file.write(chunk) if chunk else False
                if total_size > 0:
                    progress_bar.update(len(chunk))
            progress_bar.close() if total_size > 0 else False
    except requests.exceptions.ConnectionError as e:
        raise Exception(ERR_MSG.format(e))

    msg = "Download finished. Extracting database files..."
    print(f"{msg}", end="", flush=True)
    with tarfile.open(db_path, "r:gz") as tar:
        tar.extractall(path=db_dir)
    print(f"\r{msg} Done.", flush=True)

    os.remove(db_path)


def _find_latest_db_url(response: bytes, database_type: str) -> str:
    """
    Finds the latest database URL based on the database type.

    Args:
        response (bytes): HTML response containing the table with DB URLs.
        database_type (str): The target database type to filter.

    Returns:
        str: The latest database URL.
    """
    soup = BeautifulSoup(response, "html.parser")
    tables = soup.find_all("table")

    for table in tables:
        # Locate the table header
        headers = table.find_all("th")
        if headers and headers[0].get_text().strip() == "Database":
            rows = table.find_all("tr")
            for row in rows:
                cells = row.find_all("td")

                # Check if the first cell contains the required database_type
                if cells and cells[0].get_text().strip() == database_type:
                    # The next row contains the desired URLs
                    next_row = row.find_next_sibling("tr")
                    if next_row:
                        url_cell = next_row.find_all("td")[-1]
                        url = url_cell.find("a")
                        if url:
                            return url["href"]

    raise ValueError(f"URL for database type '{database_type}' not found.")


def _build_kaiju_bwt_index(fasta_file: str, output_dir: str, threads: int) -> str:
    """
    Build BWT index using kaiju-mkbwt.
    
    Args:
        fasta_file: Path to input FASTA file
        output_dir: Directory where index files will be created
        threads: Number of threads to use
        
    Returns:
        str: Path to the generated BWT file
    """
    bwt_file = os.path.join(output_dir, "kaiju_db_idx.bwt")
    
    cmd = [
        "kaiju-mkbwt",
        "-n", str(threads),
        "-a", "ACDEFGHIKLMNPQRSTVWY",
        "-o", bwt_file,
        fasta_file
    ]
    
    try:
        run_command(cmd=cmd, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            f"An error was encountered while building the BWT index, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )
    
    return bwt_file


def _build_kaiju_fmi_index(bwt_file: str) -> str:
    """
    Build FM index using kaiju-mkfmi.
    
    Args:
        bwt_file: Path to BWT file
        
    Returns:
        str: Path to the generated FMI file
    """
    cmd = [
        "kaiju-mkfmi",
        bwt_file
    ]
    
    try:
        run_command(cmd=cmd, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            f"An error was encountered while building the FM index, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )
    
    # The FMI file is created with the same base name as the BWT file
    fmi_file = bwt_file.replace('.bwt', '.fmi')
    return fmi_file


def build_kaiju_db(
    seqs: List[ProteinFASTAFormat],
    taxonomy: ReferenceDB[NCBITaxonomy],
    threads: int = 1,
) -> KaijuDBDirectoryFormat:
    """
    Build custom Kaiju database from protein sequences.
    
    Args:
        seqs: List of protein FASTA sequences to include in the database
        taxonomy: NCBI taxonomy reference database containing nodes.dmp and names.dmp
        threads: Number of threads to use for index construction
        
    Returns:
        KaijuDBDirectoryFormat: Built Kaiju database
    """
    if not seqs:
        raise ValueError("At least one sequence file must be provided.")
    
    db = KaijuDBDirectoryFormat()
    
    with tempfile.TemporaryDirectory() as tmp_dir:
        # Combine all FASTA files into a single file
        combined_fasta = os.path.join(tmp_dir, "combined_seqs.faa")
        
        with open(combined_fasta, 'w') as outfile:
            for seq_file in seqs:
                with open(str(seq_file.path), 'r') as infile:
                    content = infile.read()
                    outfile.write(content)
                    # Ensure there's a newline at the end if missing
                    if content and not content.endswith('\n'):
                        outfile.write('\n')
        
        # Build BWT index
        bwt_file = _build_kaiju_bwt_index(combined_fasta, tmp_dir, threads)
        
        # Build FM index  
        fmi_file = _build_kaiju_fmi_index(bwt_file)
        
        # Copy FMI file to final destination
        dst_fmi_path = os.path.join(str(db.path), "kaiju_db_idx.fmi")
        if os.path.exists(fmi_file):
            os.rename(fmi_file, dst_fmi_path)
        
        # Copy taxonomy files to final destination
        import shutil
        nodes_src = os.path.join(str(taxonomy.path), "nodes.dmp")
        names_src = os.path.join(str(taxonomy.path), "names.dmp")
        nodes_dst = os.path.join(str(db.path), "nodes.dmp")
        names_dst = os.path.join(str(db.path), "names.dmp")
        
        if os.path.exists(nodes_src):
            shutil.copy2(nodes_src, nodes_dst)
        if os.path.exists(names_src):
            shutil.copy2(names_src, names_dst)
    
    return db


def fetch_kaiju_db(
    database_type: str,
) -> KaijuDBDirectoryFormat:

    try:
        response = requests.get(KAIJU_SERVER_URL)
    except requests.exceptions.RequestException as e:
        raise Exception(ERR_MSG.format(e))

    download_link = _find_latest_db_url(response.content, database_type)

    db = KaijuDBDirectoryFormat()
    _fetch_and_extract_db(download_link, str(db.path))

    return db
