# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import tarfile
import tempfile
import shutil
import subprocess
from urllib.request import urlretrieve

import requests
from bs4 import BeautifulSoup
from tqdm import tqdm
import pandas as pd

from q2_types.kaiju import KaijuDBDirectoryFormat
from q2_types.genome_data import ProteinsDirectoryFormat
from qiime2 import Metadata
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


def _fetch_ncbi_taxonomy_files(target_dir: str):
    """
    Fetch NCBI taxonomy files (nodes.dmp and names.dmp) needed for Kaiju.
    
    Args:
        target_dir (str): Directory to save the taxonomy files
    """
    # NCBI taxonomy FTP URL
    base_url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/"
    taxonomy_url = f"{base_url}taxdump.tar.gz"
    
    print("Downloading NCBI taxonomy files...")
    
    # Download the taxonomy archive
    response = requests.get(taxonomy_url, stream=True)
    response.raise_for_status()
    
    total_size = int(response.headers.get("content-length", 0))
    progress_bar = None
    if total_size > 0:
        progress_bar = tqdm(
            desc="Downloading NCBI taxonomy",
            total=total_size,
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
        )
    
    # Save to temporary file
    with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
        for chunk in response.iter_content(chunk_size=CHUNK_SIZE):
            if chunk:
                tmp_file.write(chunk)
                if progress_bar:
                    progress_bar.update(len(chunk))
        tmp_filepath = tmp_file.name
    
    if progress_bar:
        progress_bar.close()
    
    print("Extracting taxonomy files...")
    
    # Extract only the files we need (nodes.dmp and names.dmp)
    with tarfile.open(tmp_filepath, "r:gz") as tar:
        for member in tar.getmembers():
            if member.name in ["nodes.dmp", "names.dmp"]:
                tar.extract(member, target_dir)
    
    # Clean up temporary file
    os.remove(tmp_filepath)
    print("NCBI taxonomy files downloaded successfully.")


def _create_protein_fasta(proteins_dir: ProteinsDirectoryFormat, metadata: Metadata, output_path: str):
    """
    Create a single FASTA file from protein files with proper headers containing NCBI taxon IDs.
    
    Args:
        proteins_dir: Directory containing protein FASTA files
        metadata: Metadata mapping genome IDs to NCBI taxon IDs
        output_path: Path to output combined FASTA file
    """
    metadata_df = metadata.to_dataframe()
    
    # Validate that metadata contains the required column
    if metadata_df.empty:
        raise ValueError("Metadata is empty")
    
    # Assume the taxon ID is in a column named 'taxon_id' or similar
    # We'll look for columns that might contain taxon IDs
    taxon_columns = [col for col in metadata_df.columns if 'taxon' in col.lower() or 'tax' in col.lower()]
    if not taxon_columns:
        # If no obvious taxon column, use the first column
        if len(metadata_df.columns) == 0:
            raise ValueError("Metadata must contain at least one column with taxon IDs")
        taxon_column = metadata_df.columns[0]
    else:
        taxon_column = taxon_columns[0]
    
    print(f"Using '{taxon_column}' column for taxon IDs")
    
    # Get all protein files
    protein_files = {}
    if hasattr(proteins_dir, 'feature_dict'):
        protein_files = proteins_dir.feature_dict()
    else:
        # If feature_dict is not available, look for .fasta files
        for filename in os.listdir(proteins_dir.path):
            if filename.endswith('.fasta') or filename.endswith('.fa'):
                genome_id = os.path.splitext(filename)[0]
                protein_files[genome_id] = os.path.join(proteins_dir.path, filename)
    
    if not protein_files:
        raise ValueError("No protein files found in the input directory")
    
    # Create combined FASTA file
    with open(output_path, 'w') as outfile:
        for genome_id, protein_file in protein_files.items():
            if genome_id not in metadata_df.index:
                print(f"Warning: Genome ID '{genome_id}' not found in metadata, skipping")
                continue
            
            # Get taxon ID for this genome
            taxon_id = metadata_df.loc[genome_id, taxon_column]
            
            # Read protein sequences and modify headers
            with open(protein_file, 'r') as infile:
                for line in infile:
                    if line.startswith('>'):
                        # Modify header to include taxon ID
                        protein_id = line[1:].strip()
                        new_header = f">{protein_id}_{genome_id}_taxid_{taxon_id}\n"
                        outfile.write(new_header)
                    else:
                        outfile.write(line)
    
    print(f"Created combined protein FASTA file with {len(protein_files)} genomes")


def build_kaiju_db(
    proteins: ProteinsDirectoryFormat,
    metadata: Metadata,
) -> KaijuDBDirectoryFormat:
    """
    Build a custom Kaiju database from protein sequences and metadata.
    
    Args:
        proteins: Directory containing protein FASTA files for each genome
        metadata: Metadata mapping genome IDs to NCBI taxon IDs
        
    Returns:
        KaijuDBDirectoryFormat: The built Kaiju database
    """
    # Create output directory
    db = KaijuDBDirectoryFormat()
    
    with tempfile.TemporaryDirectory() as tmp_dir:
        # Step 1: Fetch NCBI taxonomy files
        _fetch_ncbi_taxonomy_files(tmp_dir)
        
        # Step 2: Create combined protein FASTA file
        protein_fasta_path = os.path.join(tmp_dir, "proteins.faa")
        _create_protein_fasta(proteins, metadata, protein_fasta_path)
        
        # Step 3: Run kaiju-mkbwt
        print("Running kaiju-mkbwt to create BWT index...")
        bwt_output = os.path.join(tmp_dir, "proteins.bwt")
        sa_output = os.path.join(tmp_dir, "proteins.sa")
        
        mkbwt_cmd = [
            "kaiju-mkbwt",
            "-a", "ACDEFGHIKLMNPQRSTVWY",  # Protein alphabet
            "-o", bwt_output,
            protein_fasta_path
        ]
        
        try:
            run_command(cmd=mkbwt_cmd, verbose=True)
        except subprocess.CalledProcessError as e:
            raise Exception(
                f"An error was encountered while running kaiju-mkbwt "
                f"(return code {e.returncode}), please inspect "
                "stdout and stderr to learn more."
            )
        
        # Step 4: Run kaiju-mkfmi
        print("Running kaiju-mkfmi to create FMI index...")
        fmi_output = os.path.join(tmp_dir, "proteins.fmi")
        
        mkfmi_cmd = [
            "kaiju-mkfmi",
            bwt_output,
            sa_output,
            fmi_output
        ]
        
        try:
            run_command(cmd=mkfmi_cmd, verbose=True)
        except subprocess.CalledProcessError as e:
            raise Exception(
                f"An error was encountered while running kaiju-mkfmi "
                f"(return code {e.returncode}), please inspect "
                "stdout and stderr to learn more."
            )
        
        # Step 5: Copy required files to final database directory
        # Copy the FMI index
        shutil.copy2(fmi_output, os.path.join(db.path, "kaiju_db.fmi"))
        
        # Copy taxonomy files
        shutil.copy2(os.path.join(tmp_dir, "nodes.dmp"), db.path)
        shutil.copy2(os.path.join(tmp_dir, "names.dmp"), db.path)
        
        print("Kaiju database built successfully!")
    
    return db
