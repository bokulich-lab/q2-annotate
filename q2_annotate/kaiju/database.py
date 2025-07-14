# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import tarfile

import requests
from bs4 import BeautifulSoup
from tqdm import tqdm

from q2_types.kaiju import KaijuDBDirectoryFormat

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
