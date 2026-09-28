#!/usr/bin/env python3

"""
Download pre-built Sma3s UniProt reference databases.

Each database requires four files:

    <prefix>.fasta
    <prefix>.annot
    <prefix>.fasta.sma3s_fasta_index.sqlite
    <prefix>.annot.q0_go0_goslim0.sma3s_annot.sqlite

Examples:
    python3 sma3s_db_download.py --list-dbs
    python3 sma3s_db_download.py --db bacteria
    python3 sma3s_db_download.py --db archaea
"""

import argparse
import shutil
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional
from urllib.error import HTTPError, URLError
from urllib.parse import urljoin
from urllib.request import Request, urlopen


DEFAULT_BASE_URL = "https://sma3s-db.clinbioinfosspa.es/"
DEFAULT_OUTPUT_DIR = "db_sma3s"
DEFAULT_DB = "bacteria"
DEFAULT_WORKERS = 2
DEFAULT_MIN_FREE_GB = 20.0
DEFAULT_RETRIES = 5
DEFAULT_TIMEOUT = 60
DEFAULT_CHUNK_MB = 8


DATABASES = {
    "bacteria": {
        "description": "UniProt bacterial proteins (Swiss-Prot + TrEMBL)",
        "prefix": "uniprot_bacteria",
    },
    "archaea": {
        "description": "UniProt archaeal proteins (Swiss-Prot + TrEMBL)",
        "prefix": "uniprot_archaea",
    },
    "virus": {
        "description": "UniProt viral proteins (Swiss-Prot + TrEMBL)",
        "prefix": "uniprot_virus",
    },
    "fungi": {
        "description": "UniProt fungal proteins (Swiss-Prot + TrEMBL)",
        "prefix": "uniprot_fungi",
    },
    "plants": {
        "description": "UniProt plant proteins (Swiss-Prot + TrEMBL)",
        "prefix": "uniprot_plants",
    },
    "metazoa": {
        "description": "UniProt metazoan proteins (Swiss-Prot + TrEMBL)",
        "prefix": "uniprot_metazoa",
    },
}


ALIASES = {
    "bacterial": "bacteria",
    "archaeal": "archaea",
    "viral": "virus",
    "viruses": "virus",
    "fungal": "fungi",
    "plant": "plants",
    "metazoan": "metazoa",
}


def required_files(prefix: str) -> list[str]:
    return [
        f"{prefix}.fasta",
        f"{prefix}.annot",
        f"{prefix}.fasta.sma3s_fasta_index.sqlite",
        f"{prefix}.annot.q0_go0_goslim0.sma3s_annot.sqlite",
    ]


def normalize_db_name(name: str) -> str:
    name = name.strip().lower()
    return ALIASES.get(name, name)


def human_size(n: int) -> str:
    if n >= 1024 ** 3:
        return f"{n / 1024 ** 3:.2f} GB"
    return f"{n / 1024 ** 2:.2f} MB"


def remote_info(url: str, timeout: int) -> tuple[bool, Optional[int]]:
    """
    Check whether a remote file exists using a minimal GET request.

    A Range request is used instead of HEAD because some web servers
    do not handle HEAD requests correctly for downloadable files.

    If the server supports byte ranges, only one byte is requested.
    If it ignores the Range header, urllib opens the response but the
    file body is not consumed here.
    """
    request = Request(
        url,
        headers={
            "User-Agent": "Sma3s-db-downloader/1.0",
            "Range": "bytes=0-0",
        },
    )

    try:
        with urlopen(request, timeout=timeout) as response:
            # For HTTP 206, Content-Range looks like:
            # bytes 0-0/123456789
            content_range = response.headers.get("Content-Range")

            if content_range and "/" in content_range:
                total = content_range.rsplit("/", 1)[1]
                if total.isdigit():
                    return True, int(total)

            # If the server ignores Range and returns HTTP 200,
            # Content-Length normally contains the complete file size.
            content_length = response.headers.get("Content-Length")
            size = int(content_length) if content_length else None

            return True, size

    except HTTPError as error:
        if error.code == 404:
            return False, None
        raise


def database_status(
    db_name: str,
    base_url: str,
    timeout: int,
) -> tuple[bool, list[str]]:
    missing = []
    prefix = DATABASES[db_name]["prefix"]

    for filename in required_files(prefix):
        url = urljoin(base_url, filename)

        try:
            exists, _ = remote_info(url, timeout)
        except (HTTPError, URLError, TimeoutError) as error:
            raise RuntimeError(
                f"Could not contact the database server while checking "
                f"{filename}: {error}"
            ) from error

        if not exists:
            missing.append(filename)

    return len(missing) == 0, missing


def list_databases(base_url: str, timeout: int) -> None:
    print("Sma3s reference databases")
    print("=========================")
    print()

    for name, info in DATABASES.items():
        try:
            ready, _ = database_status(name, base_url, timeout)
            status = "AVAILABLE" if ready else "UNDER PREPARATION"
        except RuntimeError:
            status = "SERVER CHECK FAILED"

        print(f"{name:<10} {status:<20} {info['description']}")


def check_free_space(path: Path, min_free_gb: float) -> None:
    free_gb = shutil.disk_usage(path).free / 1024 ** 3

    if free_gb < min_free_gb:
        raise RuntimeError(
            f"Not enough free disk space in {path.resolve()}. "
            f"Available: {free_gb:.2f} GB; "
            f"required minimum: {min_free_gb:.2f} GB."
        )


def download_file(
    filename: str,
    base_url: str,
    output_dir: Path,
    overwrite: bool,
    retries: int,
    timeout: int,
    chunk_size: int,
) -> str:

    url = urljoin(base_url, filename)
    output_path = output_dir / filename
    part_path = output_dir / f"{filename}.part"

    _, remote_size = remote_info(url, timeout)

    if output_path.exists() and not overwrite:
        if remote_size is None or output_path.stat().st_size == remote_size:
            return f"[SKIP] {filename} already exists"

        raise RuntimeError(
            f"{filename} already exists but its size differs from "
            f"the server file. Use --overwrite to download it again."
        )

    if overwrite:
        output_path.unlink(missing_ok=True)
        part_path.unlink(missing_ok=True)

    for attempt in range(1, retries + 1):
        existing_size = part_path.stat().st_size if part_path.exists() else 0

        headers = {"User-Agent": "Sma3s-db-downloader/1.0"}
        if existing_size:
            headers["Range"] = f"bytes={existing_size}-"

        request = Request(url, headers=headers)

        # DOWNLOAD COMMAND LOCATION:
        # If wget is preferred in the future, replace the urllib download
        # block below with the wget command here.
        # Example: wget -c -O <filename>.part <url>

        try:
            with urlopen(request, timeout=timeout) as response:
                if existing_size and response.getcode() != 206:
                    part_path.unlink(missing_ok=True)
                    existing_size = 0

                mode = "ab" if existing_size else "wb"
                downloaded = existing_size
                last_report = time.time()

                with open(part_path, mode) as handle:
                    while True:
                        chunk = response.read(chunk_size)
                        if not chunk:
                            break

                        handle.write(chunk)
                        downloaded += len(chunk)

                        if time.time() - last_report >= 10:
                            if remote_size:
                                percent = 100 * downloaded / remote_size
                                print(
                                    f"[PROGRESS] {filename}: "
                                    f"{human_size(downloaded)} / "
                                    f"{human_size(remote_size)} "
                                    f"({percent:.1f}%)"
                                )
                            else:
                                print(
                                    f"[PROGRESS] {filename}: "
                                    f"{human_size(downloaded)}"
                                )
                            last_report = time.time()

            if remote_size is not None:
                final_size = part_path.stat().st_size
                if final_size != remote_size:
                    raise RuntimeError(
                        f"Incomplete download: "
                        f"{human_size(final_size)} / "
                        f"{human_size(remote_size)}"
                    )

            part_path.replace(output_path)

            return (
                f"[OK] {filename} downloaded"
                + (
                    f" ({human_size(remote_size)})"
                    if remote_size is not None
                    else ""
                )
            )

        except (HTTPError, URLError, TimeoutError, RuntimeError) as error:
            if attempt == retries:
                raise RuntimeError(
                    f"{filename} failed after {retries} attempts: {error}"
                ) from error

            wait = min(60, 5 * attempt)
            print(
                f"[RETRY] {filename}: attempt {attempt}/{retries} failed. "
                f"Retrying in {wait}s..."
            )
            time.sleep(wait)

    raise RuntimeError(f"Could not download {filename}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download a pre-built Sma3s UniProt reference database."
    )

    parser.add_argument(
        "--db",
        default=DEFAULT_DB,
        metavar="DATABASE",
        help=(
            f"Database to download. Default: {DEFAULT_DB}. "
            "Use --list-dbs to check availability."
        ),
    )

    parser.add_argument(
        "--list-dbs",
        action="store_true",
        help="List databases and their current availability.",
    )

    parser.add_argument(
        "--url",
        default=DEFAULT_BASE_URL,
        help=f"Database server URL. Default: {DEFAULT_BASE_URL}",
    )

    parser.add_argument(
        "-o",
        "--output",
        default=DEFAULT_OUTPUT_DIR,
        help=(
            "Root download directory. "
            f"Default: {DEFAULT_OUTPUT_DIR}"
        ),
    )

    parser.add_argument(
        "-w",
        "--workers",
        type=int,
        default=DEFAULT_WORKERS,
        help=(
            "Number of simultaneous downloads. "
            f"Default: {DEFAULT_WORKERS}"
        ),
    )

    parser.add_argument(
        "--min-free-gb",
        type=float,
        default=DEFAULT_MIN_FREE_GB,
        help=(
            "Minimum required free disk space in GB. "
            f"Default: {DEFAULT_MIN_FREE_GB}"
        ),
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing files.",
    )

    parser.add_argument(
        "--retries",
        type=int,
        default=DEFAULT_RETRIES,
        help=f"Retries per file. Default: {DEFAULT_RETRIES}",
    )

    parser.add_argument(
        "--timeout",
        type=int,
        default=DEFAULT_TIMEOUT,
        help=f"Connection timeout in seconds. Default: {DEFAULT_TIMEOUT}",
    )

    parser.add_argument(
        "--chunk-mb",
        type=int,
        default=DEFAULT_CHUNK_MB,
        help=f"Download chunk size in MB. Default: {DEFAULT_CHUNK_MB}",
    )

    args = parser.parse_args()

    base_url = args.url.rstrip("/") + "/"

    if args.list_dbs:
        list_databases(base_url, args.timeout)
        return

    db_name = normalize_db_name(args.db)

    if db_name not in DATABASES:
        print(f"[ERROR] Unknown database: {args.db}")
        print("Supported databases:")
        for name in DATABASES:
            print(f"  - {name}")
        sys.exit(2)

    print(f"Checking '{db_name}' database on the server...")

    try:
        ready, missing_files = database_status(
            db_name,
            base_url,
            args.timeout,
        )
    except RuntimeError as error:
        print(f"[ERROR] {error}")
        sys.exit(1)

    if not ready:
        print()
        print(
            f"[UNDER PREPARATION] The '{db_name}' database "
            "is not ready for download."
        )
        print()
        print(
            "Sma3s requires all four database files. "
            "The following file(s) are not currently available:"
        )
        for filename in missing_files:
            print(f"  - {filename}")
        sys.exit(2)

    prefix = DATABASES[db_name]["prefix"]
    files = required_files(prefix)

    # Default layout:
    # db_sma3s/bacteria/
    # db_sma3s/archaea/
    output_dir = Path(args.output) / db_name
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        check_free_space(output_dir, args.min_free_gb)
    except RuntimeError as error:
        print(f"[ERROR] {error}")
        sys.exit(1)

    workers = max(1, min(args.workers, len(files)))
    chunk_size = max(1, args.chunk_mb) * 1024 * 1024

    print()
    print("Sma3s database downloader")
    print("========================")
    print(f"Database: {db_name}")
    print(f"Server: {base_url}")
    print(f"Output directory: {output_dir.resolve()}")
    print(f"Simultaneous downloads: {workers}")
    print()
    print("Files:")
    for filename in files:
        print(f"  - {filename}")
    print()

    failed = []

    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                download_file,
                filename,
                base_url,
                output_dir,
                args.overwrite,
                args.retries,
                args.timeout,
                chunk_size,
            ): filename
            for filename in files
        }

        for future in as_completed(futures):
            filename = futures[future]
            try:
                print(future.result())
            except Exception as error:
                print(f"[ERROR] {filename}: {error}")
                failed.append(filename)

    print()

    if failed:
        print("Some files could not be downloaded:")
        for filename in failed:
            print(f"  - {filename}")
        print("Run the same command again to resume partial downloads.")
        sys.exit(1)

    print(f"Sma3s '{db_name}' database downloaded successfully.")
    print(f"Database directory: {output_dir.resolve()}")


if __name__ == "__main__":
    main()
