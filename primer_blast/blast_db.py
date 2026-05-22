"""BLAST database management — auto-build from genome FASTA with sentinel caching."""

import hashlib
import os
import subprocess
from typing import Optional


def _fasta_checksum(fasta_path: str) -> str:
    """Compute a quick checksum of the FASTA file (first + last bytes + size)."""
    size = os.path.getsize(fasta_path)
    with open(fasta_path, "rb") as f:
        head = f.read(4096)
        f.seek(max(0, size - 4096))
        tail = f.read(4096)
    return hashlib.md5(head + tail + str(size).encode()).hexdigest()


def is_blast_db_current(fasta_path: str, db_dir: str) -> bool:
    """Check if BLAST database is up-to-date with the FASTA file."""
    sentinel = os.path.join(db_dir, ".blastdb_built")
    if not os.path.exists(sentinel):
        return False

    with open(sentinel, "r") as f:
        stored_checksum = f.read().strip()

    current_checksum = _fasta_checksum(fasta_path)
    return stored_checksum == current_checksum


def build_blast_db(
    fasta_path: str,
    db_dir: str,
    db_name: str = "gmax_genome",
    makeblastdb_path: str = "makeblastdb",
) -> str:
    """Build a BLAST nucleotide database from the genome FASTA.

    Args:
        fasta_path: Path to genome FASTA.
        db_dir: Directory to store the BLAST database files.
        db_name: Base name for the database files.
        makeblastdb_path: Path to makeblastdb executable.

    Returns:
        Path prefix of the built database.

    Raises:
        RuntimeError: If makeblastdb fails.
    """
    os.makedirs(db_dir, exist_ok=True)
    db_prefix = os.path.join(db_dir, db_name)

    cmd = [
        makeblastdb_path,
        "-dbtype", "nucl",
        "-in", fasta_path,
        "-out", db_prefix,
        "-parse_seqids",
        "-title", "Glycine max Wm82 v4.0",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

    if result.returncode != 0:
        raise RuntimeError(f"makeblastdb failed: {result.stderr}")

    # Write sentinel
    sentinel = os.path.join(db_dir, ".blastdb_built")
    with open(sentinel, "w") as f:
        f.write(_fasta_checksum(fasta_path))

    return db_prefix


def get_blast_db(
    fasta_path: str,
    db_dir: str = "blastdb",
    db_name: str = "gmax_genome",
    makeblastdb_path: str = "makeblastdb",
    force_rebuild: bool = False,
) -> str:
    """Get BLAST database path, building it if needed.

    Args:
        fasta_path: Path to genome FASTA.
        db_dir: Directory to store the BLAST database.
        db_name: Base name for the database files.
        makeblastdb_path: Path to makeblastdb executable.
        force_rebuild: Force database rebuild even if current.

    Returns:
        Path prefix of the BLAST database.
    """
    if force_rebuild or not is_blast_db_current(fasta_path, db_dir):
        print(f"Building BLAST database from {fasta_path}...", flush=True)
        print("This may take several minutes for a ~1 GB genome.", flush=True)
        return build_blast_db(
            fasta_path=fasta_path,
            db_dir=db_dir,
            db_name=db_name,
            makeblastdb_path=makeblastdb_path,
        )

    return os.path.join(db_dir, db_name)
