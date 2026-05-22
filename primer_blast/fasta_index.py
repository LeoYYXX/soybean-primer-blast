"""FASTA index for O(1) random access to genome sequences.

Builds and loads a .fai-style index mapping seqid -> (offset, line_bases, line_bytes).
"""

import os
from typing import Dict, Tuple


FAIEntry = Tuple[int, int, int, int]  # (byte_offset, seq_length, line_bases, line_bytes)


def build_fasta_index(fasta_path: str) -> Dict[str, FAIEntry]:
    """Scan FASTA and build seqid -> (offset, length, line_bases, line_bytes) dict.

    Also writes a .fai file for reuse.
    """
    fai_path = fasta_path + ".fai"
    index: Dict[str, FAIEntry] = {}

    # Try loading existing .fai first
    if os.path.exists(fai_path):
        fai_mtime = os.path.getmtime(fai_path)
        fasta_mtime = os.path.getmtime(fasta_path)
        if fai_mtime >= fasta_mtime:
            return _load_fai(fai_path)

    # Build from scratch
    with open(fasta_path, "rb") as f:
        offset = 0
        current_seqid = ""
        seq_length = 0
        line_bases = 0
        line_bytes = 0
        started = False

        while True:
            line = f.readline()
            if not line:
                if started and current_seqid:
                    index[current_seqid] = (offset, seq_length, line_bases, line_bytes)
                break

            line_str = line.decode("ascii", errors="ignore").rstrip("\r\n")

            if line.startswith(b">"):
                if started and current_seqid:
                    index[current_seqid] = (offset, seq_length, line_bases, line_bytes)

                current_seqid = line_str[1:].split()[0]
                offset = f.tell()
                seq_length = 0
                line_bases = 0
                line_bytes = 0
                started = True
            elif started and line_str:
                if line_bases == 0:
                    line_bases = len(line_str)
                    line_bytes = len(line)  # includes \n
                seq_length += len(line_str)

    # Write .fai
    with open(fai_path, "w") as out:
        for seqid, (off, slen, lbases, lbytes) in index.items():
            out.write(f"{seqid}\t{slen}\t{off}\t{lbases}\t{lbytes}\n")

    return index


def _load_fai(fai_path: str) -> Dict[str, FAIEntry]:
    """Load an existing .fai index file."""
    index: Dict[str, FAIEntry] = {}
    with open(fai_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 5:
                seqid = parts[0]
                slen = int(parts[1])
                offset = int(parts[2])
                lbases = int(parts[3])
                lbytes = int(parts[4])
                index[seqid] = (offset, slen, lbases, lbytes)
    return index


def extract_sequence(
    fasta_path: str,
    seqid: str,
    start: int,
    end: int,
    strand: str = "+",
    index: Dict[str, FAIEntry] = None,
) -> str:
    """Extract a subsequence from the FASTA.

    Args:
        fasta_path: Path to FASTA file.
        seqid: Sequence ID (chromosome/contig name).
        start: 1-based start coordinate (inclusive).
        end: 1-based end coordinate (inclusive).
        strand: "+" or "-". If "-", returns reverse complement.
        index: Pre-built FASTA index. If None, builds one.

    Returns:
        DNA sequence string in uppercase.
    """
    if index is None:
        index = build_fasta_index(fasta_path)

    if seqid not in index:
        raise KeyError(f"Sequence '{seqid}' not found in FASTA. Available: {list(index.keys())[:20]}...")

    offset, seq_len, line_bases, line_bytes = index[seqid]

    # Convert 1-based genomic coords to 0-based
    s0 = max(0, start - 1)
    e0 = min(seq_len, end)

    if s0 >= seq_len or e0 <= 0 or s0 >= e0:
        return ""

    # Calculate file positions
    seq_start = s0
    seq_end = e0 - 1

    # Line number where the sequence starts
    start_line = seq_start // line_bases
    start_pos_in_line = seq_start % line_bases

    end_line = seq_end // line_bases
    end_pos_in_line = seq_end % line_bases

    # Calculate byte offset
    byte_start = offset + start_line * line_bytes + start_pos_in_line
    byte_end = offset + end_line * line_bytes + end_pos_in_line + 1

    with open(fasta_path, "rb") as f:
        f.seek(byte_start)
        raw = f.read(byte_end - byte_start)

    seq = raw.decode("ascii", errors="ignore").replace("\n", "").replace("\r", "").upper()

    if strand == "-":
        seq = reverse_complement(seq)

    return seq


def reverse_complement(seq: str) -> str:
    """Return reverse complement of a DNA sequence."""
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N",
            "a": "t", "t": "a", "c": "g", "g": "c", "n": "n",
            "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
            "B": "V", "V": "B", "D": "H", "H": "D",
            "r": "y", "y": "r", "s": "s", "w": "w", "k": "m", "m": "k",
            "b": "v", "v": "b", "d": "h", "h": "d"}
    return "".join(comp.get(c, c) for c in reversed(seq))


def load_fasta_index(fasta_path: str) -> Dict[str, FAIEntry]:
    """Load or build FASTA index. Wrapper around build_fasta_index."""
    return build_fasta_index(fasta_path)
