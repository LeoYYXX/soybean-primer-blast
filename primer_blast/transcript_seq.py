"""Transcript (cDNA) FASTA index and cDNA amplicon extraction.

The transcript FASTA contains spliced mRNA sequences (no introns).
This module indexes the file for O(1) random access and provides
genomic-to-cDNA coordinate conversion for computing cDNA amplicon products.
"""

import os
from typing import Dict, List, Optional, Tuple

from .fasta_index import reverse_complement


def build_transcript_seq_index(fasta_path: str) -> Dict[str, Tuple[int, int, int, int]]:
    """Index the transcript FASTA. Returns {transcript_id: (offset, len, line_bases, line_bytes)}."""
    index: Dict[str, Tuple[int, int, int, int]] = {}
    fai_path = fasta_path + ".fai"

    if os.path.exists(fai_path):
        if os.path.getmtime(fai_path) >= os.path.getmtime(fasta_path):
            return _load_transcript_fai(fai_path)

    with open(fasta_path, "rb") as f:
        offset = 0
        current_id = ""
        seq_length = 0
        line_bases = 0
        line_bytes = 0
        started = False

        while True:
            line = f.readline()
            if not line:
                if started and current_id:
                    index[current_id] = (offset, seq_length, line_bases, line_bytes)
                break

            line_str = line.decode("ascii", errors="ignore").rstrip("\r\n")

            if line.startswith(b">"):
                if started and current_id:
                    index[current_id] = (offset, seq_length, line_bases, line_bytes)
                # Parse: >Glyma.19G000100.3 pacid=... locus=... ID=Glyma.19G000100.3.Wm82.a4.v1 ...
                header = line_str[1:]
                current_id = _parse_transcript_id(header)
                offset = f.tell()
                seq_length = 0
                line_bases = 0
                line_bytes = 0
                started = True
            elif started and line_str:
                if line_bases == 0:
                    line_bases = len(line_str)
                    line_bytes = len(line)
                seq_length += len(line_str)

    with open(fai_path, "w") as out:
        for tid, (off, slen, lbases, lbytes) in index.items():
            out.write(f"{tid}\t{slen}\t{off}\t{lbases}\t{lbytes}\n")

    return index


def _load_transcript_fai(fai_path: str) -> Dict[str, Tuple[int, int, int, int]]:
    """Load an existing transcript .fai index."""
    index: Dict[str, Tuple[int, int, int, int]] = {}
    with open(fai_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 5:
                index[parts[0]] = (int(parts[1]), int(parts[2]),
                                   int(parts[3]), int(parts[4]))
    return index


def _parse_transcript_id(header: str) -> str:
    """Extract the full transcript ID from a FASTA header."""
    # Format: >Glyma.19G000100.3 pacid=41161166 locus=Glyma.19G000100 ID=Glyma.19G000100.3.Wm82.a4.v1 ...
    # Return the full ID= value if available, else the first token
    for part in header.split():
        if part.startswith("ID="):
            return part[3:]
    return header.split()[0]


def extract_transcript_sequence(
    fasta_path: str,
    transcript_id: str,
    index: Dict[str, Tuple[int, int, int, int]] = None,
) -> str:
    """Extract full spliced transcript (cDNA) sequence by transcript ID."""
    if index is None:
        index = build_transcript_seq_index(fasta_path)

    if transcript_id not in index:
        raise KeyError(f"Transcript '{transcript_id}' not found in transcript FASTA")

    slen, offset, lbases, lbytes = index[transcript_id]
    start_byte = offset
    end_byte = offset + ((slen - 1) // lbases) * lbytes + (slen - 1) % lbases + 1

    with open(fasta_path, "rb") as f:
        f.seek(start_byte)
        raw = f.read(end_byte - start_byte)

    return raw.decode("ascii", errors="ignore").replace("\n", "").replace("\r", "").upper()


def genomic_to_cdna_coords(
    genomic_pos: int,
    exons: List[Tuple[int, int]],
    strand: str,
) -> Optional[int]:
    """Map a genomic coordinate to a 0-based cDNA (spliced) coordinate.

    Args:
        genomic_pos: Position on the genome (1-based).
        exons: List of (start, end) genomic exon coordinates, sorted ascending (genomic order).
        strand: "+" or "-".

    Returns:
        0-based position within the spliced cDNA, or None if outside all exons.
    """
    cumulative = 0

    if strand == "+":
        for ex_start, ex_end in exons:
            if ex_start <= genomic_pos <= ex_end:
                return cumulative + (genomic_pos - ex_start)
            cumulative += (ex_end - ex_start + 1)
    else:
        # For minus strand: genomic order is 3'→5' relative to cDNA
        # Reverse exon order so we iterate 5'→3' in cDNA direction
        for ex_start, ex_end in exons:
            if ex_start <= genomic_pos <= ex_end:
                return cumulative + (ex_end - genomic_pos)
            cumulative += (ex_end - ex_start + 1)

    return None


def compute_cdna_amplicon(
    fasta_path: str,
    transcript_id: str,
    exons: List[Tuple[int, int]],
    strand: str,
    left_genomic_pos: int,   # 1-based, forward primer binding position
    right_genomic_pos: int,  # 1-based, reverse primer binding position
    index: Dict[str, Tuple[int, int, int, int]] = None,
) -> Optional[dict]:
    """Compute the cDNA (spliced) amplicon for a primer pair.

    Returns dict with cdna_sequence, cdna_length, left_cdna_pos, right_cdna_pos,
    or None if the primers don't both land in exons.
    """
    if index is None:
        index = build_transcript_seq_index(fasta_path)

    if transcript_id not in index:
        return None

    if not exons:
        return None

    # Map genomic binding positions to cDNA coordinates
    left_cdna = genomic_to_cdna_coords(left_genomic_pos, exons, strand)
    right_cdna = genomic_to_cdna_coords(right_genomic_pos, exons, strand)

    if left_cdna is None or right_cdna is None:
        return None

    # Ensure left < right in cDNA space
    if left_cdna > right_cdna:
        left_cdna, right_cdna = right_cdna, left_cdna

    full_cdna = extract_transcript_sequence(fasta_path, transcript_id, index)
    cdna_length = right_cdna - left_cdna + 1

    if cdna_length <= 0 or left_cdna >= len(full_cdna):
        return None

    amplicon_seq = full_cdna[left_cdna:right_cdna + 1]

    return {
        "cdna_sequence": amplicon_seq,
        "cdna_length": len(amplicon_seq),
        "left_cdna_position": left_cdna,
        "right_cdna_position": right_cdna,
    }
