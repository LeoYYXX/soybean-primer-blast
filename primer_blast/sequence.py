"""Sequence extraction for genes, CDS, UTR regions from soybean genome."""

from typing import Dict, Optional, Tuple

from .fasta_index import extract_sequence, FAIEntry
from .gff_index import TranscriptRecord


def extract_gene_sequence(
    fasta_path: str,
    transcript: TranscriptRecord,
    index: Dict[str, FAIEntry] = None,
) -> str:
    """Extract the full mRNA span (gene body) from the genome."""
    return extract_sequence(
        fasta_path=fasta_path,
        seqid=transcript.seqid,
        start=transcript.start,
        end=transcript.end,
        strand=transcript.strand,
        index=index,
    )


def extract_cds_sequence(
    fasta_path: str,
    transcript: TranscriptRecord,
    index: Dict[str, FAIEntry] = None,
) -> str:
    """Extract and concatenate CDS exons in coding order.

    For plus-strand genes: genomic order = coding order.
    For minus-strand genes: CDS exons are extracted from genomic coordinates
    individually and reverse-complemented, then joined in genomic order
    (which is also coding order since the exons are already reversed).
    """
    if not transcript.cds_exons:
        return ""

    parts = []
    for cds_start, cds_end in transcript.cds_exons:
        part = extract_sequence(
            fasta_path=fasta_path,
            seqid=transcript.seqid,
            start=cds_start,
            end=cds_end,
            strand=transcript.strand,
            index=index,
        )
        parts.append(part)

    # For plus strand: concatenate in genomic order (5'->3')
    # For minus strand: the individual parts above are already reverse-complemented,
    # but we extracted them in genomic order (start->end ascending).
    # For minus strand genes, genomic order is 3' to 5' relative to transcript,
    # so we need to reverse the order of exons.
    if transcript.strand == "-":
        parts.reverse()

    return "".join(parts)


def extract_region_with_padding(
    fasta_path: str,
    seqid: str,
    start: int,
    end: int,
    strand: str = "+",
    upstream_pad: int = 0,
    downstream_pad: int = 0,
    index: Dict[str, FAIEntry] = None,
) -> str:
    """Extract genomic region with flanking padding."""
    total_start = max(1, start - upstream_pad)
    total_end = end + downstream_pad
    return extract_sequence(
        fasta_path=fasta_path,
        seqid=seqid,
        start=total_start,
        end=total_end,
        strand=strand,
        index=index,
    )


def get_target_coordinates(
    transcript: TranscriptRecord,
    target: str,
) -> Optional[Tuple[int, int]]:
    """Get genomic coordinates for a named target region.

    Args:
        transcript: The transcript record.
        target: One of 'cds', 'utr', 'all', or 'start,end' coords.

    Returns:
        (start, end) genomic coordinates, or None if not applicable.
    """
    target_lower = target.lower().strip()

    # Try parsing as numeric coordinates "start,end"
    if "," in target_lower:
        parts = target_lower.split(",")
        if len(parts) == 2:
            try:
                t_start = int(parts[0].strip())
                t_end = int(parts[1].strip())
                return (min(t_start, t_end), max(t_start, t_end))
            except ValueError:
                pass

    if transcript is None:
        # For region/sequence mode with numeric targets, only "," parsing works
        return None

    if target_lower == "cds":
        if transcript.cds_start is not None:
            return (transcript.cds_start, transcript.cds_end)
        return None

    if target_lower == "utr":
        all_utrs = transcript.five_prime_utr + transcript.three_prime_utr
        if all_utrs:
            utr_start = min(u[0] for u in all_utrs)
            utr_end = max(u[1] for u in all_utrs)
            return (utr_start, utr_end)
        return None

    if target_lower in ("all", "gene", "body"):
        return (transcript.start, transcript.end)

    return None


def genomic_to_template_coords(
    genomic_start: int,
    genomic_end: int,
    template_start: int,
    template_end: int,
    strand: str,
) -> Tuple[int, int]:
    """Convert genomic coordinates to 0-based template-relative coordinates.

    Returns (start_0based, length) suitable for Primer3 SEQUENCE_TARGET.
    """
    if strand == "+":
        rel_start = genomic_start - template_start
        rel_end = genomic_end - template_start
    else:
        # For minus strand, template is reverse-complemented
        rel_start = template_end - genomic_end
        rel_end = template_end - genomic_start

    return (rel_start, rel_end - rel_start + 1)
