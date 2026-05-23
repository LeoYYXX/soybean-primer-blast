"""Primer specificity checking via BLAST against the soybean genome.

Replicates NCBI Primer-BLAST's approach:
1. BLAST each primer against the genome
2. Parse hits and check 3'-end complementarity
3. Detect off-target primer-pair amplicons
"""

import os
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

from .constants import (
    DEFAULT_BLAST_EVALUE,
    DEFAULT_BLAST_MAX_TARGET_SEQS,
    DEFAULT_BLAST_WORD_SIZE,
    DEFAULT_MAX_AMPLICON_SIZE,
    DEFAULT_MAX_OFF_TARGET_MISMATCH,
    DEFAULT_MIN_3PRIME_MATCH,
)


@dataclass
class BlastHit:
    """A single BLAST hit for a primer against a genome sequence."""
    primer_seq: str
    subject_id: str
    query_start: int
    query_end: int
    subject_start: int   # genomic coordinate
    subject_end: int     # genomic coordinate
    evalue: float
    bit_score: float
    alignment_length: int
    num_mismatch: int
    num_gap: int
    qseq: str            # aligned query portion
    sseq: str            # aligned subject portion

    @property
    def is_reverse(self) -> bool:
        """True if hit is on the reverse strand (subject_end < subject_start)."""
        return self.subject_end < self.subject_start


@dataclass
class OffTargetAmplicon:
    """A potential off-target amplicon from a primer pair."""
    left_hit: BlastHit
    right_hit: BlastHit
    product_size: int
    chromosome: str
    gene_name: str = ""


@dataclass
class SpecificityResult:
    """Specificity check result for a primer pair."""
    is_specific: bool
    off_target_count: int
    off_target_amplicons: List[OffTargetAmplicon]
    left_hits: List[BlastHit]
    right_hits: List[BlastHit]
    left_target_hit: Optional[BlastHit] = None
    right_target_hit: Optional[BlastHit] = None


def run_blast_primers(
    primer_sequences: List[str],
    blast_db: str,
    blastn_path: str = "blastn",
    evalue: float = DEFAULT_BLAST_EVALUE,
    word_size: int = DEFAULT_BLAST_WORD_SIZE,
    max_target_seqs: int = DEFAULT_BLAST_MAX_TARGET_SEQS,
    num_threads: int = 2,
    timeout: int = 120,
) -> Dict[str, List[BlastHit]]:
    """Run BLAST for all primer sequences against the genome.

    Args:
        primer_sequences: List of primer sequences to BLAST.
        blast_db: Path prefix of the BLAST database.
        blastn_path: Path to blastn executable.
        evalue: BLAST E-value threshold.
        word_size: BLAST word size.
        max_target_seqs: Max hits per query.
        num_threads: Number of CPU threads.
        timeout: Subprocess timeout in seconds.

    Returns:
        Dict mapping primer_sequence -> list of BlastHit objects.
    """
    if not primer_sequences:
        return {}

    # Write primer sequences to temp FASTA
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fa", delete=False, prefix="primers_"
    ) as f:
        for i, seq in enumerate(primer_sequences):
            f.write(f">primer_{i}\n{seq}\n")
        query_file = f.name

    try:
        # Custom outfmt capturing alignment details needed for 3' end analysis
        outfmt = "6 qseqid sseqid qstart qend sstart send evalue bitscore length mismatch gapopen qseq sseq"
        cmd = [
            blastn_path,
            "-task", "blastn-short",
            "-db", blast_db,
            "-query", query_file,
            "-outfmt", outfmt,
            "-evalue", str(evalue),
            "-word_size", str(word_size),
            "-dust", "no",
            "-max_target_seqs", str(max_target_seqs),
            "-num_threads", str(num_threads),
        ]

        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout
        )

        if result.returncode != 0:
            raise RuntimeError(f"blastn failed: {result.stderr}")

        return _parse_blast_output(result.stdout, primer_sequences)

    finally:
        os.unlink(query_file)


def _parse_blast_output(
    output: str,
    primer_sequences: List[str],
) -> Dict[str, List[BlastHit]]:
    """Parse blastn -outfmt 6 tabular output into BlastHit dict."""
    # Build mapping from query_id -> primer sequence
    id_to_seq = {f"primer_{i}": seq for i, seq in enumerate(primer_sequences)}

    hits: Dict[str, List[BlastHit]] = defaultdict(list)

    for line in output.strip().split("\n"):
        if not line.strip():
            continue
        parts = line.strip().split("\t")
        if len(parts) < 13:
            continue

        query_id = parts[0]
        primer_seq = id_to_seq.get(query_id)
        if not primer_seq:
            continue

        try:
            hit = BlastHit(
                primer_seq=primer_seq,
                subject_id=parts[1],
                query_start=int(parts[2]),
                query_end=int(parts[3]),
                subject_start=int(parts[4]),
                subject_end=int(parts[5]),
                evalue=float(parts[6]),
                bit_score=float(parts[7]),
                alignment_length=int(parts[8]),
                num_mismatch=int(parts[9]),
                num_gap=int(parts[10]),
                qseq=parts[11],
                sseq=parts[12],
            )
            hits[primer_seq].append(hit)
        except (ValueError, IndexError):
            continue

    return dict(hits)


def check_3prime_match(hit: BlastHit, min_match_bases: int = DEFAULT_MIN_3PRIME_MATCH) -> int:
    """Count consecutive matching bases at the 3' end of the primer.

    Uses the BLAST hit orientation (is_reverse) to determine which end of
    the query alignment corresponds to the primer's 3' end.

    Returns:
        Number of consecutive matches from the 3' end.
    """
    qseq = hit.qseq.upper()
    sseq = hit.sseq.upper()

    if not qseq or not sseq:
        return 0

    # For forward hits (subject_start < subject_end): 3' end = right side of alignment (query_end)
    # For reverse hits (subject_end < subject_start): 3' end = left side of alignment (query_start)
    if not hit.is_reverse:
        matches = 0
        for i in range(1, min(len(qseq), len(sseq)) + 1):
            if qseq[-i] == sseq[-i] or qseq[-i] == "-" or sseq[-i] == "-":
                matches += 1
            else:
                break
        return matches
    else:
        matches = 0
        for i in range(min(len(qseq), len(sseq))):
            if qseq[i] == sseq[i] or qseq[i] == "-" or sseq[i] == "-":
                matches += 1
            else:
                break
        return matches


def _alignment_cover_fraction(hit: BlastHit) -> float:
    """Return what fraction of the full primer length is covered by this alignment.

    Uses query_start/query_end to determine alignment span vs full primer length.
    """
    primer_len = len(hit.primer_seq)
    if primer_len == 0:
        return 0.0
    qspan = abs(hit.query_end - hit.query_start) + 1
    return qspan / primer_len


def filter_specific_hits(
    hits: List[BlastHit],
    target_seqid: str,
    target_coords: Tuple[int, int],
    max_mismatches: int = DEFAULT_MAX_OFF_TARGET_MISMATCH,
    min_3prime_match: int = DEFAULT_MIN_3PRIME_MATCH,
    min_alignment_cover: float = 0.75,
    tolerance: int = 50,
) -> List[BlastHit]:
    """Filter hits to only those that are potential off-target binding sites.

    A hit is considered an off-target risk if:
    1. It is NOT the intended target site
    2. Total mismatches + gaps <= max_mismatches
    3. The 3' end has >= min_3prime_match consecutive matches
    4. The alignment covers >= min_alignment_cover fraction of the primer

    Args:
        hits: All BLAST hits for a primer.
        target_seqid: Expected chromosome/contig of the intended binding site.
        target_coords: (start, end) of the intended binding site on the genome.
        max_mismatches: Max allowed mismatches+gaps for off-target concern.
        min_3prime_match: Min 3' end consecutive matches required.
        min_alignment_cover: Min fraction of primer covered by alignment (0.0-1.0).
        tolerance: bp tolerance for matching the intended target coordinates.

    Returns:
        List of hits that are potential off-target concerns.
    """
    tgt_start, tgt_end = target_coords
    off_targets = []

    for hit in hits:
        # Check if this is the intended target
        if hit.subject_id == target_seqid:
            hit_mid = (hit.subject_start + hit.subject_end) // 2
            tgt_mid = (tgt_start + tgt_end) // 2
            if abs(hit_mid - tgt_mid) < tolerance:
                continue  # This is the intended target

        # Check mismatch threshold
        if hit.num_mismatch + hit.num_gap > max_mismatches:
            continue

        # Check alignment coverage (skip partial short alignments)
        if _alignment_cover_fraction(hit) < min_alignment_cover:
            continue

        # Check 3' end match
        if check_3prime_match(hit, min_3prime_match) < min_3prime_match:
            continue

        off_targets.append(hit)

    return off_targets


def find_target_hit(
    hits: List[BlastHit],
    target_seqid: str,
    target_coords: Tuple[int, int],
    tolerance: int = 100,
) -> Optional[BlastHit]:
    """Find the BLAST hit corresponding to the intended target site."""
    tgt_start, tgt_end = target_coords
    tgt_mid = (tgt_start + tgt_end) // 2

    best: Optional[BlastHit] = None
    best_dist = float("inf")

    for hit in hits:
        if hit.subject_id == target_seqid:
            hit_mid = (hit.subject_start + hit.subject_end) // 2
            dist = abs(hit_mid - tgt_mid)
            if dist < tolerance and dist < best_dist:
                best = hit
                best_dist = dist

    return best


def check_pair_specificity(
    left_hits: List[BlastHit],
    right_hits: List[BlastHit],
    max_amplicon_size: int = DEFAULT_MAX_AMPLICON_SIZE,
    min_off_target_product: int = 50,
    genomic_gene_map: dict = None,
) -> List[OffTargetAmplicon]:
    """Check if any left+right hit combination could produce an off-target amplicon.

    An off-target amplicon is possible when:
    1. Both primers hit the SAME chromosome/contig
    2. The primers bind in opposite orientations (one forward, one reverse)
    3. Their binding sites are within max_amplicon_size and above min_off_target_product
    4. Their 3' ends face each other

    If genomic_gene_map is provided (from gff_index.build_genomic_gene_map),
    each OffTargetAmplicon will be annotated with the overlapping gene name.

    Returns:
        List of OffTargetAmplicon objects.
    """
    # Group hits by chromosome
    left_by_chr: Dict[str, List[BlastHit]] = defaultdict(list)
    right_by_chr: Dict[str, List[BlastHit]] = defaultdict(list)

    for hit in left_hits:
        left_by_chr[hit.subject_id].append(hit)
    for hit in right_hits:
        right_by_chr[hit.subject_id].append(hit)

    off_targets: List[OffTargetAmplicon] = []

    # Only check chromosomes where both primers have hits
    common_chrs = set(left_by_chr.keys()) & set(right_by_chr.keys())

    for chr_id in common_chrs:
        for l_hit in left_by_chr[chr_id]:
            for r_hit in right_by_chr[chr_id]:
                # Determine strand orientation for each hit
                l_on_reverse = l_hit.is_reverse
                r_on_reverse = r_hit.is_reverse

                # For a valid PCR amplicon, primers must bind to opposite strands
                if l_on_reverse == r_on_reverse:
                    continue

                # Calculate the 3' end positions of each primer
                if l_on_reverse:
                    l_3prime = l_hit.subject_start
                else:
                    l_3prime = l_hit.subject_end

                if r_on_reverse:
                    r_3prime = r_hit.subject_start
                else:
                    r_3prime = r_hit.subject_end

                # The two 3' ends must face each other
                if l_on_reverse and not r_on_reverse:
                    if l_3prime < r_3prime:
                        product_size = r_3prime - l_3prime
                        if min_off_target_product <= product_size <= max_amplicon_size:
                            gene_name = _lookup_gene(l_hit, r_hit, chr_id, genomic_gene_map)
                            off_targets.append(OffTargetAmplicon(
                                left_hit=l_hit, right_hit=r_hit,
                                product_size=product_size, chromosome=chr_id,
                                gene_name=gene_name,
                            ))
                elif not l_on_reverse and r_on_reverse:
                    if r_3prime < l_3prime:
                        product_size = l_3prime - r_3prime
                        if min_off_target_product <= product_size <= max_amplicon_size:
                            gene_name = _lookup_gene(l_hit, r_hit, chr_id, genomic_gene_map)
                            off_targets.append(OffTargetAmplicon(
                                left_hit=l_hit, right_hit=r_hit,
                                product_size=product_size, chromosome=chr_id,
                                gene_name=gene_name,
                            ))

    return off_targets


def _lookup_gene(l_hit, r_hit, chr_id, genomic_gene_map):
    """Look up gene name for an off-target amplicon position."""
    if not genomic_gene_map:
        return ""
    from .gff_index import find_gene_at_position_fast
    mid_pos = (l_hit.subject_start + r_hit.subject_end) // 2
    gene_name = find_gene_at_position_fast(chr_id, mid_pos, genomic_gene_map)
    return gene_name or ""


def evaluate_specificity(
    left_hits: List[BlastHit],
    right_hits: List[BlastHit],
    target_seqid: str,
    target_left_coords: Tuple[int, int],
    target_right_coords: Tuple[int, int],
    max_mismatches: int = DEFAULT_MAX_OFF_TARGET_MISMATCH,
    min_3prime_match: int = DEFAULT_MIN_3PRIME_MATCH,
    min_alignment_cover: float = 0.75,
    max_amplicon_size: int = DEFAULT_MAX_AMPLICON_SIZE,
    genomic_gene_map: dict = None,
) -> SpecificityResult:
    """Perform full specificity evaluation for a primer pair.

    Args:
        left_hits: All BLAST hits for the left (forward) primer.
        right_hits: All BLAST hits for the right (reverse) primer.
        target_seqid: Expected chromosome/contig.
        target_left_coords: (start, end) of left primer's intended binding site.
        target_right_coords: (start, end) of right primer's intended binding site.
        max_mismatches: Max mismatches allowed for off-target concern.
        min_3prime_match: Min consecutive 3' end matches required.
        min_alignment_cover: Min fraction of primer covered by alignment.
        max_amplicon_size: Max amplicon size for off-target detection.

    Returns:
        SpecificityResult with full details.
    """
    # Find intended target hits
    left_target = find_target_hit(left_hits, target_seqid, target_left_coords)
    right_target = find_target_hit(right_hits, target_seqid, target_right_coords)

    # Filter for potential off-target hits
    left_off = filter_specific_hits(
        left_hits, target_seqid, target_left_coords,
        max_mismatches, min_3prime_match, min_alignment_cover,
    )
    right_off = filter_specific_hits(
        right_hits, target_seqid, target_right_coords,
        max_mismatches, min_3prime_match, min_alignment_cover,
    )

    # Check pair-level off-target amplicons
    off_amplicons = check_pair_specificity(
        left_off, right_off, max_amplicon_size,
        genomic_gene_map=genomic_gene_map,
    )

    is_specific = len(off_amplicons) == 0

    return SpecificityResult(
        is_specific=is_specific,
        off_target_count=len(off_amplicons),
        off_target_amplicons=off_amplicons,
        left_hits=left_hits,
        right_hits=right_hits,
        left_target_hit=left_target,
        right_target_hit=right_target,
    )
