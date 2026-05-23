"""Check user-provided primers for quality and specificity."""

from typing import Any, Dict, List, Optional

import primer3

from .constants import (
    DEFAULT_MONO_CATIONS,
    DEFAULT_DIVA_CATIONS,
    DEFAULT_CON_DNTPS,
)
from .primer_score import (
    score_single_primer,
    score_delta_tm,
    score_cross_dimer,
)


def calc_primer_tm(seq: str) -> float:
    """Calculate Tm for a primer sequence using Primer3's SantaLucia 1998 method."""
    return primer3.calc_tm(
        seq,
        mv_conc=DEFAULT_MONO_CATIONS,
        dv_conc=DEFAULT_DIVA_CATIONS,
        dntp_conc=DEFAULT_CON_DNTPS,
        dmso_conc=0,
        salt_corrections_method='santalucia',
        tm_method='santalucia',
    )


def calc_primer_gc(seq: str) -> float:
    """Calculate GC percentage for a primer sequence."""
    seq = seq.upper()
    gc = sum(1 for b in seq if b in ('G', 'C'))
    return round(gc / len(seq) * 100, 1) if seq else 0.0


def check_primers(
    forward_seq: str,
    reverse_seq: str,
    blast_hits: Optional[Dict[str, List]] = None,
    specificity_fn=None,
    **specificity_kwargs,
) -> Dict[str, Any]:
    """Evaluate a user-provided primer pair for quality and specificity.

    Args:
        forward_seq: Forward (left) primer sequence 5'→3'.
        reverse_seq: Reverse (right) primer sequence 5'→3'.
        blast_hits: Pre-computed BLAST hits dict {seq: [BlastHit, ...]}.
        specificity_fn: Callable to evaluate specificity (evaluate_specificity).
        **specificity_kwargs: Passed to specificity_fn.

    Returns:
        Dict with 'forward', 'reverse', 'pair_scores', and optionally 'specificity'.
    """
    fwd = forward_seq.strip().upper()
    rev = reverse_seq.strip().upper()

    if not fwd or not rev:
        raise ValueError("Both forward and reverse primer sequences are required.")

    for label, seq in [("Forward", fwd), ("Reverse", rev)]:
        if not all(b in 'ATGCRYSWKMBDHVN' for b in seq):
            raise ValueError(f"{label} primer contains invalid nucleotide characters.")
        if len(seq) < 12:
            raise ValueError(f"{label} primer is too short ({len(seq)} nt, minimum 12).")
        if len(seq) > 40:
            raise ValueError(f"{label} primer is too long ({len(seq)} nt, maximum 40).")

    fwd_tm = calc_primer_tm(fwd)
    rev_tm = calc_primer_tm(rev)
    fwd_gc = calc_primer_gc(fwd)
    rev_gc = calc_primer_gc(rev)

    fwd_scores = score_single_primer(fwd, fwd_tm, fwd_gc)
    rev_scores = score_single_primer(rev, rev_tm, rev_gc)

    cross = {
        'delta_tm': score_delta_tm(fwd_tm, rev_tm),
        'cross_dimer': score_cross_dimer(fwd, rev),
    }

    per_primer_avg = (fwd_scores['total'] + rev_scores['total']) / 2.0
    cross_total = sum(cross.values())
    overall = round(per_primer_avg + cross_total, 1)

    result = {
        'forward': {
            'sequence': fwd,
            'length': len(fwd),
            'tm': round(fwd_tm, 1),
            'gc_percent': fwd_gc,
            'scores': fwd_scores,
        },
        'reverse': {
            'sequence': rev,
            'length': len(rev),
            'tm': round(rev_tm, 1),
            'gc_percent': rev_gc,
            'scores': rev_scores,
        },
        'pair': {
            'total_score': overall,
            'per_primer_avg': round(per_primer_avg, 1),
            'cross_total': cross_total,
            'cross_scores': cross,
            'delta_tm': round(abs(fwd_tm - rev_tm), 1),
        },
    }

    if blast_hits and specificity_fn:
        left_hits = blast_hits.get(fwd, [])
        right_hits = blast_hits.get(rev, [])
        try:
            spec_result = specificity_fn(
                left_hits=left_hits,
                right_hits=right_hits,
                **specificity_kwargs,
            )
            result['specificity'] = {
                'is_specific': spec_result.is_specific,
                'off_target_count': spec_result.off_target_count,
                'off_target_amplicons': [
                    {
                        'chromosome': a.chromosome,
                        'product_size': a.product_size,
                        'gene_name': a.gene_name,
                        'left_hit': {
                            'subject_start': a.left_hit.subject_start,
                            'subject_end': a.left_hit.subject_end,
                            'num_mismatch': a.left_hit.num_mismatch,
                        },
                        'right_hit': {
                            'subject_start': a.right_hit.subject_start,
                            'subject_end': a.right_hit.subject_end,
                            'num_mismatch': a.right_hit.num_mismatch,
                        },
                    }
                    for a in spec_result.off_target_amplicons
                ],
            }
        except Exception:
            result['specificity'] = {'error': 'Specificity evaluation failed'}

    return result
