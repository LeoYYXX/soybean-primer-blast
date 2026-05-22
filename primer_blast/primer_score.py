"""Primer scoring system — 100-point scale across 7 categories.

Categories:
  1. 3' terminal base (10 pts) — G/C preference, GC clamp penalty
  2. 3' end 8-base stability (10 pts) — ΔG of last 8 nt
  3. 3' hairpin (10 pts) — hairpin involvement at 3' end
  4. Tm + ΔTm (20 pts) — optimal Tm range + F/R Tm difference
  5. GC content (15 pts) — optimal GC% range
  6. Base distribution (15 pts) — sliding window GC%, homopolymer runs
  7. Dimers (20 pts) — self-dimers + cross-dimer 3' complement
"""

import primer3


def _parse_structure_3prime_involvement(structure_lines, seq_len):
    """Count how many 3'-end bases participate in a secondary structure.

    structure_lines: list of strings from primer3 ascii_structure_lines.
    Returns (num_3prime_bases, total_paired_bases).
    """
    if not structure_lines or len(structure_lines) < 2:
        return 0, 0

    # The last line shows structure symbols; the second-to-last shows the sequence
    struct_line = structure_lines[-1]
    paired_positions = 0
    three_prime_paired = 0

    for i, ch in enumerate(reversed(struct_line)):
        if i >= seq_len:
            break
        if ch != ' ':
            paired_positions += 1
            three_prime_paired += 1
        else:
            # Stop counting 3' involvement at first unpaired base from 3' end
            break

    return three_prime_paired, paired_positions


def _max_consecutive_pairs(structure_lines):
    """Find maximum consecutive paired bases from structure_lines."""
    if not structure_lines or len(structure_lines) < 2:
        return 0
    struct_line = structure_lines[-1]
    max_run = 0
    cur_run = 0
    for ch in struct_line:
        if ch != ' ':
            cur_run += 1
            if cur_run > max_run:
                max_run = cur_run
        else:
            cur_run = 0
    return max_run


def _count_polyx_runs(seq):
    """Count homopolymer runs. Returns max_run and list of (base, run_length)."""
    if not seq:
        return 0, []
    runs = []
    cur_base = seq[0]
    cur_len = 1
    max_run = 1
    for i in range(1, len(seq)):
        if seq[i] == cur_base:
            cur_len += 1
        else:
            if cur_len >= 2:
                runs.append((cur_base, cur_len))
            max_run = max(max_run, cur_len)
            cur_base = seq[i]
            cur_len = 1
    if cur_len >= 2:
        runs.append((cur_base, cur_len))
    max_run = max(max_run, cur_len)
    return max_run, runs


# ---------------------------------------------------------------------------
# Category 1: 3' terminal base (10 pts)
# ---------------------------------------------------------------------------
def score_terminal_base(seq):
    """Score the 3'-most base. G/C = 10, A/T = 0, GC clamp >=4 in a row = −5."""
    last = seq[-1].upper() if seq else ''
    score = 10 if last in ('G', 'C') else 0

    # Check for excessive GC clamp at 3' end (>=4 consecutive G or C)
    g_count = 0
    c_count = 0
    for base in reversed(seq):
        if base == 'G':
            g_count += 1
            c_count = 0
        elif base == 'C':
            c_count += 1
            g_count = 0
        else:
            break
    if g_count >= 4 or c_count >= 4:
        score -= 5

    return max(0, score)


# ---------------------------------------------------------------------------
# Category 2: 3' end 8-base stability (10 pts)
# ---------------------------------------------------------------------------
def score_end_stability(seq):
    """Score based on ΔG of last 8 bases. Lower ΔG (more stable) = worse."""
    last8 = seq[-8:] if len(seq) >= 8 else seq
    if len(last8) < 4:
        return 5

    try:
        dg = primer3.calc_end_stability(last8, last8).dg
    except Exception:
        return 5

    dg_kcal = abs(dg) / 1000.0 if dg else 0

    if dg_kcal > 9:
        return 10
    elif dg_kcal >= 6:
        return 6
    else:
        return 0


# ---------------------------------------------------------------------------
# Category 3: 3' hairpin (10 pts)
# ---------------------------------------------------------------------------
def score_hairpin(seq):
    """Score based on hairpin formation at 3' end."""
    try:
        hp = primer3.calc_hairpin(seq)
    except Exception:
        return 5

    if not hp.structure_found or hp.dg > -2000:
        return 10

    three_paired, _ = _parse_structure_3prime_involvement(
        hp.ascii_structure_lines, len(seq)
    )
    dg_kcal = abs(hp.dg) / 1000.0 if hp.dg else 0

    if three_paired == 0 or dg_kcal < 2:
        return 10
    elif three_paired <= 2 or dg_kcal < 3:
        return 5
    else:
        return 0


# ---------------------------------------------------------------------------
# Category 4a: Single primer Tm (10 pts)
# ---------------------------------------------------------------------------
def score_tm(tm):
    """Score based on primer Tm."""
    if 58 <= tm <= 62:
        return 10
    elif 55 <= tm < 58 or 62 < tm <= 65:
        return 7
    elif 53 <= tm < 55 or 65 < tm <= 67:
        return 3
    else:
        return 0


# ---------------------------------------------------------------------------
# Category 4b: F/R Tm difference (10 pts)
# ---------------------------------------------------------------------------
def score_delta_tm(tm_f, tm_r):
    """Score based on absolute Tm difference between forward and reverse."""
    delta = abs(tm_f - tm_r)
    if delta <= 1:
        return 10
    elif delta <= 2:
        return 6
    elif delta <= 3:
        return 2
    else:
        return 0


# ---------------------------------------------------------------------------
# Category 5: GC content (15 pts)
# ---------------------------------------------------------------------------
def score_gc(gc_percent):
    """Score based on GC percentage."""
    if 45 <= gc_percent <= 55:
        return 15
    elif 40 <= gc_percent < 45 or 55 < gc_percent <= 60:
        return 11
    elif 35 <= gc_percent < 40 or 60 < gc_percent <= 65:
        return 5
    else:
        return 0


# ---------------------------------------------------------------------------
# Category 6: Base distribution uniformity (15 pts)
# ---------------------------------------------------------------------------
def score_base_distribution(seq):
    """Score based on 5-nt sliding window GC% and homopolymer runs."""
    seq = seq.upper()
    n = len(seq)
    bad_windows = 0  # windows with GC=0% or 100%

    for i in range(n - 4):
        window = seq[i:i + 5]
        gc_count = window.count('G') + window.count('C')
        if gc_count == 0 or gc_count == 5:
            bad_windows += 1

    if bad_windows == 0:
        score = 15
    elif bad_windows == 1:
        score = 8
    else:
        score = 0

    # Homopolymer run penalties
    max_run, _ = _count_polyx_runs(seq)
    if max_run >= 6:
        score = 0
    elif max_run == 5:
        score = max(0, score - 8)
    elif max_run == 4:
        score = max(0, score - 3)

    return score


# ---------------------------------------------------------------------------
# Category 7a: Self-dimer (10 pts per primer, averaged)
# ---------------------------------------------------------------------------
def score_self_dimer(seq):
    """Score based on maximum consecutive complementary bases in homodimer."""
    try:
        hd = primer3.calc_homodimer(seq)
    except Exception:
        return 5

    if not hd.structure_found:
        return 10

    max_consec = _max_consecutive_pairs(hd.ascii_structure_lines)
    if max_consec <= 3:
        return 10
    elif max_consec == 4:
        return 5
    else:
        return 0


# ---------------------------------------------------------------------------
# Category 7b: F/R cross-dimer 3' complement (10 pts)
# ---------------------------------------------------------------------------
def score_cross_dimer(seq_f, seq_r):
    """Score based on 3'-end complementarity between forward and reverse primers."""
    try:
        het = primer3.calc_heterodimer(seq_f, seq_r)
    except Exception:
        return 5

    if not het.structure_found:
        return 10

    # Count 3'-end paired bases (check both primers' 3' ends)
    lines = het.ascii_structure_lines
    if not lines or len(lines) < 2:
        return 10

    struct = lines[-1]
    # Check how many bases at the 3' end of each primer are paired
    # The structure line shows paired bases; check the right end
    three_prime_paired = 0
    for ch in reversed(struct):
        if ch != ' ':
            three_prime_paired += 1
        else:
            break

    if three_prime_paired <= 2:
        return 10
    elif three_prime_paired == 3:
        return 4
    else:
        return 0


# ---------------------------------------------------------------------------
# Aggregate scoring
# ---------------------------------------------------------------------------
def score_single_primer(seq, tm, gc_percent):
    """Score a single primer on all per-primer categories.

    Returns dict of individual category scores and total.
    """
    scores = {
        'terminal_base': score_terminal_base(seq),
        'end_stability': score_end_stability(seq),
        'hairpin': score_hairpin(seq),
        'tm': score_tm(tm),
        'gc_content': score_gc(gc_percent),
        'base_distribution': score_base_distribution(seq),
        'self_dimer': score_self_dimer(seq),
    }
    scores['total'] = sum(scores.values())
    return scores


def score_primer_pair(pair):
    """Score a PrimerPair on the 100-point scale.

    Args:
        pair: PrimerPair dataclass with left_primer, right_primer, left_tm,
              right_tm, left_gc, right_gc.

    Returns:
        dict with 'total', 'forward', 'reverse', 'cross', and breakdown details.
    """
    fwd = score_single_primer(pair.left_primer, pair.left_tm, pair.left_gc)
    rev = score_single_primer(pair.right_primer, pair.right_tm, pair.right_gc)

    cross = {
        'delta_tm': score_delta_tm(pair.left_tm, pair.right_tm),
        'cross_dimer': score_cross_dimer(pair.left_primer, pair.right_primer),
    }

    # Per-primer items are averaged between F and R
    per_primer_total = (fwd['total'] + rev['total']) / 2.0
    cross_total = sum(cross.values())
    overall = round(per_primer_total + cross_total, 1)

    return {
        'total': overall,
        'forward': fwd,
        'reverse': rev,
        'cross': cross,
        'per_primer_avg': round(per_primer_total, 1),
        'cross_total': cross_total,
    }
