"""Primer3 wrapper for primer design with soybean-appropriate defaults."""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import primer3

from .constants import (
    DEFAULT_CON_DNTPS,
    DEFAULT_DIVA_CATIONS,
    DEFAULT_GC_CLAMP,
    DEFAULT_MAX_END_GC,
    DEFAULT_MAX_END_STABILITY,
    DEFAULT_MAX_POLY_X,
    DEFAULT_MONO_CATIONS,
    DEFAULT_NUM_RETURN,
    DEFAULT_PRIMER_GC_MAX,
    DEFAULT_PRIMER_GC_MIN,
    DEFAULT_PRIMER_MAX_DIFF_TM,
    DEFAULT_PRIMER_MAX_SIZE,
    DEFAULT_PRIMER_MAX_TM,
    DEFAULT_PRIMER_MIN_SIZE,
    DEFAULT_PRIMER_MIN_TM,
    DEFAULT_PRIMER_OPT_SIZE,
    DEFAULT_PRIMER_OPT_TM,
    DEFAULT_PRIMER_PRODUCT_MAX,
    DEFAULT_PRIMER_PRODUCT_MIN,
    DEFAULT_SALT_CORRECTION,
    DEFAULT_TM_METHOD,
)


@dataclass
class PrimerPair:
    """A single primer pair designed by Primer3."""
    left_primer: str
    right_primer: str
    left_tm: float
    right_tm: float
    left_gc: float
    right_gc: float
    product_size: int
    left_position: int          # 0-based within template
    right_position: int         # 0-based within template
    penalty: float
    pair_product_tm: float = 0.0
    left_self_any_th: float = 0.0
    left_self_end_th: float = 0.0
    right_self_any_th: float = 0.0
    right_self_end_th: float = 0.0
    pair_compl_any_th: float = 0.0
    pair_compl_end_th: float = 0.0
    left_hairpin_th: float = 0.0
    right_hairpin_th: float = 0.0
    left_stability: float = 0.0
    right_stability: float = 0.0


def design_primers(
    template: str,
    target_start: Optional[int] = None,
    target_length: Optional[int] = None,
    product_size_min: int = DEFAULT_PRIMER_PRODUCT_MIN,
    product_size_max: int = DEFAULT_PRIMER_PRODUCT_MAX,
    primer_opt_size: int = DEFAULT_PRIMER_OPT_SIZE,
    primer_min_size: int = DEFAULT_PRIMER_MIN_SIZE,
    primer_max_size: int = DEFAULT_PRIMER_MAX_SIZE,
    primer_opt_tm: float = DEFAULT_PRIMER_OPT_TM,
    primer_min_tm: float = DEFAULT_PRIMER_MIN_TM,
    primer_max_tm: float = DEFAULT_PRIMER_MAX_TM,
    primer_max_diff_tm: float = DEFAULT_PRIMER_MAX_DIFF_TM,
    primer_gc_min: float = DEFAULT_PRIMER_GC_MIN,
    primer_gc_max: float = DEFAULT_PRIMER_GC_MAX,
    max_poly_x: int = DEFAULT_MAX_POLY_X,
    gc_clamp: int = DEFAULT_GC_CLAMP,
    max_end_stability: float = DEFAULT_MAX_END_STABILITY,
    max_end_gc: int = DEFAULT_MAX_END_GC,
    mono_cations: float = DEFAULT_MONO_CATIONS,
    diva_cations: float = DEFAULT_DIVA_CATIONS,
    con_dntps: float = DEFAULT_CON_DNTPS,
    salt_correction: int = DEFAULT_SALT_CORRECTION,
    tm_method: int = DEFAULT_TM_METHOD,
    num_return: int = DEFAULT_NUM_RETURN,
    left_primer_input: Optional[str] = None,
    right_primer_input: Optional[str] = None,
    excluded_regions: Optional[List[Tuple[int, int]]] = None,
) -> List[PrimerPair]:
    """Design PCR primers using Primer3.

    Args:
        template: Template DNA sequence (uppercase ACGT).
        target_start: 0-based start of target region within template.
        target_length: Length of target region.
        product_size_min/max: Allowed PCR product size range.
        primer_opt_size: Optimal primer length.
        primer_min_size/max_size: Primer length range.
        primer_opt_tm: Optimal melting temperature.
        primer_min_tm/max_tm: Tm range.
        primer_max_diff_tm: Max Tm difference between primers in a pair.
        primer_gc_min/max: GC content range (percent).
        max_poly_x: Max mononucleotide repeat length.
        gc_clamp: Number of consecutive G/C at 3' end.
        max_end_stability: Max stability of last five 3' bases.
        max_end_gc: Max G/C in last five 3' bases.
        mono_cations: mM monovalent cation concentration.
        diva_cations: mM divalent cation concentration.
        con_dntps: mM dNTP concentration.
        salt_correction: Salt correction formula (0=Schildkraut, 1=SantaLucia, 2=Owczarzy).
        tm_method: Tm method (0=Breslauer, 1=SantaLucia).
        num_return: Number of primer pairs to return.
        left_primer_input: User-specified forward primer sequence.
        right_primer_input: User-specified reverse primer sequence.
        excluded_regions: List of (start, length) regions to exclude.

    Returns:
        List of PrimerPair objects. Empty list if design fails.
    """
    seq_args: Dict[str, Any] = {
        "SEQUENCE_TEMPLATE": template.upper(),
    }

    if target_start is not None and target_length is not None:
        seq_args["SEQUENCE_TARGET"] = [target_start, target_length]

    if left_primer_input:
        seq_args["SEQUENCE_PRIMER"] = left_primer_input.upper()
    if right_primer_input:
        seq_args["SEQUENCE_PRIMER_REVCOMP"] = right_primer_input.upper()

    if excluded_regions:
        excl_str = " ".join(f"{s},{l}" for s, l in excluded_regions)
        seq_args["SEQUENCE_EXCLUDED_REGION"] = excl_str

    # When a target region is specified, ensure the minimum product size
    # can physically accommodate the target + flanking primers.
    # Do NOT auto-expand max_product — the user's product_size_max is a hard constraint.
    min_product = product_size_min
    max_product = product_size_max
    if target_start is not None and target_length is not None:
        min_needed = target_length + 2 * primer_min_size
        if min_product < min_needed:
            min_product = min_needed
        # Ensure valid Primer3 input (min must be <= max)
        if max_product < min_product:
            max_product = min_product

    global_args: Dict[str, Any] = {
        "PRIMER_TASK": "generic",
        "PRIMER_PICK_LEFT_PRIMER": 1,
        "PRIMER_PICK_RIGHT_PRIMER": 1,
        "PRIMER_PICK_INTERNAL_OLIGO": 0,
        "PRIMER_NUM_RETURN": num_return,
        "PRIMER_OPT_SIZE": primer_opt_size,
        "PRIMER_MIN_SIZE": primer_min_size,
        "PRIMER_MAX_SIZE": primer_max_size,
        "PRIMER_OPT_TM": primer_opt_tm,
        "PRIMER_MIN_TM": primer_min_tm,
        "PRIMER_MAX_TM": primer_max_tm,
        "PRIMER_MAX_DIFF_TM": primer_max_diff_tm,
        "PRIMER_PRODUCT_SIZE_RANGE": [[min_product, max_product]],
        "PRIMER_MIN_GC": primer_gc_min,
        "PRIMER_MAX_GC": primer_gc_max,
        "PRIMER_MAX_POLY_X": max_poly_x,
        "PRIMER_GC_CLAMP": gc_clamp,
        "PRIMER_MAX_END_STABILITY": max_end_stability,
        "PRIMER_MAX_END_GC": max_end_gc,
        "PRIMER_SALT_MONOVALENT": mono_cations,
        "PRIMER_SALT_DIVALENT": diva_cations,
        "PRIMER_DNTP_CONC": con_dntps,
        "PRIMER_SALT_CORRECTIONS": salt_correction,
        "PRIMER_TM_FORMULA": tm_method,
        "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT": 1,
    }

    tpl_len = len(template)
    if tpl_len < product_size_min:
        return []

    try:
        result = primer3.design_primers(seq_args, global_args)
    except OSError:
        return []

    pairs = _parse_primer3_result(result, template)
    return pairs


def _parse_primer3_result(result: Dict[str, Any], template: str) -> List[PrimerPair]:
    """Parse Primer3 Boulder-IO output into PrimerPair objects."""
    num_returned = result.get("PRIMER_PAIR_NUM_RETURNED", 0)
    if num_returned == 0:
        return []

    pairs = []
    for i in range(num_returned):
        left_seq = result.get(f"PRIMER_LEFT_{i}_SEQUENCE", "")
        right_seq = result.get(f"PRIMER_RIGHT_{i}_SEQUENCE", "")

        if not left_seq or not right_seq:
            continue

        pair = PrimerPair(
            left_primer=left_seq,
            right_primer=right_seq,
            left_tm=float(result.get(f"PRIMER_LEFT_{i}_TM", 0)),
            right_tm=float(result.get(f"PRIMER_RIGHT_{i}_TM", 0)),
            left_gc=float(result.get(f"PRIMER_LEFT_{i}_GC_PERCENT", 0)),
            right_gc=float(result.get(f"PRIMER_RIGHT_{i}_GC_PERCENT", 0)),
            product_size=int(result.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE", 0)),
            left_position=int(result.get(f"PRIMER_LEFT_{i}", [0, 0])[0]),
            right_position=int(result.get(f"PRIMER_RIGHT_{i}", [0, 0])[0]),
            penalty=float(result.get(f"PRIMER_PAIR_{i}_PENALTY", 999)),
            pair_product_tm=float(result.get(f"PRIMER_PAIR_{i}_PRODUCT_TM", 0)),
            left_self_any_th=float(result.get(f"PRIMER_LEFT_{i}_SELF_ANY_TH", 0)),
            left_self_end_th=float(result.get(f"PRIMER_LEFT_{i}_SELF_END_TH", 0)),
            right_self_any_th=float(result.get(f"PRIMER_RIGHT_{i}_SELF_ANY_TH", 0)),
            right_self_end_th=float(result.get(f"PRIMER_RIGHT_{i}_SELF_END_TH", 0)),
            pair_compl_any_th=float(result.get(f"PRIMER_PAIR_{i}_COMPL_ANY_TH", 0)),
            pair_compl_end_th=float(result.get(f"PRIMER_PAIR_{i}_COMPL_END_TH", 0)),
            left_hairpin_th=float(result.get(f"PRIMER_LEFT_{i}_HAIRPIN_TH", 0)),
            right_hairpin_th=float(result.get(f"PRIMER_RIGHT_{i}_HAIRPIN_TH", 0)),
            left_stability=float(result.get(f"PRIMER_LEFT_{i}_END_STABILITY", 0)),
            right_stability=float(result.get(f"PRIMER_RIGHT_{i}_END_STABILITY", 0)),
        )
        pairs.append(pair)

    return pairs


def run_thermo_analysis(primer_seq: str) -> Dict[str, Any]:
    """Run supplementary thermodynamic analysis on a single primer."""
    results = {}
    try:
        results["hairpin"] = primer3.calc_hairpin(primer_seq)
    except Exception:
        pass
    try:
        results["homodimer"] = primer3.calc_homodimer(primer_seq)
    except Exception:
        pass
    return results


def enrich_pair_thermo(pair: PrimerPair) -> PrimerPair:
    """Run thermo analysis and update the pair."""
    try:
        hetero = primer3.calc_heterodimer(pair.left_primer, pair.right_primer)
        pair.pair_compl_any_th = hetero.tm if hasattr(hetero, "tm") else 0
    except Exception:
        pass
    return pair
