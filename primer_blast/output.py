"""Output formatters for primer design results."""

import json
from typing import List, Optional

from .primer_design import PrimerPair
from .specificity import SpecificityResult


def format_text(
    pairs: List[PrimerPair],
    specificity_results: Optional[List[SpecificityResult]] = None,
    template_info: str = "",
) -> str:
    """Format primer pairs as a human-readable text table."""
    lines = []

    if template_info:
        lines.append(f"Template: {template_info}")
        lines.append("=" * 60)

    for i, pair in enumerate(pairs):
        lines.append(f"")
        lines.append(f"=== Primer Pair {i + 1} ===")
        lines.append(f"Forward:  {pair.left_primer}  ({len(pair.left_primer)} nt)")
        lines.append(f"          Tm={pair.left_tm:.1f}C  GC={pair.left_gc:.1f}%  Pos={pair.left_position}")
        lines.append(f"Reverse:  {pair.right_primer}  ({len(pair.right_primer)} nt)")
        lines.append(f"          Tm={pair.right_tm:.1f}C  GC={pair.right_gc:.1f}%  Pos={pair.right_position}")
        lines.append(f"DNA amplicon (genomic):  {pair.product_size} bp  Tm={pair.pair_product_tm:.1f}C  Penalty={pair.penalty:.3f}")

        # Self-complementarity
        if any([
            pair.left_self_any_th, pair.left_self_end_th,
            pair.right_self_any_th, pair.right_self_end_th,
            pair.left_hairpin_th, pair.right_hairpin_th,
        ]):
            lines.append(f"  L self: any={pair.left_self_any_th:.1f} end={pair.left_self_end_th:.1f} hp={pair.left_hairpin_th:.1f}")
            lines.append(f"  R self: any={pair.right_self_any_th:.1f} end={pair.right_self_end_th:.1f} hp={pair.right_hairpin_th:.1f}")

        # Specificity
        if specificity_results and i < len(specificity_results):
            sr = specificity_results[i]
            if sr.is_specific:
                lines.append(f"Specific: PASS")
            else:
                lines.append(f"Specific: FAIL ({sr.off_target_count} off-target amplicons)")
                for amplicon in sr.off_target_amplicons[:5]:
                    gene_label = f" gene={amplicon.gene_name}" if amplicon.gene_name else ""
                    lines.append(
                        f"  Off-target: {amplicon.chromosome} "
                        f"({amplicon.left_hit.subject_start}-{amplicon.right_hit.subject_end}) "
                        f"size={amplicon.product_size}bp{gene_label}"
                    )
            if sr.left_target_hit:
                lines.append(
                    f"  L target: {sr.left_target_hit.subject_id}:"
                    f"{sr.left_target_hit.subject_start}-{sr.left_target_hit.subject_end}"
                    f" mm={sr.left_target_hit.num_mismatch}"
                )
            if sr.right_target_hit:
                lines.append(
                    f"  R target: {sr.right_target_hit.subject_id}:"
                    f"{sr.right_target_hit.subject_start}-{sr.right_target_hit.subject_end}"
                    f" mm={sr.right_target_hit.num_mismatch}"
                )

        lines.append("-" * 60)

    return "\n".join(lines)


def format_json(
    pairs: List[PrimerPair],
    specificity_results: Optional[List[SpecificityResult]] = None,
    template_info: str = "",
) -> str:
    """Format results as structured JSON."""
    output = {
        "template": template_info,
        "pairs": [],
    }

    for i, pair in enumerate(pairs):
        pair_dict = {
            "index": i + 1,
            "forward_primer": pair.left_primer,
            "reverse_primer": pair.right_primer,
            "forward_length": len(pair.left_primer),
            "reverse_length": len(pair.right_primer),
            "forward_tm": round(pair.left_tm, 1),
            "reverse_tm": round(pair.right_tm, 1),
            "forward_gc_percent": round(pair.left_gc, 1),
            "reverse_gc_percent": round(pair.right_gc, 1),
            "product_size": pair.product_size,
            "product_tm": round(pair.pair_product_tm, 1),
            "penalty": round(pair.penalty, 3),
            "forward_position": pair.left_position,
            "reverse_position": pair.right_position,
            "self_complementarity": {
                "forward_self_any": round(pair.left_self_any_th, 1),
                "forward_self_end": round(pair.left_self_end_th, 1),
                "forward_hairpin": round(pair.left_hairpin_th, 1),
                "reverse_self_any": round(pair.right_self_any_th, 1),
                "reverse_self_end": round(pair.right_self_end_th, 1),
                "reverse_hairpin": round(pair.right_hairpin_th, 1),
            },
        }

        if specificity_results and i < len(specificity_results):
            sr = specificity_results[i]
            pair_dict["specificity"] = {
                "is_specific": sr.is_specific,
                "off_target_count": sr.off_target_count,
                "off_target_amplicons": [
                    {
                        "chromosome": a.chromosome,
                        "gene_name": a.gene_name or "",
                        "left_position": a.left_hit.subject_start,
                        "right_position": a.right_hit.subject_end,
                        "product_size": a.product_size,
                    }
                    for a in sr.off_target_amplicons
                ],
            }
            if sr.left_target_hit:
                pair_dict["specificity"]["left_target"] = {
                    "chromosome": sr.left_target_hit.subject_id,
                    "start": sr.left_target_hit.subject_start,
                    "end": sr.left_target_hit.subject_end,
                    "mismatches": sr.left_target_hit.num_mismatch,
                }
            if sr.right_target_hit:
                pair_dict["specificity"]["right_target"] = {
                    "chromosome": sr.right_target_hit.subject_id,
                    "start": sr.right_target_hit.subject_start,
                    "end": sr.right_target_hit.subject_end,
                    "mismatches": sr.right_target_hit.num_mismatch,
                }

        output["pairs"].append(pair_dict)

    return json.dumps(output, indent=2, ensure_ascii=False)


def format_tsv(
    pairs: List[PrimerPair],
    specificity_results: Optional[List[SpecificityResult]] = None,
    template_info: str = "",
) -> str:
    """Format results as tab-separated values."""
    header = [
        "Pair", "Forward", "Forward_len", "Forward_Tm", "Forward_GC",
        "Reverse", "Reverse_len", "Reverse_Tm", "Reverse_GC",
        "Product_bp", "Product_Tm", "Penalty",
        "L_self_any", "L_self_end", "L_hairpin",
        "R_self_any", "R_self_end", "R_hairpin",
    ]
    if specificity_results:
        header.append("Specific")
        header.append("Off_targets")

    lines = ["\t".join(header)]

    for i, pair in enumerate(pairs):
        row = [
            str(i + 1),
            pair.left_primer, str(len(pair.left_primer)),
            f"{pair.left_tm:.1f}", f"{pair.left_gc:.1f}",
            pair.right_primer, str(len(pair.right_primer)),
            f"{pair.right_tm:.1f}", f"{pair.right_gc:.1f}",
            str(pair.product_size), f"{pair.pair_product_tm:.1f}",
            f"{pair.penalty:.3f}",
            f"{pair.left_self_any_th:.1f}", f"{pair.left_self_end_th:.1f}",
            f"{pair.left_hairpin_th:.1f}",
            f"{pair.right_self_any_th:.1f}", f"{pair.right_self_end_th:.1f}",
            f"{pair.right_hairpin_th:.1f}",
        ]
        if specificity_results and i < len(specificity_results):
            sr = specificity_results[i]
            row.append("PASS" if sr.is_specific else "FAIL")
            row.append(str(sr.off_target_count))

        lines.append("\t".join(row))

    return "\n".join(lines)
