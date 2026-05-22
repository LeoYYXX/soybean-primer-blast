"""Batch test: run primer design on ~100 soybean genes, collect results and errors."""
import json
import os
import sys
import time
import traceback
from collections import defaultdict
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from primer_blast.fasta_index import load_fasta_index
from primer_blast.gff_index import load_or_build_index
from primer_blast.sequence import extract_gene_sequence
from primer_blast.primer_design import design_primers
from primer_blast.primer_score import score_primer_pair
from primer_blast.constants import (
    DEFAULT_GENOME_FA, DEFAULT_GFF3,
)

# --- Config ---
GENES_TO_TEST = 1000
SPECIFICITY_SUBSET = 20  # number of genes for full specificity test
OUTPUT_DIR = "test_results_v1.0_1000genes"

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=" * 70)
    print("Soybean Primer-BLAST v1.0 — Batch Test (1000 genes)")
    print(f"Started: {datetime.now().isoformat()}")
    print("=" * 70)

    # --- Load indices ---
    print("\n[1/4] Loading genome indices...")
    t0 = time.time()
    if not os.path.exists(DEFAULT_GENOME_FA):
        print(f"FATAL: Genome FASTA not found: {DEFAULT_GENOME_FA}")
        sys.exit(1)
    fai = load_fasta_index(DEFAULT_GENOME_FA)
    gene_index, transcript_index = load_or_build_index(DEFAULT_GFF3)
    print(f"  FASTA: {len(fai)} sequences")
    print(f"  Genes: {len(gene_index)}, Transcripts: {len(transcript_index)}")
    print(f"  Time: {time.time() - t0:.1f}s")

    # --- Select genes ---
    print(f"\n[2/4] Selecting {GENES_TO_TEST} genes...")
    all_gene_names = sorted(transcript_index.keys())
    step = max(1, len(all_gene_names) // GENES_TO_TEST)
    test_genes = [all_gene_names[i] for i in range(0, len(all_gene_names), step)][:GENES_TO_TEST]
    print(f"  Selected {len(test_genes)} genes across {len(all_gene_names)} total")

    # --- Run primer design ---
    print(f"\n[3/4] Running primer design on {len(test_genes)} genes...")
    results = []
    errors = defaultdict(list)
    stats = {
        "total": len(test_genes),
        "success": 0,
        "failed": 0,
        "no_primers": 0,
        "total_primer_pairs": 0,
        "total_time_s": 0,
    }

    t_start = time.time()
    for i, gene_name in enumerate(test_genes):
        gene_start = time.time()
        transcript = transcript_index.get(gene_name)
        status = {"gene": gene_name, "error": None, "pairs": 0}

        try:
            template = extract_gene_sequence(DEFAULT_GENOME_FA, transcript, fai)
            if not template:
                status["error"] = "empty_template"
                errors["empty_template"].append(gene_name)
                stats["failed"] += 1
                results.append(status)
                continue

            status["template_length"] = len(template)
            status["strand"] = transcript.strand
            status["seqid"] = transcript.seqid

            pairs = design_primers(
                template=template,
                target_start=None,
                target_length=None,
                product_size_min=100,
                product_size_max=400,
                num_return=5,
            )

            if not pairs:
                status["error"] = "no_primers"
                errors["no_primers"].append(gene_name)
                stats["no_primers"] += 1
                results.append(status)
                continue

            status["pairs"] = len(pairs)
            stats["total_primer_pairs"] += len(pairs)
            stats["success"] += 1

            # Score each pair
            pair_details = []
            for p in pairs:
                detail = {
                    "left_primer": p.left_primer,
                    "right_primer": p.right_primer,
                    "left_tm": round(p.left_tm, 1),
                    "right_tm": round(p.right_tm, 1),
                    "left_gc": round(p.left_gc, 1),
                    "right_gc": round(p.right_gc, 1),
                    "product_size": p.product_size,
                    "penalty": round(p.penalty, 3),
                    "left_position": p.left_position,
                    "right_position": p.right_position,
                }
                try:
                    score_result = score_primer_pair(p)
                    detail["score"] = score_result.get("total") if isinstance(score_result, dict) else score_result
                except Exception:
                    detail["score"] = None
                pair_details.append(detail)

            status["pair_details"] = pair_details

        except Exception as e:
            status["error"] = type(e).__name__
            error_key = f"{type(e).__name__}: {str(e)[:80]}"
            errors[error_key].append(gene_name)
            stats["failed"] += 1

        elapsed = time.time() - gene_start
        status["elapsed_s"] = round(elapsed, 2)
        results.append(status)

        if (i + 1) % 100 == 0:
            print(f"  Progress: {i+1}/{len(test_genes)} genes "
                  f"({stats['success']} ok, {stats['failed']} fail, "
                  f"{stats['no_primers']} no-primers)")

    stats["total_time_s"] = round(time.time() - t_start, 1)
    print(f"  Done. Total time: {stats['total_time_s']:.1f}s")
    print(f"  Success: {stats['success']}, Failed: {stats['failed']}, "
          f"No-primers: {stats['no_primers']}")

    # --- Specificity test on subset ---
    print(f"\n[4/4] Specificity check on {SPECIFICITY_SUBSET} genes...")
    specificity_results = []
    spec_stats = {"tested": 0, "specific": 0, "off_target": 0, "blast_fail": 0}

    # Pick first SPECIFICITY_SUBSET successful genes for specificity test
    spec_genes = [r for r in results if not r["error"]][:SPECIFICITY_SUBSET]

    if spec_genes:
        from primer_blast.blast_db import get_blast_db
        from primer_blast.specificity import run_blast_primers, evaluate_specificity

        blast_db = get_blast_db(
            fasta_path=DEFAULT_GENOME_FA,
            db_dir=".primer_blast_cache/blastdb",
            force_rebuild=False,
        )

        for r in spec_genes:
            gene_name = r["gene"]
            transcript = transcript_index[gene_name]

            # Get all primers from designed pairs
            all_primers = set()
            for pd in r.get("pair_details", []):
                all_primers.add(pd["left_primer"])
                all_primers.add(pd["right_primer"])

            if not all_primers:
                continue

            try:
                t0_spec = time.time()
                hits = run_blast_primers(
                    primer_sequences=list(all_primers),
                    blast_db=blast_db,
                )

                spec_pairs = []
                for pd in r.get("pair_details", []):
                    left_hits = hits.get(pd["left_primer"], [])
                    right_hits = hits.get(pd["right_primer"], [])

                    # Genomic coords for left/right primers
                    if transcript.strand == "+":
                        left_gen_start = transcript.start + pd["left_position"]
                        left_gen_end = left_gen_start + len(pd["left_primer"]) - 1
                        right_gen_end = transcript.start + pd["right_position"]
                        right_gen_start = right_gen_end - len(pd["right_primer"]) + 1
                    else:
                        left_gen_end = transcript.end - pd["left_position"]
                        left_gen_start = left_gen_end - len(pd["left_primer"]) + 1
                        right_gen_start = transcript.end - pd["right_position"]
                        right_gen_end = right_gen_start + len(pd["right_primer"]) - 1

                    result = evaluate_specificity(
                        left_hits=left_hits,
                        right_hits=right_hits,
                        target_seqid=transcript.seqid,
                        target_left_coords=(min(left_gen_start, left_gen_end), max(left_gen_start, left_gen_end)),
                        target_right_coords=(min(right_gen_start, right_gen_end), max(right_gen_start, right_gen_end)),
                    )
                    spec_pairs.append({
                        "is_specific": result.is_specific,
                        "off_target_count": result.off_target_count,
                    })

                    if result.is_specific:
                        spec_stats["specific"] += 1
                    else:
                        spec_stats["off_target"] += 1

                spec_stats["tested"] += len(spec_pairs)
                r["specificity_pairs"] = spec_pairs
                r["specificity_time_s"] = round(time.time() - t0_spec, 1)

            except Exception as e:
                r["specificity_error"] = f"{type(e).__name__}: {str(e)[:100]}"
                spec_stats["blast_fail"] += 1

        print(f"  Specificity tested: {spec_stats['tested']} primer pairs "
              f"({spec_stats['specific']} specific, {spec_stats['off_target']} off-target, "
              f"{spec_stats['blast_fail']} BLAST failures)")

    # --- Save detailed results ---
    print(f"\n{'=' * 70}")
    print("Saving results...")

    report = {
        "version": "1.0",
        "timestamp": datetime.now().isoformat(),
        "config": {
            "genome": DEFAULT_GENOME_FA,
            "gff3": DEFAULT_GFF3,
            "genes_tested": GENES_TO_TEST,
            "specificity_subset": SPECIFICITY_SUBSET,
        },
        "stats": stats,
        "specificity_stats": spec_stats,
        "errors": {k: v for k, v in errors.items()},
        "error_summary": {k: len(v) for k, v in errors.items()},
        "results": results,
    }

    report_path = os.path.join(OUTPUT_DIR, "batch_test_report.json")
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    print(f"  Full report: {report_path}")

    # --- Summary ---
    summary_lines = [
        "",
        "=" * 70,
        "BATCH TEST SUMMARY",
        "=" * 70,
        f"Genes tested:       {stats['total']}",
        f"Primer design OK:   {stats['success']} ({100*stats['success']/max(1,stats['total']):.1f}%)",
        f"  No primers found: {stats['no_primers']}",
        f"  Errors:           {stats['failed']}",
        f"Total primer pairs: {stats['total_primer_pairs']}",
        f"Avg pairs/gene:     {stats['total_primer_pairs']/max(1,stats['success']):.1f}",
        f"Total time:         {stats['total_time_s']:.1f}s",
        f"Avg time/gene:      {stats['total_time_s']/stats['total']:.2f}s",
        "",
        f"Specificity subset ({spec_stats['tested']} pairs):",
        f"  Specific:   {spec_stats['specific']}",
        f"  Off-target: {spec_stats['off_target']}",
        f"  BLAST fail: {spec_stats['blast_fail']}",
        "",
        "ERROR BREAKDOWN:",
    ]

    for error_type, genes in sorted(errors.items(), key=lambda x: -len(x[1])):
        summary_lines.append(f"  [{len(genes):3d}] {error_type}")
        if len(genes) <= 5:
            summary_lines.append(f"         Genes: {', '.join(genes)}")

    summary_lines.append("")
    summary_lines.append("=" * 70)

    summary_text = "\n".join(summary_lines)
    print(summary_text)

    summary_path = os.path.join(OUTPUT_DIR, "summary.txt")
    with open(summary_path, "w") as f:
        f.write(summary_text)
    print(f"\nSummary saved to: {summary_path}")

    return report


if __name__ == "__main__":
    main()
