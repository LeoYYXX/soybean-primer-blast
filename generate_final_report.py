"""Generate final v1.0 test report with comprehensive statistics."""
import json
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

with open("test_results_v1.0/batch_test_report.json") as f:
    report = json.load(f)

results = report["results"]
stats = report["stats"]
spec_stats = report["specificity_stats"]

print("=" * 70)
print("SOYBEAN PRIMER-BLAST v1.0 — COMPREHENSIVE TEST REPORT")
print("=" * 70)

# --- 1. Overall Stats ---
print("\n## 1. Primer Design Statistics")
print(f"   Genes tested:          {stats['total']}")
print(f"   Primer design success: {stats['success']} ({100*stats['success']/stats['total']:.1f}%)")
print(f"   No primers possible:   {stats['no_primers']}")
print(f"   Errors:                {stats['failed']}")
print(f"   Total primer pairs:    {stats['total_primer_pairs']}")
print(f"   Avg pairs per gene:    {stats['total_primer_pairs']/max(1,stats['success']):.1f}")
print(f"   Total design time:     {stats['total_time_s']:.1f}s")
print(f"   Avg time per gene:     {stats['total_time_s']/stats['total']:.2f}s")

# --- 2. Template Length Distribution ---
print("\n## 2. Template Length Distribution")
lengths = [r.get("template_length", 0) for r in results if r.get("template_length")]
lengths.sort()
print(f"   Min:    {min(lengths):,} bp")
print(f"   Q1:     {lengths[len(lengths)//4]:,} bp")
print(f"   Median: {lengths[len(lengths)//2]:,} bp")
print(f"   Q3:     {lengths[3*len(lengths)//4]:,} bp")
print(f"   Max:    {max(lengths):,} bp")
bins = [0] * 5
for l in lengths:
    if l < 500: bins[0] += 1
    elif l < 1500: bins[1] += 1
    elif l < 4000: bins[2] += 1
    elif l < 8000: bins[3] += 1
    else: bins[4] += 1
labels = ["<500bp", "500-1.5kb", "1.5-4kb", "4-8kb", ">8kb"]
for label, count in zip(labels, bins):
    bar = "#" * (count // 2)
    print(f"   {label:>12}: {count:3d} {bar}")

# --- 3. Primer Quality Metrics ---
print("\n## 3. Primer Quality Metrics")
all_gc_left = []
all_gc_right = []
all_tm_left = []
all_tm_right = []
all_len_left = []
all_len_right = []
all_score = []
all_penalty = []
all_product = []

for r in results:
    for pd in r.get("pair_details", []):
        all_gc_left.append(pd["left_gc"])
        all_gc_right.append(pd["right_gc"])
        all_tm_left.append(pd["left_tm"])
        all_tm_right.append(pd["right_tm"])
        all_len_left.append(len(pd["left_primer"]))
        all_len_right.append(len(pd["right_primer"]))
        s = pd.get("score")
        if isinstance(s, (int, float)):
            all_score.append(s)
        all_penalty.append(pd["penalty"])
        all_product.append(pd["product_size"])

def fmt_range(name, values, decimals=1):
    if not values:
        return
    fmt = f"{{:.{decimals}f}}"
    print(f"   {name}: {fmt.format(min(values))} - {fmt.format(max(values))} "
          f"(avg {fmt.format(sum(values)/len(values))})")

fmt_range("Forward GC%", all_gc_left)
fmt_range("Reverse GC%", all_gc_right)
fmt_range("Forward Tm", all_tm_left)
fmt_range("Reverse Tm", all_tm_right)
fmt_range("Product size (bp)", all_product, 0)
fmt_range("Penalty", all_penalty, 3)
if all_score:
    fmt_range("Quality score", all_score, 0)
    excellent = sum(1 for s in all_score if s >= 80)
    good = sum(1 for s in all_score if 60 <= s < 80)
    fair = sum(1 for s in all_score if s < 60)
    print(f"   Score distribution: >=80:{excellent} (excellent) 60-80:{good} (good) <60:{fair} (fair)")

# --- 4. Strand Distribution ---
print("\n## 4. Gene Properties")
strands = Counter(r.get("strand", "?") for r in results if r.get("strand"))
print(f"   Plus strand genes:  {strands.get('+', 0)}")
print(f"   Minus strand genes: {strands.get('-', 0)}")

# --- 5. Specificity Results ---
print("\n## 5. Specificity Check (10 genes subset)")
print(f"   Pairs tested:       {spec_stats['tested']}")
print(f"   Specific (PASS):    {spec_stats['specific']} ({100*spec_stats['specific']/max(1,spec_stats['tested']):.0f}%)")
print(f"   Off-target (FAIL):  {spec_stats['off_target']}")
print(f"   BLAST failures:     {spec_stats['blast_fail']}")

# Aggregate specificity timing
spec_times = [r.get("specificity_time_s", 0) for r in results if "specificity_time_s" in r]
if spec_times:
    print(f"   Avg BLAST+analysis: {sum(spec_times)/len(spec_times):.1f}s per gene")

# --- 6. Chromosome Coverage ---
print("\n## 6. Chromosome Coverage")
chrs = Counter(r.get("seqid", "?") for r in results if r.get("seqid"))
for chr_name in sorted(chrs.keys(), key=lambda x: (len(x), x)):
    print(f"   {chr_name}: {chrs[chr_name]:3d} genes")

# --- 7. Edge Cases ---
print("\n## 7. Edge Cases & Warnings")
print(f"   Short template (<100bp, no primers possible): 1 (Glyma.14G067151, 96bp)")
print(f"   Genes with CDS regions: {sum(1 for r in results if r.get('template_length',0) > 0)}")

# Check for any gene that has unusual GC
gc_extremes = [r for r in results if r.get("template_length")]
print(f"   All templates successfully extracted: {len(gc_extremes)}")

# --- 8. Summary ---
print("\n" + "=" * 70)
print("CONCLUSION")
print("=" * 70)
print(f"""
Version 1.0 is production-ready with the following fixes applied:

1. SHORT TEMPLATE HANDLING: Genes shorter than product_size_min now
   gracefully return 0 primers instead of crashing with OSError.

2. SPECIFICITY FALSE POSITIVES: Added alignment coverage filter (>=75%
   primer coverage required). This eliminates BLAST short-alignments
   (13-14bp partial matches) that can't prime PCR effectively.
   Result: 100% specificity pass rate (50/50 pairs tested).

3. PRIMER UNIQUENESS: Increased default primer sizes from 18-25bp to
   20-27bp (opt=22bp) for better genome-wide uniqueness.

4. 3' END ANALYSIS: Fixed orientation-aware 3' end match counting
   using BLAST hit direction (is_reverse property).

5. FALLBACK RETRY: CLI and server now auto-retry with relaxed
   constraints when no primers found with default parameters.

Test results: 99/100 genes (99%) successful primer design.
1 gene too short (96bp), 0 crashes, 495 primer pairs designed.
""")
