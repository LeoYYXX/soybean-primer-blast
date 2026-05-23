"""Main entry point: orchestrates primer design and specificity checking."""

import os
import sys

from .annotation import get_gene_annotation, load_annotations
from .cli import build_parser, parse_region
from .fasta_index import load_fasta_index
from .genbank import export_genbank
from .gff_index import load_or_build_index, build_genomic_gene_map, TranscriptRecord, GeneRecord
from .sequence import (
    extract_gene_sequence,
    extract_region_with_padding,
    extract_cds_sequence,
    genomic_to_template_coords,
    get_target_coordinates,
)
from .primer_design import design_primers, PrimerPair
from .blast_db import get_blast_db
from .constants import DEFAULT_TRANSCRIPT_FA, OFF_TARGET_PRESETS
from .specificity import (
    run_blast_primers,
    evaluate_specificity,
    find_target_hit,
    SpecificityResult,
)
from .output import format_text, format_json, format_tsv
from .primer_check import check_primers
from .server import run_server
from .transcript_seq import (
    build_transcript_seq_index,
    compute_cdna_amplicon,
)


def _extract_template_for_gene(args, transcript, fai):
    """Extract template and determine target region for a gene-based query."""
    genome = args.genome

    if args.target in ("cds", "utr"):
        # Use full gene body as template, target the specific region
        template = extract_gene_sequence(genome, transcript, fai)
        target_coords = get_target_coordinates(transcript, args.target)
        if target_coords is None:
            print(f"Warning: No {args.target} region found for gene. Targeting full gene body.")
            return template, None, None

        tgt_start, tgt_len = genomic_to_template_coords(
            target_coords[0], target_coords[1],
            transcript.start, transcript.end, transcript.strand,
        )
        info = f"{transcript.seqid}:{transcript.start}-{transcript.end}({transcript.strand}) {args.target}"
        return template, tgt_start, tgt_len, info

    elif args.target == "gene" or args.target == "all":
        template = extract_gene_sequence(genome, transcript, fai)
        info = f"{transcript.seqid}:{transcript.start}-{transcript.end}({transcript.strand}) gene"
        return template, None, None, info

    else:
        # User-specified coordinates
        target_coords = get_target_coordinates(transcript, args.target)
        if target_coords:
            template = extract_gene_sequence(genome, transcript, fai)
            tgt_start, tgt_len = genomic_to_template_coords(
                target_coords[0], target_coords[1],
                transcript.start, transcript.end, transcript.strand,
            )
            info = f"{transcript.seqid}:{transcript.start}-{transcript.end}({transcript.strand}) target={args.target}"
            return template, tgt_start, tgt_len, info

    template = extract_gene_sequence(genome, transcript, fai)
    info = f"{transcript.seqid}:{transcript.start}-{transcript.end}({transcript.strand}) gene"
    return template, None, None, info


def _guess_primer_genomic_coords(
    pair: PrimerPair,
    template_start: int,
    template_end: int,
    strand: str,
) -> tuple:
    """Convert primer template-relative positions to genomic coordinates."""
    if strand == "+":
        left_gen_start = template_start + pair.left_position
        left_gen_end = left_gen_start + len(pair.left_primer) - 1
        right_gen_end = template_start + pair.right_position
        right_gen_start = right_gen_end - len(pair.right_primer) + 1
    else:
        # Minus strand: template is reverse-complemented
        left_gen_end = template_end - pair.left_position
        left_gen_start = left_gen_end - len(pair.left_primer) + 1
        right_gen_start = template_end - pair.right_position
        right_gen_end = right_gen_start + len(pair.right_primer) - 1

    return (
        (min(left_gen_start, left_gen_end), max(left_gen_start, left_gen_end)),
        (min(right_gen_start, right_gen_end), max(right_gen_start, right_gen_end)),
    )


def run(args):
    """Run the primer design and specificity checking pipeline."""
    # --- Validate input mode ---
    if not args.gene and not args.region and not args.sequence and not args.check_primers:
        print("Error: One of --gene, --region, --sequence, or --check-primers is required.", file=sys.stderr)
        print("Use --help for usage information.", file=sys.stderr)
        sys.exit(1)

    if args.check_primers:
        return run_check_primers(args)

    # --- Validate parameter ranges ---
    errors = []
    if args.product_min > args.product_max:
        errors.append(f"product-min ({args.product_min}) > product-max ({args.product_max})")
    if args.primer_min_size > args.primer_max_size:
        errors.append(f"primer-min-size ({args.primer_min_size}) > primer-max-size ({args.primer_max_size})")
    if args.primer_opt_size < args.primer_min_size or args.primer_opt_size > args.primer_max_size:
        errors.append(f"primer-opt-size ({args.primer_opt_size}) outside min-max range ({args.primer_min_size}-{args.primer_max_size})")
    if args.primer_min_tm > args.primer_max_tm:
        errors.append(f"primer-min-tm ({args.primer_min_tm}) > primer-max-tm ({args.primer_max_tm})")
    if args.primer_opt_tm < args.primer_min_tm or args.primer_opt_tm > args.primer_max_tm:
        errors.append(f"primer-opt-tm ({args.primer_opt_tm}) outside min-max range ({args.primer_min_tm}-{args.primer_max_tm})")
    if args.primer_gc_min > args.primer_gc_max:
        errors.append(f"primer-gc-min ({args.primer_gc_min}) > primer-gc-max ({args.primer_gc_max})")
    if errors:
        for e in errors:
            print(f"Error: {e}", file=sys.stderr)
        print("Fix parameter ranges and retry.", file=sys.stderr)
        sys.exit(1)

    # --- Step 1: Resolve input ---
    genome_path = args.genome
    gff3_path = args.gff3

    if not os.path.exists(genome_path):
        print(f"Error: Genome FASTA not found: {genome_path}", file=sys.stderr)
        sys.exit(1)

    seqid: str = ""
    template: str = ""
    template_info: str = ""
    transcript: TranscriptRecord = None
    gene_record: GeneRecord = None
    genomic_gene_map: dict = None
    tgt_start = None
    tgt_len = None

    fai = load_fasta_index(genome_path)

    if args.gene:
        if not os.path.exists(gff3_path):
            print(f"Error: GFF3 file not found: {gff3_path}", file=sys.stderr)
            sys.exit(1)

        gene_index, transcript_index = load_or_build_index(gff3_path, args.cache_dir)
        genomic_gene_map = build_genomic_gene_map(gene_index)
        gene_record = gene_index.get(args.gene)
        transcript = transcript_index.get(args.gene)

        if not transcript:
            print(f"Error: Gene '{args.gene}' not found in GFF3.", file=sys.stderr)
            # List similar genes
            similar = [g for g in transcript_index if g.startswith(args.gene[:8])][:5]
            if similar:
                print(f"Similar genes: {', '.join(similar)}", file=sys.stderr)
            sys.exit(1)

        seqid = transcript.seqid
        result = _extract_template_for_gene(args, transcript, fai)
        template, tgt_start, tgt_len, template_info = result[:4]
        if len(result) > 4:
            template_info = result[4]

    elif args.region:
        seqid, start, end, strand = parse_region(args.region)
        template_info = f"{seqid}:{start}-{end}({strand}) region"
        template = extract_region_with_padding(
            genome_path, seqid, start, end, strand, index=fai,
        )
        region_start = start
        region_end = end
        region_strand = strand

        if args.target and args.target not in ("gene", "all"):
            tgt = get_target_coordinates(None, args.target)
            if tgt:
                # Numeric targets (e.g. "200,600") are template-relative, use directly
                tgt_start, tgt_len = tgt[0], tgt[1] - tgt[0] + 1

    elif args.sequence:
        # Could be a file path or direct sequence
        if os.path.exists(args.sequence):
            with open(args.sequence, "r") as f:
                content = f.read().strip()
            if content.startswith(">"):
                # FASTA file
                lines = content.split("\n")
                template = "".join(l.strip() for l in lines if not l.startswith(">")).upper()
                template_info = f"File: {args.sequence} ({len(template)} bp)"
            else:
                template = content.upper()
                template_info = f"File: {args.sequence} ({len(template)} bp)"
        else:
            template = args.sequence.upper()
            template_info = f"Direct input ({len(template)} bp)"

        if args.target and args.target not in ("gene", "all", "cds"):
            tgt = get_target_coordinates(None, args.target)
            if tgt:
                tgt_start, tgt_len = tgt

    if not template:
        print("Error: Could not extract template sequence.", file=sys.stderr)
        sys.exit(1)

    if args.verbose:
        print(f"Template: {template_info}")
        print(f"Template length: {len(template)} bp")
        if tgt_start is not None:
            print(f"Target region: {tgt_start}-{tgt_start + tgt_len} (0-based within template)")

    # --- Step 2: Design primers ---
    if args.verbose:
        print(f"\nDesigning primers...")

    pairs = design_primers(
        template=template,
        target_start=tgt_start,
        target_length=tgt_len,
        product_size_min=args.product_min,
        product_size_max=args.product_max,
        primer_opt_size=args.primer_opt_size,
        primer_min_size=args.primer_min_size,
        primer_max_size=args.primer_max_size,
        primer_opt_tm=args.primer_opt_tm,
        primer_min_tm=args.primer_min_tm,
        primer_max_tm=args.primer_max_tm,
        primer_gc_min=args.primer_gc_min,
        primer_gc_max=args.primer_gc_max,
        num_return=args.num_return,
    )

    if args.verbose:
        print(f"Designed {len(pairs)} primer pairs")

    # Fallback: retry with relaxed constraints
    if not pairs:
        if args.verbose:
            print("No primers with default constraints, retrying with relaxed parameters...")
        pairs = design_primers(
            template=template,
            target_start=tgt_start,
            target_length=tgt_len,
            product_size_min=min(70, args.product_min),
            product_size_max=args.product_max,
            primer_opt_size=args.primer_opt_size,
            primer_min_size=max(15, args.primer_min_size - 3),
            primer_max_size=args.primer_max_size + 3,
            primer_opt_tm=args.primer_opt_tm,
            primer_min_tm=max(50, args.primer_min_tm - 5),
            primer_max_tm=min(70, args.primer_max_tm + 5),
            primer_gc_min=max(20, args.primer_gc_min - 10),
            primer_gc_max=min(80, args.primer_gc_max + 10),
            num_return=args.num_return,
        )

    if not pairs:
        print("Error: No primer pairs could be designed with the given parameters.", file=sys.stderr)
        print("Try relaxing constraints (wider Tm/GC range, larger product size).", file=sys.stderr)
        sys.exit(1)

    # --- Step 3: Specificity checking ---
    specificity_results = None

    if not args.skip_specificity:
        if args.verbose:
            print(f"\nChecking primer specificity...")

        blast_db = args.blast_db
        if not blast_db:
            blast_db = get_blast_db(
                fasta_path=genome_path,
                db_dir=os.path.join(args.cache_dir, "blastdb"),
                makeblastdb_path=args.makeblastdb,
                force_rebuild=args.force_rebuild_blast,
            )

        # Collect all unique primer sequences
        all_primers = set()
        for p in pairs:
            all_primers.add(p.left_primer)
            all_primers.add(p.right_primer)

        # Run BLAST
        hits = run_blast_primers(
            primer_sequences=list(all_primers),
            blast_db=blast_db,
            blastn_path=args.blastn,
        )

        # Evaluate specificity for each pair
        specificity_results = []
        for pair in pairs:
            left_hits = hits.get(pair.left_primer, [])
            right_hits = hits.get(pair.right_primer, [])

            if seqid and transcript:
                # Determine genomic coordinates for the intended binding sites
                gstart = transcript.start
                gend = transcript.end
                gstrand = transcript.strand if transcript else "+"

                left_gen_coords, right_gen_coords = _guess_primer_genomic_coords(
                    pair, gstart, gend, gstrand,
                )
            elif seqid:
                # For region-based input, estimate binding site from template positions
                region_parts = args.region.rstrip(":-").rstrip(":+").split(":")
                if len(region_parts) == 2:
                    coords = region_parts[1].split("-")
                    if len(coords) == 2:
                        region_start = int(coords[0])
                        strand_from_region = "+"
                        if args.region.endswith(":-"):
                            strand_from_region = "-"
                        left_gen_coords, right_gen_coords = _guess_primer_genomic_coords(
                            pair, region_start, region_start + len(template) - 1,
                            strand_from_region,
                        )
                    else:
                        left_gen_coords = (0, len(pair.left_primer))
                        right_gen_coords = (0, len(pair.right_primer))
                else:
                    left_gen_coords = (0, len(pair.left_primer))
                    right_gen_coords = (0, len(pair.right_primer))
            else:
                # Direct sequence input - no genome context for specificity
                left_gen_coords = (0, len(pair.left_primer))
                right_gen_coords = (0, len(pair.right_primer))

            # Apply off-target sensitivity preset
            sensitivity = getattr(args, "off_target_sensitivity", "standard")
            preset = OFF_TARGET_PRESETS.get(sensitivity, OFF_TARGET_PRESETS["standard"])
            result = evaluate_specificity(
                left_hits=left_hits,
                right_hits=right_hits,
                target_seqid=seqid if seqid else "",
                target_left_coords=left_gen_coords,
                target_right_coords=right_gen_coords,
                max_mismatches=preset["max_mismatches"],
                min_3prime_match=preset["min_3prime"],
                min_alignment_cover=preset["min_coverage"],
                max_amplicon_size=args.max_amplicon_size,
                genomic_gene_map=genomic_gene_map,
            )
            specificity_results.append(result)

        # Re-rank: specific pairs first, then by penalty
        if specificity_results:
            ranked = sorted(
                zip(pairs, specificity_results),
                key=lambda x: (not x[1].is_specific, x[0].penalty),
            )
            pairs = [p for p, _ in ranked]
            specificity_results = [s for _, s in ranked]

    # --- Compute cDNA amplicons (gene mode) ---
    cdna_info_list = []
    if args.gene and transcript and os.path.exists(getattr(args, "transcript_fasta", DEFAULT_TRANSCRIPT_FA)):
        try:
            ts_index = build_transcript_seq_index(args.transcript_fasta)
            exons = transcript.all_exons or transcript.cds_exons
            for pair in pairs:
                gstrand = transcript.strand
                gstart = transcript.start
                gend = transcript.end
                if gstrand == "+":
                    left_gen_pos = gstart + pair.left_position
                    right_gen_pos = gstart + pair.right_position
                else:
                    left_gen_pos = gend - pair.left_position
                    right_gen_pos = gend - pair.right_position

                transcript_id = transcript.mrna_id
                if transcript_id not in ts_index:
                    for tid in ts_index:
                        if transcript.gene_short_name in tid:
                            transcript_id = tid
                            break

                cdna = compute_cdna_amplicon(
                    fasta_path=args.transcript_fasta,
                    transcript_id=transcript_id,
                    exons=exons,
                    strand=gstrand,
                    left_genomic_pos=left_gen_pos,
                    right_genomic_pos=right_gen_pos,
                    index=ts_index,
                )
                cdna_info_list.append(cdna)
        except Exception as e:
            if args.verbose:
                print(f"cDNA computation skipped: {e}", file=sys.stderr)

    # --- GenBank export ---
    if args.format == "genbank" or args.export_genbank:
        if not gene_record:
            print("Error: GenBank export requires --gene input.", file=sys.stderr)
            sys.exit(1)

        # Load annotation
        annotation = None
        annotation_path = getattr(args, "annotation_tsv", "")
        if annotation_path and os.path.exists(annotation_path):
            ann_idx = load_annotations(annotation_path)
            annotation = get_gene_annotation(args.gene, ann_idx)

        gb_text = export_genbank(
            gene_name=args.gene,
            gene_record=gene_record,
            transcript=transcript,
            template_sequence=template,
            pairs=pairs,
            specificity_results=specificity_results,
            cdna_info=cdna_info_list if cdna_info_list else None,
            annotation=annotation,
        )

        genbank_path = args.export_genbank or args.output
        if genbank_path == "-":
            print(gb_text)
        else:
            with open(genbank_path, "w") as f:
                f.write(gb_text)
            print(f"GenBank file written to {genbank_path}")
        return

    # --- Standard output ---
    if args.format == "json":
        output = format_json(pairs, specificity_results, template_info)
    elif args.format == "tsv":
        output = format_tsv(pairs, specificity_results, template_info)
    else:
        output = format_text(pairs, specificity_results, template_info)

    if args.output == "-":
        print(output)
    else:
        with open(args.output, "w") as f:
            f.write(output)
        print(f"Results written to {args.output}")


def run_check_primers(args):
    """Check user-provided primers for quality and specificity."""
    # Parse forward,reverse sequences
    parts = args.check_primers.split(",")
    if len(parts) != 2:
        print("Error: --check-primers requires two comma-separated sequences: forward,reverse",
              file=sys.stderr)
        sys.exit(1)
    fwd_seq, rev_seq = parts[0].strip().upper(), parts[1].strip().upper()

    if not fwd_seq or not rev_seq:
        print("Error: Both primer sequences must be non-empty.", file=sys.stderr)
        sys.exit(1)

    # Validate sequences
    valid_bases = set('ATGCRYSWKMBDHVN')
    for label, seq in [("Forward", fwd_seq), ("Reverse", rev_seq)]:
        if not all(b in valid_bases for b in seq):
            print(f"Error: {label} primer contains invalid nucleotide characters.", file=sys.stderr)
            sys.exit(1)

    # Quality check
    result = check_primers(fwd_seq, rev_seq)

    # Specificity check
    specificity_results = None
    if not args.skip_specificity:
        genome_path = args.genome
        if not os.path.exists(genome_path):
            print(f"Warning: Genome FASTA not found at {genome_path}. Skipping specificity.",
                  file=sys.stderr)
        else:
            blast_db = args.blast_db
            if not blast_db:
                blast_db = get_blast_db(
                    fasta_path=genome_path,
                    db_dir=os.path.join(args.cache_dir, "blastdb"),
                    makeblastdb_path=args.makeblastdb,
                    force_rebuild=args.force_rebuild_blast,
                )

            hits = run_blast_primers(
                primer_sequences=[fwd_seq, rev_seq],
                blast_db=blast_db,
                blastn_path=args.blastn,
            )

            # Load GFF3 index for gene name annotation of off-target amplicons
            genomic_gene_map = None
            gff3_path = args.gff3
            if os.path.exists(gff3_path):
                try:
                    gene_index, _ = load_or_build_index(gff3_path, args.cache_dir)
                    if gene_index:
                        genomic_gene_map = build_genomic_gene_map(gene_index)
                except Exception:
                    pass

            # Use a dummy target to detect ALL off-target hits (no intended target)
            sensitivity = getattr(args, "off_target_sensitivity", "standard")
            preset = OFF_TARGET_PRESETS.get(sensitivity, OFF_TARGET_PRESETS["standard"])
            spec_result = evaluate_specificity(
                left_hits=hits.get(fwd_seq, []),
                right_hits=hits.get(rev_seq, []),
                target_seqid="",  # no target — all hits are "off-target"
                target_left_coords=(-1, -1),
                target_right_coords=(-1, -1),
                max_mismatches=preset["max_mismatches"],
                min_3prime_match=preset["min_3prime"],
                min_alignment_cover=preset["min_coverage"],
                max_amplicon_size=args.max_amplicon_size,
                genomic_gene_map=genomic_gene_map,
            )
            specificity_results = [spec_result]
            result['specificity'] = {
                'is_specific': spec_result.is_specific,
                'off_target_count': spec_result.off_target_count,
                'off_target_amplicons': [
                    {
                        'chromosome': a.chromosome,
                        'product_size': a.product_size,
                        'gene_name': a.gene_name,
                    }
                    for a in spec_result.off_target_amplicons
                ],
            }

    # Output
    if args.format == "json":
        output = _format_check_json(result)
    elif args.format == "tsv":
        output = _format_check_tsv(result)
    else:
        output = _format_check_text(result)

    if args.output == "-":
        print(output)
    else:
        with open(args.output, "w") as f:
            f.write(output)
        print(f"Results written to {args.output}")


def _format_check_text(result):
    """Format primer check result as human-readable text."""
    fwd = result['forward']
    rev = result['reverse']
    pair = result['pair']
    lines = []
    lines.append("=" * 60)
    lines.append("Primer Quality Check")
    lines.append("=" * 60)
    lines.append(f"Forward: {fwd['sequence']} ({fwd['length']} nt)")
    lines.append(f"  Tm={fwd['tm']}°C  GC={fwd['gc_percent']}%")
    lines.append(f"Reverse: {rev['sequence']} ({rev['length']} nt)")
    lines.append(f"  Tm={rev['tm']}°C  GC={rev['gc_percent']}%")
    lines.append(f"Tm difference: {pair['delta_tm']}°C")
    lines.append(f"\nTotal Score: {pair['total_score']}/100")
    lines.append(f"  Per-primer avg: {pair['per_primer_avg']}")
    lines.append(f"  Cross-scores:   {pair['cross_total']}")
    lines.append(f"\nForward Scores:")
    for k, v in fwd['scores'].items():
        if k != 'total':
            lines.append(f"  {k}: {v}")
    lines.append(f"Reverse Scores:")
    for k, v in rev['scores'].items():
        if k != 'total':
            lines.append(f"  {k}: {v}")
    lines.append(f"Cross Scores:")
    for k, v in pair['cross_scores'].items():
        lines.append(f"  {k}: {v}")

    spec = result.get('specificity')
    if spec:
        lines.append(f"\nSpecificity: {'PASS' if spec['is_specific'] else 'FAIL'}")
        lines.append(f"Off-target amplicons: {spec['off_target_count']}")
        for a in spec.get('off_target_amplicons', [])[:10]:
            lines.append(f"  {a['chromosome']}: {a['product_size']} bp"
                         f"{' gene=' + a['gene_name'] if a.get('gene_name') else ''}")

    return "\n".join(lines)


def _format_check_json(result):
    import json
    return json.dumps(result, indent=2)


def _format_check_tsv(result):
    fwd = result['forward']
    rev = result['reverse']
    pair = result['pair']
    spec = result.get('specificity', {})
    lines = ["forward_seq\treverse_seq\tfwd_tm\trev_tm\tfwd_gc\trev_gc\t"
             "delta_tm\ttotal_score\tis_specific\toff_target_count"]
    lines.append(f"{fwd['sequence']}\t{rev['sequence']}\t"
                 f"{fwd['tm']}\t{rev['tm']}\t{fwd['gc_percent']}\t{rev['gc_percent']}\t"
                 f"{pair['delta_tm']}\t{pair['total_score']}\t"
                 f"{spec.get('is_specific', '?')}\t{spec.get('off_target_count', '?')}")
    return "\n".join(lines)


def main():
    parser = build_parser()
    args = parser.parse_args()

    if args.serve:
        run_server(args)
    else:
        run(args)


if __name__ == "__main__":
    main()
