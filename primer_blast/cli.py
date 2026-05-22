"""Command-line interface for Soybean Primer-BLAST."""

import argparse
import os
import sys
from typing import Optional

from .constants import (
    DEFAULT_BLASTN_PATH,
    DEFAULT_GENOME_FA,
    DEFAULT_GFF3,
    DEFAULT_MAKEBLASTDB_PATH,
    DEFAULT_MAX_AMPLICON_SIZE,
    DEFAULT_MAX_OFF_TARGET_MISMATCH,
    DEFAULT_MIN_3PRIME_MATCH,
    DEFAULT_NUM_RETURN,
    DEFAULT_PRIMER_GC_MAX,
    DEFAULT_PRIMER_GC_MIN,
    DEFAULT_PRIMER_MAX_SIZE,
    DEFAULT_PRIMER_MAX_TM,
    DEFAULT_PRIMER_MIN_SIZE,
    DEFAULT_PRIMER_MIN_TM,
    DEFAULT_PRIMER_OPT_SIZE,
    DEFAULT_PRIMER_OPT_TM,
    DEFAULT_PRIMER_PRODUCT_MAX,
    DEFAULT_PRIMER_PRODUCT_MIN,
)


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser."""
    parser = argparse.ArgumentParser(
        description="Soybean Primer-BLAST — local primer design with specificity checking",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python -m primer_blast --gene Glyma.01G000100
  python -m primer_blast --gene Glyma.01G000100 --target cds --num-return 10
  python -m primer_blast --region Gm01:27344-28430
  python -m primer_blast --region Gm01:27344-28430:- --product-min 150 --product-max 500
  python -m primer_blast --gene Glyma.01G000100 --skip-specificity
  python -m primer_blast --gene Glyma.01G000100 -f json -o primers.json
        """,
    )

    # Input
    input_group = parser.add_mutually_exclusive_group()
    input_group.add_argument(
        "--gene", "-g", type=str,
        help="Gene short name, e.g. Glyma.01G000100",
    )
    input_group.add_argument(
        "--region", "-r", type=str,
        help="Genomic region: chr:start-end[:strand], e.g. Gm01:27344-28430 or Gm01:27344-28430:-",
    )
    input_group.add_argument(
        "--sequence", "-s", type=str,
        help="Direct FASTA sequence or path to FASTA file containing the template",
    )

    # Target
    parser.add_argument(
        "--target", "-t", type=str, default="gene",
        help="Sub-region to amplify: 'cds', 'utr', 'gene', or 'start,end' coords. "
             "For --gene input, defaults to full gene body. Use 'cds' for CDS-specific primers.",
    )

    # Primer3 parameters
    p3_group = parser.add_argument_group("Primer parameters (Primer3)")
    p3_group.add_argument("--product-min", type=int, default=DEFAULT_PRIMER_PRODUCT_MIN)
    p3_group.add_argument("--product-max", type=int, default=DEFAULT_PRIMER_PRODUCT_MAX)
    p3_group.add_argument("--primer-opt-size", type=int, default=DEFAULT_PRIMER_OPT_SIZE)
    p3_group.add_argument("--primer-min-size", type=int, default=DEFAULT_PRIMER_MIN_SIZE)
    p3_group.add_argument("--primer-max-size", type=int, default=DEFAULT_PRIMER_MAX_SIZE)
    p3_group.add_argument("--primer-opt-tm", type=float, default=DEFAULT_PRIMER_OPT_TM)
    p3_group.add_argument("--primer-min-tm", type=float, default=DEFAULT_PRIMER_MIN_TM)
    p3_group.add_argument("--primer-max-tm", type=float, default=DEFAULT_PRIMER_MAX_TM)
    p3_group.add_argument("--primer-gc-min", type=float, default=DEFAULT_PRIMER_GC_MIN)
    p3_group.add_argument("--primer-gc-max", type=float, default=DEFAULT_PRIMER_GC_MAX)
    p3_group.add_argument("--num-return", "-n", type=int, default=DEFAULT_NUM_RETURN)

    # Specificity
    spec_group = parser.add_argument_group("Specificity checking")
    spec_group.add_argument(
        "--skip-specificity", action="store_true",
        help="Skip BLAST specificity check (faster, but no filtering)",
    )
    spec_group.add_argument(
        "--max-off-target-mismatch", type=int, default=DEFAULT_MAX_OFF_TARGET_MISMATCH,
        help="Max mismatches allowed for off-target consideration",
    )
    spec_group.add_argument(
        "--min-3prime-match", type=int, default=DEFAULT_MIN_3PRIME_MATCH,
        help="Min consecutive matches required at primer 3' end",
    )
    spec_group.add_argument(
        "--max-amplicon-size", type=int, default=DEFAULT_MAX_AMPLICON_SIZE,
        help="Max amplicon size for off-target detection",
    )

    # File paths
    path_group = parser.add_argument_group("File paths")
    path_group.add_argument("--genome", type=str, default=DEFAULT_GENOME_FA)
    path_group.add_argument("--gff3", type=str, default=DEFAULT_GFF3)
    path_group.add_argument("--blast-db", type=str, help="Pre-built BLAST DB (skip auto-build)")
    path_group.add_argument("--blastn", type=str, default=DEFAULT_BLASTN_PATH)
    path_group.add_argument("--makeblastdb", type=str, default=DEFAULT_MAKEBLASTDB_PATH)
    path_group.add_argument("--cache-dir", type=str, default=".primer_blast_cache")

    # Output
    parser.add_argument("--output", "-o", type=str, default="-", help="Output file (default: stdout)")
    parser.add_argument(
        "--format", "-f", type=str, choices=["text", "json", "tsv"], default="text",
        help="Output format",
    )

    # Misc
    parser.add_argument("--verbose", "-v", action="store_true")
    parser.add_argument(
        "--force-rebuild-blast", action="store_true",
        help="Force rebuild of BLAST database",
    )

    # Web server
    parser.add_argument(
        "--serve", action="store_true",
        help="Start as a local web server instead of CLI mode",
    )
    parser.add_argument("--port", type=int, default=5000, help="Server port (default: 5000)")
    parser.add_argument("--host", type=str, default="127.0.0.1", help="Server host (default: 127.0.0.1)")
    parser.add_argument("--debug", action="store_true", help="Enable Flask debug mode")

    return parser


def parse_region(region_str: str):
    """Parse a region string like 'Gm01:27344-28430' or 'Gm01:27344-28430:-'.

    Returns (seqid, start, end, strand).
    """
    strand = "+"
    if region_str.endswith(":-"):
        strand = "-"
        region_str = region_str[:-2]
    elif region_str.endswith(":+"):
        region_str = region_str[:-2]

    if ":" not in region_str:
        raise ValueError(f"Invalid region format: {region_str}. Expected chr:start-end")

    seqid, coords = region_str.split(":", 1)
    if "-" not in coords:
        raise ValueError(f"Invalid region format: {region_str}. Expected chr:start-end")

    start_str, end_str = coords.split("-", 1)
    return seqid, int(start_str), int(end_str), strand
