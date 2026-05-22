"""Flask web server for Soybean Primer-BLAST.

Replicates the NCBI Primer-BLAST interactive experience as a local web app.
"""

import os
import sys
import time
import traceback
from typing import Dict, List, Optional

from flask import Flask, jsonify, render_template, request

from .blast_db import get_blast_db
from .cli import build_parser, parse_region
from .constants import (
    DEFAULT_BLASTN_PATH,
    DEFAULT_GENOME_FA,
    DEFAULT_GFF3,
    DEFAULT_MAKEBLASTDB_PATH,
    DEFAULT_MAX_AMPLICON_SIZE,
    DEFAULT_MAX_OFF_TARGET_MISMATCH,
    DEFAULT_MIN_3PRIME_MATCH,
)
from .fasta_index import load_fasta_index, extract_sequence
from .gff_index import load_or_build_index
from .primer_design import design_primers
from .primer_score import score_primer_pair
from .sequence import (
    extract_gene_sequence,
    extract_region_with_padding,
    genomic_to_template_coords,
    get_target_coordinates,
)
from .specificity import evaluate_specificity, run_blast_primers

# Module-level globals initialized at startup
_fai = None
_gene_index = None
_transcript_index = None
_blast_db_path: str = ""
_server_config: dict = {}

app = Flask(__name__)


def init_server(
    genome_path: str = DEFAULT_GENOME_FA,
    gff3_path: str = DEFAULT_GFF3,
    blast_db: str = "",
    blastn_path: str = DEFAULT_BLASTN_PATH,
    makeblastdb_path: str = DEFAULT_MAKEBLASTDB_PATH,
    cache_dir: str = ".primer_blast_cache",
) -> dict:
    """Initialize indices and BLAST database. Call once at startup.

    Returns a dict of status messages.
    """
    global _fai, _gene_index, _transcript_index, _blast_db_path, _server_config

    status = {}

    if not os.path.exists(genome_path):
        msg = f"Genome FASTA not found: {genome_path}"
        status["genome"] = msg
        print(msg, file=sys.stderr)
        return status

    if not os.path.exists(gff3_path):
        msg = f"GFF3 file not found: {gff3_path}"
        status["gff3"] = msg
        print(msg, file=sys.stderr)
        return status

    # Load FASTA index
    try:
        _fai = load_fasta_index(genome_path)
        status["fasta_index"] = f"Loaded {len(_fai)} sequences"
        print(status["fasta_index"])
    except Exception as e:
        status["fasta_index"] = f"Failed: {e}"
        print(f"WARNING: {e}", file=sys.stderr)

    # Load GFF3 index (uses pickle cache after first parse)
    try:
        _gene_index, _transcript_index = load_or_build_index(gff3_path, cache_dir)
        status["gff_index"] = f"Loaded {len(_transcript_index)} transcripts"
        print(status["gff_index"])
    except Exception as e:
        status["gff_index"] = f"Failed: {e}"
        print(f"WARNING: {e}", file=sys.stderr)

    # Resolve BLAST DB
    if blast_db:
        _blast_db_path = blast_db
    else:
        try:
            _blast_db_path = get_blast_db(
                fasta_path=genome_path,
                db_dir=os.path.join(cache_dir, "blastdb"),
                makeblastdb_path=makeblastdb_path,
                force_rebuild=False,
            )
            status["blast_db"] = f"Using: {_blast_db_path}"
        except Exception as e:
            status["blast_db"] = f"Failed: {e}"
            print(f"WARNING: BLAST DB not available: {e}", file=sys.stderr)

    _server_config = {
        "genome_path": genome_path,
        "gff3_path": gff3_path,
        "blast_db": _blast_db_path,
        "blastn_path": blastn_path,
        "makeblastdb_path": makeblastdb_path,
        "cache_dir": cache_dir,
    }
    return status


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/api/genes")
def api_genes():
    """Autocomplete endpoint: search gene names by partial match."""
    q = request.args.get("q", "").strip()
    if not q or len(q) < 2:
        return jsonify({"genes": []})

    if _transcript_index is None:
        return jsonify({"genes": [], "error": "GFF3 index not loaded"})

    # Case-insensitive prefix/suffix match, limit to 20
    q_lower = q.lower()
    matches = []
    for gene_id in _transcript_index:
        if q_lower in gene_id.lower():
            matches.append(gene_id)
        if len(matches) >= 20:
            break

    matches.sort()
    return jsonify({"genes": matches})


@app.route("/api/gene-structure/<gene_id>")
def api_gene_structure(gene_id):
    """Return expanded gene structure: promoter blocks, exons, introns, UTR, downstream blocks.

    Query params: upstream_kb (default 1), downstream_kb (default 4). Each kb = one clickable block.
    All coordinates are view-relative (1-based from 5' end).
    """
    if _transcript_index is None:
        return jsonify({"error": "GFF3 index not loaded"}), 500

    transcript = _transcript_index.get(gene_id)
    if not transcript:
        return jsonify({"error": f"Gene '{gene_id}' not found"}), 404

    strand = transcript.strand
    gs = transcript.start  # gene genomic start (smaller coord)
    ge = transcript.end    # gene genomic end (larger coord)
    seqid = transcript.seqid

    # Upstream/downstream flanking from query params (in kb, converted to bp)
    try:
        upstream_bp = int(request.args.get("upstream_kb", 1)) * 1000
    except (ValueError, TypeError):
        upstream_bp = 1000
    try:
        downstream_bp = int(request.args.get("downstream_kb", 4)) * 1000
    except (ValueError, TypeError):
        downstream_bp = 4000

    # Determine chromosome bounds from FASTA index
    chr_len = None
    if _fai and seqid in _fai:
        chr_len = _fai[seqid][1]

    # View boundaries: gene ± flanking, clamped to chromosome
    view_start_gen = max(1, gs - upstream_bp)
    view_end_gen = ge + downstream_bp
    if chr_len:
        view_end_gen = min(view_end_gen, chr_len)
    view_length = view_end_gen - view_start_gen + 1

    # Coordinate converter: genomic -> view (1-based)
    def g2v(gstart, gend):
        """Convert a genomic interval to view-relative coords.
        gstart/gend are in genomic coordinate space (gstart <= gend).
        For minus strand, higher genomic coords map to lower view coords.
        """
        if strand == "+":
            v1 = gstart - view_start_gen + 1
            v2 = gend - view_start_gen + 1
        else:
            v1 = view_end_gen - gend + 1
            v2 = view_end_gen - gstart + 1
        return (v1, v2)

    # Gene body in view coords
    gene_vs, gene_ve = g2v(gs, ge)

    # --- Promoter blocks (N × 1kb upstream of gene 5' end) ---
    num_promoter = max(1, upstream_bp // 1000)
    promoter_blocks = []
    for i in range(num_promoter):
        if strand == "+":
            pb_start = gs - upstream_bp + i * 1000
            pb_end = gs - upstream_bp + (i + 1) * 1000 - 1
        else:
            pb_start = ge + 1 + i * 1000
            pb_end = ge + (i + 1) * 1000

        pb_start = max(1, pb_start)
        pb_end = min(pb_end, view_end_gen) if i == num_promoter - 1 else pb_end
        if chr_len:
            pb_start = min(pb_start, chr_len)
            pb_end = min(pb_end, chr_len)
        if pb_start > view_end_gen or pb_end < view_start_gen:
            continue

        vs, ve = g2v(pb_start, pb_end)
        promoter_blocks.append({
            "index": len(promoter_blocks) + 1,
            "genomic_start": pb_start,
            "genomic_end": pb_end,
            "view_start": vs,
            "view_end": ve,
            "label": f"P{len(promoter_blocks) + 1}",
        })

    # --- Downstream blocks (N × 1kb downstream of gene 3' end) ---
    num_downstream = max(1, downstream_bp // 1000)
    downstream_blocks = []
    for i in range(num_downstream):
        if strand == "+":
            db_start = ge + 1 + i * 1000
            db_end = ge + (i + 1) * 1000
        else:
            db_start = gs - downstream_bp + i * 1000
            db_end = gs - downstream_bp + (i + 1) * 1000 - 1

        db_start = max(1, db_start)
        db_end = min(db_end, view_end_gen) if i == num_downstream - 1 else db_end
        if chr_len:
            db_start = min(db_start, chr_len)
            db_end = min(db_end, chr_len)
        if db_start > view_end_gen or db_end < view_start_gen:
            continue

        vs, ve = g2v(db_start, db_end)
        downstream_blocks.append({
            "index": len(downstream_blocks) + 1,
            "genomic_start": db_start,
            "genomic_end": db_end,
            "view_start": vs,
            "view_end": ve,
            "label": f"D{len(downstream_blocks) + 1}",
        })

    # --- Exons (use all_exons if available, otherwise reconstruct from CDS+UTR) ---
    exons = []
    if transcript.all_exons:
        for i, (ex_start, ex_end) in enumerate(transcript.all_exons):
            vs, ve = g2v(ex_start, ex_end)
            exons.append({
                "index": i + 1,
                "genomic_start": ex_start,
                "genomic_end": ex_end,
                "view_start": vs,
                "view_end": ve,
            })
    else:
        # Reconstruct exon boundaries from CDS + UTR regions (merge overlapping)
        all_parts = list(transcript.cds_exons) + list(transcript.five_prime_utr) + list(transcript.three_prime_utr)
        if all_parts:
            all_parts.sort(key=lambda c: c[0])
            merged = [list(all_parts[0])]
            for part in all_parts[1:]:
                if part[0] <= merged[-1][1] + 1:
                    merged[-1][1] = max(merged[-1][1], part[1])
                else:
                    merged.append(list(part))
            for i, (ex_start, ex_end) in enumerate(merged):
                vs, ve = g2v(ex_start, ex_end)
                exons.append({
                    "index": i + 1,
                    "genomic_start": ex_start,
                    "genomic_end": ex_end,
                    "view_start": vs,
                    "view_end": ve,
                })

    # --- Introns (gaps between consecutive exons) ---
    introns = []
    sorted_exons = sorted(exons, key=lambda e: e["view_start"])
    for i in range(1, len(sorted_exons)):
        gap_start = sorted_exons[i - 1]["view_end"] + 1
        gap_end = sorted_exons[i]["view_start"] - 1
        if gap_start <= gap_end:
            introns.append({
                "index": i,
                "view_start": gap_start,
                "view_end": gap_end,
            })

    # --- UTR regions (in view coords) ---
    utr5 = []
    for u_start, u_end in transcript.five_prime_utr:
        vs, ve = g2v(u_start, u_end)
        utr5.append({"view_start": vs, "view_end": ve})

    utr3 = []
    for u_start, u_end in transcript.three_prime_utr:
        vs, ve = g2v(u_start, u_end)
        utr3.append({"view_start": vs, "view_end": ve})

    return jsonify({
        "gene": gene_id,
        "seqid": seqid,
        "strand": strand,
        "view_start": view_start_gen,
        "view_end": view_end_gen,
        "view_length": view_length,
        "gene_view_start": gene_vs,
        "gene_view_end": gene_ve,
        "promoter_blocks": promoter_blocks,
        "downstream_blocks": downstream_blocks,
        "exons": exons,
        "introns": introns,
        "utr5_regions": utr5,
        "utr3_regions": utr3,
    })


def _compute_template_pattern(primer_seq, hit):
    """Build NCBI-style match_pattern showing template bases at mismatch positions.

    Extracts the full genomic context around the off-target binding site and
    compares base-by-base with the primer. Returns a string where:
      - '.' means primer base == template base (match)
      - A/C/G/T means primer base != template base (shows the *template* base)
    """
    if _fai is None or not _server_config.get("genome_path"):
        # Fallback: use BLAST alignment data only
        pattern = []
        qseq = (hit.qseq or "").upper()
        sseq = (hit.sseq or "").upper()
        qs = hit.query_start or 1
        qe = hit.query_end or len(primer_seq)
        for i, p in enumerate(primer_seq.upper()):
            idx = i + 1
            if idx < qs or idx > qe:
                pattern.append(p)
            else:
                j = idx - qs
                if j < len(qseq) and j < len(sseq):
                    pattern.append('.') if qseq[j] == sseq[j] else pattern.append(sseq[j])
                else:
                    pattern.append(p)
        return ''.join(pattern)

    is_rev = hit.subject_end < hit.subject_start
    qs = hit.query_start or 1
    plen = len(primer_seq)

    if is_rev:
        genomic_end = hit.subject_start + qs - 1
        genomic_start = genomic_end - plen + 1
        strand = "-"
    else:
        genomic_start = hit.subject_start - (qs - 1)
        genomic_end = genomic_start + plen - 1
        strand = "+"

    try:
        template_seq = extract_sequence(
            fasta_path=_server_config["genome_path"],
            seqid=hit.subject_id,
            start=max(1, genomic_start),
            end=genomic_end,
            strand=strand,
            index=_fai,
        )
    except (KeyError, Exception):
        template_seq = ""

    pattern = []
    for i, p in enumerate(primer_seq.upper()):
        if i < len(template_seq):
            t = template_seq[i].upper()
            pattern.append('.') if p == t else pattern.append(t)
        else:
            pattern.append(p)
    return ''.join(pattern)


@app.route("/api/design", methods=["POST"])
def api_design():
    """Core endpoint: design primers and check specificity."""
    t_start = time.time()

    try:
        data = request.get_json(force=True)
    except Exception:
        return jsonify({"success": False, "error": "Invalid JSON body"}), 400

    # --- Validate required fields ---
    input_mode = data.get("input_mode", "gene")
    if input_mode not in ("gene", "region", "sequence"):
        return jsonify({"success": False, "error": f"Unknown input_mode: {input_mode}"}), 400

    # --- Extract parameters ---
    params = {
        "target": data.get("target", "gene"),
        "cds_from": int(data.get("cds_from", 0)),
        "cds_to": int(data.get("cds_to", 0)),
        "flanking": int(data.get("flanking", 500)),
        "product_size_min": int(data.get("product_size_min", data.get("product_min", 100))),
        "product_size_max": int(data.get("product_size_max", data.get("product_max", 400))),
        "forward_region": data.get("forward_region"),  # [view_start, view_end] for visual mode
        "reverse_region": data.get("reverse_region"),  # [view_start, view_end] for visual mode
        "primer_opt_size": int(data.get("primer_opt_size", 20)),
        "primer_min_size": int(data.get("primer_min_size", 18)),
        "primer_max_size": int(data.get("primer_max_size", 25)),
        "primer_opt_tm": float(data.get("primer_opt_tm", 60.0)),
        "primer_min_tm": float(data.get("primer_min_tm", 57.0)),
        "primer_max_tm": float(data.get("primer_max_tm", 63.0)),
        "primer_gc_min": float(data.get("primer_gc_min", 35.0)),
        "primer_gc_max": float(data.get("primer_gc_max", 65.0)),
        "num_return": int(data.get("num_return", 5)),
        "skip_specificity": bool(data.get("skip_specificity", False)),
        "max_off_target_mismatch": int(data.get("max_off_target_mismatch", DEFAULT_MAX_OFF_TARGET_MISMATCH)),
        "min_3prime_match": int(data.get("min_3prime_match", DEFAULT_MIN_3PRIME_MATCH)),
        "max_amplicon_size": int(data.get("max_amplicon_size", DEFAULT_MAX_AMPLICON_SIZE)),
    }

    if params["product_size_min"] > params["product_size_max"]:
        return jsonify({"success": False, "error": "product_min must be <= product_max"}), 400
    if params["primer_min_size"] > params["primer_max_size"]:
        return jsonify({"success": False, "error": "primer_min_size must be <= primer_max_size"}), 400
    if params["primer_min_tm"] > params["primer_max_tm"]:
        return jsonify({"success": False, "error": "primer_min_tm must be <= primer_max_tm"}), 400

    warnings = []
    template = ""
    template_info = ""
    seqid = ""
    tgt_start = None
    tgt_len = None
    transcript = None

    # --- Step 1: Resolve input ---
    if input_mode == "gene":
        gene_name = (data.get("gene") or "").strip()
        if not gene_name:
            return jsonify({"success": False, "error": "Gene name is required for gene input mode"}), 400

        if _transcript_index is None:
            return jsonify({"success": False, "error": "GFF3 index not loaded. Check server startup."}), 500

        transcript = _transcript_index.get(gene_name)
        if not transcript:
            similar = [g for g in _transcript_index if gene_name[:8].lower() in g.lower()][:5]
            return jsonify({
                "success": False,
                "error": f"Gene '{gene_name}' not found in GFF3.",
                "suggestions": similar,
            }), 404

        seqid = transcript.seqid

        # --- Visual target mode: extract gene ± flanking, use view coords as template coords ---
        if params["target"] == "visual":
            up_bp = int(data.get("upstream_bp", 4000))
            dn_bp = int(data.get("downstream_bp", 4000))
            chr_len = None
            if _fai and seqid in _fai:
                chr_len = _fai[seqid][1]
            tpl_gen_start = max(1, transcript.start - up_bp)
            tpl_gen_end = transcript.end + dn_bp
            if chr_len:
                tpl_gen_end = min(tpl_gen_end, chr_len)

            template = extract_sequence(
                fasta_path=_server_config["genome_path"],
                seqid=seqid,
                start=tpl_gen_start,
                end=tpl_gen_end,
                strand=transcript.strand,
                index=_fai,
            )
            if not template:
                return jsonify({"success": False, "error": "Failed to extract expanded template"}), 500

            template_info = f"{seqid}:{tpl_gen_start}-{tpl_gen_end}({transcript.strand}) visual"
            view_length = tpl_gen_end - tpl_gen_start + 1

            # forward_region and reverse_region are view-relative [start, end] pairs
            fwd_reg = params["forward_region"]
            rev_reg = params["reverse_region"]

            if not fwd_reg or not rev_reg or len(fwd_reg) != 2 or len(rev_reg) != 2:
                return jsonify({
                    "success": False,
                    "error": "Visual mode requires forward_region and reverse_region as [start, end] pairs.",
                }), 400

            fwd_start, fwd_end = int(fwd_reg[0]), int(fwd_reg[1])
            rev_start, rev_end = int(rev_reg[0]), int(rev_reg[1])

            # Ensure forward is left of reverse in template coords
            if fwd_start > rev_start:
                fwd_start, fwd_end, rev_start, rev_end = rev_start, rev_end, fwd_start, fwd_end

            # Target = region between the two selected binding regions (0-based template coords)
            # fwd_start, fwd_end, rev_start, rev_end are 1-based view coords
            # Template position = view_coord - 1
            tgt_start = fwd_end  # 0-based, right after forward region
            tgt_len = max(1, rev_start - fwd_end - 1)  # gap between regions

            template_info += f" Fwd:{fwd_start}-{fwd_end} Rev:{rev_start}-{rev_end}"

        else:
            template = extract_gene_sequence(_server_config["genome_path"], transcript, _fai)
            if not template:
                return jsonify({"success": False, "error": "Failed to extract gene sequence"}), 500

        if params["target"] != "visual":
            template_info = f"{transcript.seqid}:{transcript.start}-{transcript.end}({transcript.strand})"
            template_info += f" {params['target']}"

        # Determine target coords
        if params["target"] == "visual":
            pass  # tgt_start, tgt_len already set above

        elif params["target"] == "cds_range":
            cds_from = params["cds_from"]
            cds_to = params["cds_to"]
            if (cds_from < 1 or cds_to < 1 or
                cds_from > len(transcript.cds_exons) or
                cds_to > len(transcript.cds_exons)):
                return jsonify({
                    "success": False,
                    "error": f"Invalid CDS range: {cds_from}-{cds_to}. Gene has {len(transcript.cds_exons)} CDS regions.",
                }), 400

            gs_start = transcript.cds_exons[cds_from - 1][0]
            gs_end = transcript.cds_exons[cds_to - 1][1]
            tgt_start, tgt_len = genomic_to_template_coords(
                gs_start, gs_end,
                transcript.start, transcript.end, transcript.strand,
            )
            template_info += f" CDS{cds_from}-CDS{cds_to}"

            # Auto-expand product size to accommodate target + flanking + primers
            flanking = params["flanking"]
            # Absolute minimum: target + room for one primer on each side
            abs_min = tgt_len + 2 * params["primer_min_size"]
            # Ideal: target + full flanking on both sides + optimal primers
            ideal = tgt_len + 2 * flanking + 2 * params["primer_opt_size"]

            if abs_min >= len(template):
                abs_min = max(tgt_len + 10, len(template) - 20)
                warnings.append(f"Target region ({tgt_len}bp) nearly fills the {len(template)}bp template.")

            params["product_size_min"] = max(params["product_size_min"], abs_min)
            params["product_size_max"] = max(params["product_size_max"], min(ideal, len(template)))
            # Ensure max never exceeds template length (Primer3 hard constraint)
            if params["product_size_max"] > len(template):
                params["product_size_max"] = len(template)
            if params["product_size_max"] < params["product_size_min"]:
                params["product_size_max"] = min(params["product_size_min"] + 200, len(template))

        elif params["target"] not in ("gene", "all"):
            target_coords = get_target_coordinates(transcript, params["target"])
            if target_coords:
                tgt_start, tgt_len = genomic_to_template_coords(
                    target_coords[0], target_coords[1],
                    transcript.start, transcript.end, transcript.strand,
                )
            else:
                warnings.append(f"No {params['target']} region found. Targeting full gene body.")

    elif input_mode == "region":
        region_raw = (data.get("region") or "").strip()
        if not region_raw:
            return jsonify({"success": False, "error": "Region string is required for region input mode"}), 400

        try:
            seqid, rstart, rend, rstrand = parse_region(region_raw)
        except ValueError as e:
            return jsonify({"success": False, "error": f"Invalid region: {e}"}), 400

        template = extract_region_with_padding(
            _server_config["genome_path"], seqid, rstart, rend, rstrand, index=_fai,
        )
        template_info = f"{seqid}:{rstart}-{rend}({rstrand}) region"

        # Numeric target (e.g. "200,600") is template-relative; named targets need transcript
        if params["target"] not in ("gene", "all"):
            target_coords = get_target_coordinates(None, params["target"])
            if target_coords:
                tgt_start, tgt_len = target_coords[0], target_coords[1] - target_coords[0] + 1

    elif input_mode == "sequence":
        raw_seq = (data.get("sequence") or "").strip()
        if not raw_seq:
            return jsonify({"success": False, "error": "Sequence is required for sequence input mode"}), 400

        # Could be file path or raw sequence
        if os.path.exists(raw_seq):
            with open(raw_seq, "r") as f:
                content = f.read().strip()
            if content.startswith(">"):
                lines = content.split("\n")
                template = "".join(l.strip() for l in lines if not l.startswith(">")).upper()
            else:
                template = content.upper()
            template_info = f"File: {raw_seq} ({len(template)} bp)"
        else:
            template = raw_seq.upper()
            # Validate DNA
            if not all(c in "ACGTN" for c in template):
                return jsonify({"success": False, "error": "Sequence contains non-DNA characters (only A,C,G,T,N allowed)"}), 400
            template_info = f"Direct input ({len(template)} bp)"

        if len(template) < 50:
            return jsonify({"success": False, "error": "Template sequence too short (minimum 50 bp)"}), 400

        # For sequence mode, numeric target is template-relative
        if params["target"] not in ("gene", "all"):
            target_coords = get_target_coordinates(None, params["target"])
            if target_coords:
                tgt_start, tgt_len = target_coords[0], target_coords[1] - target_coords[0] + 1

    if not template:
        return jsonify({"success": False, "error": "Could not extract template sequence"}), 500

    # --- Step 2: Design primers ---
    try:
        pairs = design_primers(
            template=template,
            target_start=tgt_start,
            target_length=tgt_len,
            product_size_min=params["product_size_min"],
            product_size_max=params["product_size_max"],
            primer_opt_size=params["primer_opt_size"],
            primer_min_size=params["primer_min_size"],
            primer_max_size=params["primer_max_size"],
            primer_opt_tm=params["primer_opt_tm"],
            primer_min_tm=params["primer_min_tm"],
            primer_max_tm=params["primer_max_tm"],
            primer_gc_min=params["primer_gc_min"],
            primer_gc_max=params["primer_gc_max"],
            num_return=params["num_return"],
        )
    except Exception as e:
        return jsonify({
            "success": False,
            "error": f"Primer3 design failed: {e}",
            "suggestion": "Try relaxing constraints (wider Tm/GC range, larger product size).",
        }), 500

    if not pairs:
        # Fallback: retry with relaxed constraints
        try:
            pairs = design_primers(
                template=template,
                target_start=tgt_start,
                target_length=tgt_len,
                product_size_min=min(70, params["product_size_min"]),
                product_size_max=max(params["product_size_max"], 1000),
                primer_opt_size=params["primer_opt_size"],
                primer_min_size=max(15, params["primer_min_size"] - 3),
                primer_max_size=params["primer_max_size"] + 3,
                primer_opt_tm=params["primer_opt_tm"],
                primer_min_tm=max(50, params["primer_min_tm"] - 5),
                primer_max_tm=min(70, params["primer_max_tm"] + 5),
                primer_gc_min=max(20, params["primer_gc_min"] - 10),
                primer_gc_max=min(80, params["primer_gc_max"] + 10),
                num_return=params["num_return"],
            )
        except Exception:
            pass

    if not pairs:
        return jsonify({
            "success": False,
            "error": "No primer pairs could be designed with the given parameters.",
            "suggestion": "Try relaxing constraints (wider Tm/GC range, larger product size).",
        }), 200

    # --- Post-filter for visual mode: enforce primer binding regions ---
    if params["target"] == "visual":
        filtered = []
        for pair in pairs:
            # fwd_start/end are 1-based view coords; pair positions are 0-based template coords
            fwd_in = (fwd_start - 1) <= pair.left_position <= (fwd_end - 1)
            rev_in = (rev_start - 1) <= pair.right_position <= (rev_end - 1)
            if fwd_in and rev_in:
                filtered.append(pair)
        pairs = filtered
        if not pairs:
            return jsonify({
                "success": False,
                "error": "No primer pairs landed within the selected binding regions.",
                "suggestion": "Try selecting wider regions or relaxing primer parameters.",
            }), 200

    # --- Step 3: Specificity checking ---
    specificity_results = None

    if not params["skip_specificity"] and _blast_db_path:
        try:
            all_primers = set()
            for p in pairs:
                all_primers.add(p.left_primer)
                all_primers.add(p.right_primer)

            hits = run_blast_primers(
                primer_sequences=list(all_primers),
                blast_db=_blast_db_path,
                blastn_path=_server_config["blastn_path"],
            )

            specificity_results = []
            for pair in pairs:
                left_hits = hits.get(pair.left_primer, [])
                right_hits = hits.get(pair.right_primer, [])

                # Determine genomic coordinates for target site identification
                if seqid and transcript:
                    gstart = transcript.start
                    gend = transcript.end
                    gstrand = transcript.strand

                    if gstrand == "+":
                        left_gen_start = gstart + pair.left_position
                        left_gen_end = left_gen_start + len(pair.left_primer) - 1
                        right_gen_end = gstart + pair.right_position
                        right_gen_start = right_gen_end - len(pair.right_primer) + 1
                    else:
                        left_gen_end = gend - pair.left_position
                        left_gen_start = left_gen_end - len(pair.left_primer) + 1
                        right_gen_start = gend - pair.right_position
                        right_gen_end = right_gen_start + len(pair.right_primer) - 1

                    left_coords = (min(left_gen_start, left_gen_end), max(left_gen_start, left_gen_end))
                    right_coords = (min(right_gen_start, right_gen_end), max(right_gen_start, right_gen_end))
                elif seqid:
                    # Region mode: approximate coords from template positions
                    left_coords = (pair.left_position, pair.left_position + len(pair.left_primer))
                    right_coords = (pair.right_position, pair.right_position + len(pair.right_primer))
                else:
                    left_coords = (0, len(pair.left_primer))
                    right_coords = (0, len(pair.right_primer))

                result = evaluate_specificity(
                    left_hits=left_hits,
                    right_hits=right_hits,
                    target_seqid=seqid if seqid else "",
                    target_left_coords=left_coords,
                    target_right_coords=right_coords,
                    max_mismatches=params["max_off_target_mismatch"],
                    min_3prime_match=params["min_3prime_match"],
                    max_amplicon_size=params["max_amplicon_size"],
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

        except Exception as e:
            warnings.append(f"Specificity check failed: {e}")

    # --- Step 4: Build response ---
    elapsed = round(time.time() - t_start, 1)

    response = {
        "success": True,
        "template": template_info,
        "template_length": len(template),
        "target_region": [tgt_start, tgt_len] if tgt_start is not None else None,
        "pairs": [],
        "warnings": warnings,
        "elapsed_seconds": elapsed,
        "parameters": {
            "input_mode": input_mode,
            **params,
        },
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

        # Primer quality score (100-point scale)
        try:
            pair_dict["score"] = score_primer_pair(pair)
        except Exception:
            pair_dict["score"] = None

        if specificity_results and i < len(specificity_results):
            sr = specificity_results[i]
            pair_dict["specificity"] = {
                "is_specific": sr.is_specific,
                "off_target_count": sr.off_target_count,
                "off_target_amplicons": [
                    {
                        "chromosome": a.chromosome,
                        "left_position": a.left_hit.subject_start,
                        "right_position": a.right_hit.subject_end,
                        "product_size": a.product_size,
                        "forward_primer": a.left_hit.primer_seq,
                        "reverse_primer": a.right_hit.primer_seq,
                        "forward_alignment": {
                            "qseq": a.left_hit.qseq,
                            "sseq": a.left_hit.sseq,
                            "query_start": a.left_hit.query_start,
                            "query_end": a.left_hit.query_end,
                            "subject_start": a.left_hit.subject_start,
                            "subject_end": a.left_hit.subject_end,
                            "num_mismatch": a.left_hit.num_mismatch,
                            "num_gap": a.left_hit.num_gap,
                        },
                        "reverse_alignment": {
                            "qseq": a.right_hit.qseq,
                            "sseq": a.right_hit.sseq,
                            "query_start": a.right_hit.query_start,
                            "query_end": a.right_hit.query_end,
                            "subject_start": a.right_hit.subject_start,
                            "subject_end": a.right_hit.subject_end,
                            "num_mismatch": a.right_hit.num_mismatch,
                            "num_gap": a.right_hit.num_gap,
                        },
                        "forward_template_pattern": _compute_template_pattern(
                            a.left_hit.primer_seq, a.left_hit,
                        ),
                        "reverse_template_pattern": _compute_template_pattern(
                            a.right_hit.primer_seq, a.right_hit,
                        ),
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

        response["pairs"].append(pair_dict)

    return jsonify(response)


def run_server(args):
    """Initialize server and start Flask app. Called when --serve is passed."""
    init_server(
        genome_path=args.genome,
        gff3_path=args.gff3,
        blast_db=args.blast_db or "",
        blastn_path=args.blastn,
        makeblastdb_path=args.makeblastdb,
        cache_dir=args.cache_dir,
    )

    print(f"\nStarting Soybean Primer-BLAST web server...")
    print(f"Open http://{args.host}:{args.port} in your browser\n")

    # Ensure the templates directory is findable
    import os as _os
    template_dir = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "templates")
    app.template_folder = template_dir

    app.run(
        host=args.host,
        port=args.port,
        debug=args.debug,
    )
