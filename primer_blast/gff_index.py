"""GFF3 parser for soybean Phytozome v13 gene annotations.

Builds in-memory indices of gene coordinates and longest transcript structures.
Results are cached as pickle files to avoid re-parsing the 125 MB GFF3 on every run.
"""

import os
import pickle
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


@dataclass
class GeneRecord:
    """A single gene from the GFF3."""
    gene_id: str          # e.g. "Glyma.01G000100.Wm82.a4.v1"
    short_name: str       # e.g. "Glyma.01G000100"
    seqid: str            # e.g. "Gm01"
    start: int            # 1-based inclusive
    end: int              # 1-based inclusive
    strand: str           # "+" or "-"


@dataclass
class TranscriptRecord:
    """Longest mRNA transcript for a gene, with CDS and UTR positions."""
    mrna_id: str
    gene_short_name: str
    seqid: str
    start: int
    end: int
    strand: str
    cds_exons: List[Tuple[int, int]] = field(default_factory=list)
    all_exons: List[Tuple[int, int]] = field(default_factory=list)
    five_prime_utr: List[Tuple[int, int]] = field(default_factory=list)
    three_prime_utr: List[Tuple[int, int]] = field(default_factory=list)

    @property
    def cds_start(self) -> Optional[int]:
        if not self.cds_exons:
            return None
        return min(e[0] for e in self.cds_exons)

    @property
    def cds_end(self) -> Optional[int]:
        if not self.cds_exons:
            return None
        return max(e[1] for e in self.cds_exons)


def _parse_gff3_attributes(attr_str: str) -> Dict[str, str]:
    """Parse GFF3 column-9 attribute string into a dict."""
    result = {}
    for part in attr_str.split(";"):
        part = part.strip()
        if "=" in part:
            key, val = part.split("=", 1)
            result[key.strip()] = val.strip()
    return result


def build_gene_index(gff3_path: str) -> Tuple[Dict[str, GeneRecord], Dict[str, TranscriptRecord]]:
    """Parse GFF3 and return (gene_index, transcript_index).

    gene_index: keyed by short gene name (e.g. "Glyma.01G000100")
    transcript_index: keyed by short gene name, value is the longest transcript
    """
    gene_index: Dict[str, GeneRecord] = {}
    transcript_index: Dict[str, TranscriptRecord] = {}
    gene_id_to_short: Dict[str, str] = {}  # gene full ID -> short name

    # First pass: collect all gene records
    with open(gff3_path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue

            seqid, source, feature, start, end, score, strand, phase, attrs = cols

            if feature == "gene":
                attr = _parse_gff3_attributes(attrs)
                gene_id = attr.get("ID", "")
                short_name = attr.get("Name", "")
                if not short_name and gene_id:
                    short_name = gene_id.rsplit(".", 2)[0]
                gene_index[short_name] = GeneRecord(
                    gene_id=gene_id,
                    short_name=short_name,
                    seqid=seqid,
                    start=int(start),
                    end=int(end),
                    strand=strand,
                )
                gene_id_to_short[gene_id] = short_name

    # Second pass: collect mRNA features, picking the longest transcript per gene
    gene_mrnas: Dict[str, List[str]] = {}  # short_gene_name -> list of mrna_ids
    mrna_info: Dict[str, dict] = {}        # mrna_id -> {start, end, strand, gene_short_name, is_longest}

    with open(gff3_path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue

            seqid, source, feature, start, end, score, strand, phase, attrs = cols

            if feature == "mRNA":
                attr = _parse_gff3_attributes(attrs)
                mrna_id = attr.get("ID", "")
                parent = attr.get("Parent", "")
                is_longest = attr.get("longest", "0") == "1"

                # Find parent gene short name via gene_id lookup
                parent_short = gene_id_to_short.get(parent, "")
                if not parent_short and parent:
                    # Fallback: try rsplit logic for edge cases
                    parent_short = parent.rsplit(".", 2)[0]
                if parent_short not in gene_mrnas:
                    gene_mrnas[parent_short] = []
                gene_mrnas[parent_short].append(mrna_id)

                mrna_info[mrna_id] = {
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "seqid": seqid,
                    "gene_short_name": parent_short,
                    "is_longest": is_longest,
                }

    # Third pass: collect CDS and UTR features, keyed by mRNA parent
    mrna_children: Dict[str, List[dict]] = {}

    with open(gff3_path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue

            seqid, source, feature, start, end, score, strand, phase, attrs = cols

            if feature in ("CDS", "five_prime_UTR", "three_prime_UTR", "exon"):
                attr = _parse_gff3_attributes(attrs)
                parent = attr.get("Parent", "")
                if parent not in mrna_children:
                    mrna_children[parent] = []
                mrna_children[parent].append({
                    "feature": feature,
                    "start": int(start),
                    "end": int(end),
                })

    # Build transcript records for the longest mRNA of each gene
    for short_name, gene in gene_index.items():
        mrna_ids = gene_mrnas.get(short_name, [])

        if not mrna_ids:
            continue

        # Pick longest transcript: prefer is_longest=1, then longest genomic span
        best_mrna_id = None
        best_span = 0
        for mid in mrna_ids:
            info = mrna_info.get(mid)
            if not info:
                continue
            span = info["end"] - info["start"]
            if info["is_longest"]:
                # explicitly flagged as longest — always prefer
                if best_mrna_id and mrna_info.get(best_mrna_id, {}).get("is_longest"):
                    if span > best_span:
                        best_mrna_id = mid
                        best_span = span
                else:
                    best_mrna_id = mid
                    best_span = span
            elif not best_mrna_id:
                best_mrna_id = mid
                best_span = span
            elif not mrna_info.get(best_mrna_id, {}).get("is_longest") and span > best_span:
                best_mrna_id = mid
                best_span = span

        if not best_mrna_id:
            continue

        info = mrna_info[best_mrna_id]
        tr = TranscriptRecord(
            mrna_id=best_mrna_id,
            gene_short_name=short_name,
            seqid=info["seqid"],
            start=info["start"],
            end=info["end"],
            strand=info["strand"],
        )

        # Collect CDS and UTR features
        children = mrna_children.get(best_mrna_id, [])
        for child in children:
            coord = (child["start"], child["end"])
            if child["feature"] == "CDS":
                tr.cds_exons.append(coord)
            elif child["feature"] == "five_prime_UTR":
                tr.five_prime_utr.append(coord)
            elif child["feature"] == "three_prime_UTR":
                tr.three_prime_utr.append(coord)
            elif child["feature"] == "exon":
                tr.all_exons.append(coord)

        # Sort exons by genomic position (ascending)
        tr.cds_exons.sort(key=lambda c: c[0])
        tr.all_exons.sort(key=lambda c: c[0])
        tr.five_prime_utr.sort(key=lambda c: c[0])
        tr.three_prime_utr.sort(key=lambda c: c[0])

        transcript_index[short_name] = tr

    return gene_index, transcript_index


def load_or_build_index(
    gff3_path: str,
    cache_dir: str = ".primer_blast_cache",
) -> Tuple[Dict[str, GeneRecord], Dict[str, TranscriptRecord]]:
    """Load cached gene/transcript indices, or rebuild if GFF3 is newer."""
    os.makedirs(cache_dir, exist_ok=True)

    gff3_mtime = os.path.getmtime(gff3_path)
    cache_file = os.path.join(cache_dir, "gene_index.pkl")

    if os.path.exists(cache_file):
        cache_mtime = os.path.getmtime(cache_file)
        if cache_mtime >= gff3_mtime:
            with open(cache_file, "rb") as f:
                return pickle.load(f)

    gene_index, transcript_index = build_gene_index(gff3_path)

    with open(cache_file, "wb") as f:
        pickle.dump((gene_index, transcript_index), f, protocol=pickle.HIGHEST_PROTOCOL)

    return gene_index, transcript_index
