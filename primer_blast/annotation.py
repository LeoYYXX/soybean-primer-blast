"""Soybean-Arabidopsis annotation index.

Loads Soybean_Arabidopsis_Complete_Annotation.tsv and provides
per-gene lookup of gene name, Arabidopsis homolog, function description,
GO annotations, and BLAST/Homology notes.
"""

import os
from typing import Dict, Optional


ANNOTATION_COLS = [
    "soybean_gene_id",    # 大豆基因号
    "gene_name",          # 基因名
    "arabidopsis_homolog",# 拟南芥同源基因
    "function_desc",      # 功能描述
    "go_annotation",      # GO注释
    "notes",              # 备注
]


def _normalize_gene_id(gene_id: str) -> str:
    """Normalize gene ID for lookup: Glyma.01G000100.1.p -> Glyma.01G000100"""
    # Strip version suffixes like .1.p, .2.p, .Wm82.a4.v1
    # GFF3 short names are like Glyma.01G000100
    # Annotation IDs are like Glyma.01G000100.1.p
    parts = gene_id.split(".")
    if len(parts) >= 3:
        return ".".join(parts[:2])  # Glyma.01G000100
    return gene_id


def load_annotations(tsv_path: str) -> Dict[str, dict]:
    """Load the soybean-Arabidopsis annotation TSV into a lookup dict.

    Keys are normalized soybean gene short names (e.g. Glyma.01G000100).
    Values are dicts with all annotation fields.

    Also builds secondary index by full gene ID.
    """
    if not os.path.exists(tsv_path):
        return {}

    index: Dict[str, dict] = {}
    full_id_index: Dict[str, str] = {}  # full ID -> normalized short name

    with open(tsv_path, "r", encoding="utf-8-sig") as f:
        header = f.readline()  # skip header

        for line in f:
            line = line.rstrip("\r\n")
            if not line.strip():
                continue

            parts = line.split("\t")
            if len(parts) < 4:
                continue

            raw_gene_id = parts[0].strip() if len(parts) > 0 else ""
            gene_name = parts[1].strip() if len(parts) > 1 else ""
            ath_homolog = parts[2].strip() if len(parts) > 2 else ""
            function_desc = parts[3].strip() if len(parts) > 3 else ""
            go_annotation = parts[4].strip() if len(parts) > 4 else ""
            notes = parts[5].strip() if len(parts) > 5 else ""

            normalized = _normalize_gene_id(raw_gene_id)
            full_id_index[raw_gene_id] = normalized

            # Prefer the first entry if there are duplicates; merge GO/desc if already exists
            if normalized in index:
                existing = index[normalized]
                if go_annotation and not existing["go_annotation"]:
                    existing["go_annotation"] = go_annotation
                if function_desc and not existing["function_desc"]:
                    existing["function_desc"] = function_desc
            else:
                index[normalized] = {
                    "soybean_gene_id": raw_gene_id,
                    "gene_name": gene_name,
                    "arabidopsis_homolog": ath_homolog,
                    "function_desc": function_desc,
                    "go_annotation": go_annotation,
                    "notes": notes,
                }

    return index


def parse_go_terms(go_text: str) -> dict:
    """Parse GO annotation text into structured dict.

    Input: "P: photosystem II assembly (GO:0010207) | F: molecular_function (GO:0003674) | C: chloroplast (GO:0009507)"
    Returns: {"biological_process": [...], "molecular_function": [...], "cellular_component": [...]}
    """
    result = {
        "biological_process": [],
        "molecular_function": [],
        "cellular_component": [],
    }

    if not go_text:
        return result

    # Split by category
    for segment in go_text.split("|"):
        segment = segment.strip()
        if not segment:
            continue
        if segment.startswith("P:"):
            result["biological_process"].append(segment[2:].strip())
        elif segment.startswith("F:"):
            result["molecular_function"].append(segment[2:].strip())
        elif segment.startswith("C:"):
            result["cellular_component"].append(segment[2:].strip())

    return result


def get_gene_annotation(
    gene_short_name: str,
    annotation_index: Dict[str, dict],
) -> Optional[dict]:
    """Look up annotation for a gene by its short name or full ID.

    Returns None if not found.
    """
    if not annotation_index:
        return None

    # Direct lookup
    if gene_short_name in annotation_index:
        return annotation_index[gene_short_name]

    # Try normalized
    normalized = _normalize_gene_id(gene_short_name)
    if normalized in annotation_index:
        return annotation_index[normalized]

    return None
