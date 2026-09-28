#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Functional enrichment from an annotation TSV.

Features
--------
1) GO enrichment using the annotation TSV as universe.
2) KEYWORD enrichment using the annotation TSV as universe.
3) Optional GO namespace split (BP/MF/CC) if go-basic.obo is provided.
4) Optional genus abundance summaries from a secondary hit/taxonomy TSV containing
   per-gene hit species, for example viral hit annotations.
5) Dotplots and barplots for GO / KEYWORD.
6) Optional treemap for KEYWORD enrichment.
7) Fixed dotplot bubble size legend at 1, 5, 10, 50 and 100 query genes.
8) Optional genus abundance plots and gene-level genus-richness plots.
9) Gene-level genus-richness labels use: GENENAME (ID).
10) All outputs are written into an output directory; default: enrichment.

Main universe for GO/KEYWORD
---------------------------
Expected main annotation file format by default:
#ID    GENENAME    DESCRIPTION    ENZYME    GO    KEYWORD    PATHWAY ...

Optional taxonomy/genus universe
--------------------------------
Expected optional taxa file format by default:
ID    GENENAME    DESCRIPTION    HIT_ID    HIT_IDENTITY    HIT_QCOV    HIT_TCOV
HIT_ALNLEN    HIT_EVALUE    HIT_SPECIES    HIT_TAXID

The taxa file is used only if --taxa-hits is supplied. Genus is extracted from
HIT_SPECIES using the first meaningful token, with special handling for
"Candidatus <Genus>". The script does not perform genus enrichment; instead it
reports genus abundance in the query set and which query genes show the highest
genus diversity across their hits.

Example:
    python3 sma3s_enrichment.py -a annotation.tsv -i genes_of_interest.txt \
        -d enrichment_results --go-obo go-basic.obo --min-term-size 2 \
        --alpha 0.05 --top 25 --use-namespace-fdr-for-namespace-plots \
        --taxa-hits perfect_hits.perfect_hits.long.tsv --make-keyword-treemap
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import pandas as pd
from scipy.stats import hypergeom

# Non-interactive backend: useful on clusters/servers without display.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import squarify  # type: ignore
    HAS_SQUARIFY = True
except Exception:
    squarify = None
    HAS_SQUARIFY = False


MISSING_VALUES = {"", ".", "-", "NA", "NaN", "nan", "None", "none", "null", "NULL"}
GO_NAMESPACES = ["biological_process", "molecular_function", "cellular_component"]
GO_NAMESPACE_LABELS = {
    "biological_process": "GO Biological Process",
    "molecular_function": "GO Molecular Function",
    "cellular_component": "GO Cellular Component",
    "unknown": "GO unknown namespace",
}
FIXED_BUBBLE_LEGEND_COUNTS = [1, 5, 10, 50, 100]


def eprint(*args, **kwargs) -> None:
    print(*args, file=sys.stderr, **kwargs)


def clean_header_name(x: str) -> str:
    return x.replace("\ufeff", "").strip()


def split_terms(value: Optional[str], delimiter_regex: str = r";") -> List[str]:
    if value is None:
        return []
    value = str(value).strip()
    if value in MISSING_VALUES:
        return []
    terms: List[str] = []
    for term in re.split(delimiter_regex, value):
        term = term.strip()
        if term and term not in MISSING_VALUES:
            terms.append(term)
    return terms


def read_query_ids(path: str) -> Set[str]:
    """
    Read query IDs from a text file.

    One ID per line is expected. If a line contains multiple columns, only the
    first column is used. A header-like first token is ignored.
    """
    ids: Set[str] = set()
    header_like = {"#ID", "ID", "Gene", "GeneID", "GENE", "GENENAME", "gene", "gene_id"}
    with open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            first = re.split(r"\t|,|\s+", line, maxsplit=1)[0].strip()
            if first in header_like:
                continue
            if first:
                ids.add(first)
    return ids


def parse_go_obo_metadata(path: str, include_obsolete: bool = False) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Minimal OBO parser: GO ID -> name and GO ID -> namespace."""
    go_names: Dict[str, str] = {}
    go_namespaces: Dict[str, str] = {}

    current_id: Optional[str] = None
    current_name: Optional[str] = None
    current_namespace: Optional[str] = None
    current_obsolete = False

    def flush_term() -> None:
        nonlocal current_id, current_name, current_namespace, current_obsolete
        if current_id and (include_obsolete or not current_obsolete):
            if current_name:
                go_names[current_id] = current_name
            if current_namespace:
                go_namespaces[current_id] = current_namespace
        current_id = None
        current_name = None
        current_namespace = None
        current_obsolete = False

    with open(path, "rt", encoding="utf-8", errors="replace") as handle:
        in_term = False
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "[Term]":
                if in_term:
                    flush_term()
                in_term = True
                continue
            if line.startswith("[") and line.endswith("]"):
                if in_term:
                    flush_term()
                in_term = False
                continue
            if not in_term:
                continue
            if line.startswith("id: "):
                current_id = line.split("id: ", 1)[1].strip()
            elif line.startswith("name: "):
                current_name = line.split("name: ", 1)[1].strip()
            elif line.startswith("namespace: "):
                current_namespace = line.split("namespace: ", 1)[1].strip()
            elif line.startswith("is_obsolete: true"):
                current_obsolete = True
        if in_term:
            flush_term()

    return go_names, go_namespaces


def read_go_names_tsv(path: str) -> Dict[str, str]:
    """
    Read GO names from a TSV/CSV-like file.

    A two-column file works: GO:0005525<TAB>GTP binding.
    A header with GO/name-like columns also works.
    """
    names: Dict[str, str] = {}
    sep = "\t" if path.lower().endswith((".tsv", ".txt")) else ","

    with open(path, "rt", encoding="utf-8", errors="replace") as handle:
        sample = handle.readline()
        if not sample:
            return names
        handle.seek(0)

        header = [clean_header_name(x) for x in sample.rstrip("\n").split(sep)]
        has_header = any(
            h.lower() in {"go", "go_id", "goid", "id", "term_id", "name", "term", "description"}
            for h in header
        )

        if has_header:
            reader = csv.DictReader(handle, delimiter=sep)
            reader.fieldnames = [clean_header_name(x) for x in (reader.fieldnames or [])]
            fields = reader.fieldnames or []
            id_candidates = ["GO", "GO_ID", "GOID", "ID", "TERM_ID", "go", "go_id", "goid", "id", "term_id"]
            name_candidates = ["NAME", "TERM", "DESCRIPTION", "GO_NAME", "name", "term", "description", "go_name"]
            id_col = next((c for c in id_candidates if c in fields), None)
            name_col = next((c for c in name_candidates if c in fields), None)
            if not id_col or not name_col:
                raise ValueError(f"Could not autodetect GO/name columns in {path}. Fields: {fields}")
            for row in reader:
                go_id = str(row.get(id_col, "")).strip()
                name = str(row.get(name_col, "")).strip()
                if go_id.startswith("GO:") and name:
                    names[go_id] = name
        else:
            reader2 = csv.reader(handle, delimiter=sep)
            for row in reader2:
                if len(row) < 2:
                    continue
                go_id = row[0].strip()
                name = row[1].strip()
                if go_id.startswith("GO:") and name:
                    names[go_id] = name

    return names


def read_annotation_data(
    annotations_tsv: str,
    id_col: str,
    gene_name_col: str,
    go_col: str,
    keyword_col: str,
    sep: str,
    term_delimiter_regex: str,
) -> Tuple[Set[str], Dict[str, Set[str]], Dict[str, Set[str]], Dict[str, str]]:
    """
    Read main annotation TSV.

    Returns:
      universe_ids
      go_term_to_ids
      keyword_term_to_ids
      gene_name_map, used for labels like GENENAME (ID)
    """
    universe_ids: Set[str] = set()
    go_term_to_ids: Dict[str, Set[str]] = defaultdict(set)
    keyword_term_to_ids: Dict[str, Set[str]] = defaultdict(set)
    gene_name_map: Dict[str, str] = {}

    with open(annotations_tsv, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=sep)
        if reader.fieldnames is None:
            raise ValueError("The annotation TSV appears to have no header.")

        reader.fieldnames = [clean_header_name(x) for x in reader.fieldnames]
        fields = set(reader.fieldnames)

        id_col = clean_header_name(id_col)
        gene_name_col = clean_header_name(gene_name_col)
        go_col = clean_header_name(go_col)
        keyword_col = clean_header_name(keyword_col)

        required = [id_col, go_col, keyword_col]
        missing = [c for c in required if c not in fields]
        if missing:
            raise ValueError(
                "Missing required column(s): "
                + ", ".join(missing)
                + f". Available columns: {', '.join(reader.fieldnames)}"
            )

        if gene_name_col not in fields:
            eprint(
                f"WARNING: gene name column '{gene_name_col}' was not found in the annotation TSV. "
                "Gene-level genus-richness labels will use IDs only."
            )

        for row in reader:
            gene_id = str(row.get(id_col, "")).strip()
            if not gene_id or gene_id in MISSING_VALUES:
                continue

            universe_ids.add(gene_id)

            if gene_name_col in fields:
                gene_name = str(row.get(gene_name_col, "")).strip()
                if gene_name and gene_name not in MISSING_VALUES:
                    # Keep the first non-empty name found for duplicated IDs.
                    gene_name_map.setdefault(gene_id, gene_name)

            for go in split_terms(row.get(go_col), term_delimiter_regex):
                go_term_to_ids[go].add(gene_id)
            for kw in split_terms(row.get(keyword_col), term_delimiter_regex):
                keyword_term_to_ids[kw].add(gene_id)

    if not universe_ids:
        raise ValueError("No valid IDs were found in the annotation TSV.")

    return universe_ids, dict(go_term_to_ids), dict(keyword_term_to_ids), gene_name_map


def extract_genus_from_species(species_name: str) -> Optional[str]:
    """
    Extract a genus-like label from a species string.

    Rules:
      - Trim whitespace and surrounding punctuation.
      - Handle 'Candidatus <Genus> ...' by returning '<Genus>'.
      - Otherwise use the first token.
    """
    if species_name is None:
        return None
    s = str(species_name).strip()
    if not s or s in MISSING_VALUES:
        return None
    s = re.sub(r"\s+", " ", s)
    s = s.strip("[](){}\"'")
    parts = s.split(" ")
    if not parts:
        return None
    if parts[0].lower() == "candidatus" and len(parts) >= 2:
        genus = parts[1]
    else:
        genus = parts[0]
    genus = genus.strip("[](){}\"',;:")
    return genus or None


def read_taxa_terms(
    taxa_tsv: str,
    id_col: str,
    species_col: str,
    sep: str,
) -> Tuple[Set[str], Dict[str, Set[str]], Dict[str, Set[str]], Dict[str, Dict[str, int]]]:
    """
    Read optional taxa hit file.

    Returns:
      taxa_universe_ids
      genus_to_ids
      species_to_ids
      gene_to_genus_hit_counts, where each value is genus -> hit row count
    """
    universe_ids: Set[str] = set()
    genus_to_ids: Dict[str, Set[str]] = defaultdict(set)
    species_to_ids: Dict[str, Set[str]] = defaultdict(set)
    gene_to_genus_hit_counts: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))

    with open(taxa_tsv, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=sep)
        if reader.fieldnames is None:
            raise ValueError("The taxa TSV appears to have no header.")

        reader.fieldnames = [clean_header_name(x) for x in reader.fieldnames]
        fields = set(reader.fieldnames)
        id_col = clean_header_name(id_col)
        species_col = clean_header_name(species_col)

        missing = [c for c in [id_col, species_col] if c not in fields]
        if missing:
            raise ValueError(
                "Missing required taxa column(s): "
                + ", ".join(missing)
                + f". Available columns: {', '.join(reader.fieldnames)}"
            )

        for row in reader:
            gene_id = str(row.get(id_col, "")).strip()
            species = str(row.get(species_col, "")).strip()
            if not gene_id or gene_id in MISSING_VALUES:
                continue

            universe_ids.add(gene_id)

            if species and species not in MISSING_VALUES:
                species_to_ids[species].add(gene_id)
                genus = extract_genus_from_species(species)
                if genus:
                    genus_to_ids[genus].add(gene_id)
                    gene_to_genus_hit_counts[gene_id][genus] += 1

    if not universe_ids:
        raise ValueError("No valid IDs were found in the taxa TSV.")

    return (
        universe_ids,
        dict(genus_to_ids),
        dict(species_to_ids),
        {gene_id: dict(counts) for gene_id, counts in gene_to_genus_hit_counts.items()},
    )


def benjamini_hochberg(pvalues: Sequence[float]) -> List[float]:
    """Return Benjamini-Hochberg adjusted p-values preserving input order."""
    m = len(pvalues)
    if m == 0:
        return []

    indexed = sorted(enumerate(pvalues), key=lambda x: x[1])
    adjusted = [1.0] * m
    prev = 1.0

    for rank_from_end, (idx, pval) in enumerate(reversed(indexed), start=1):
        rank = m - rank_from_end + 1
        qval = min(prev, (float(pval) * m) / rank)
        qval = min(max(qval, 0.0), 1.0)
        adjusted[idx] = qval
        prev = qval

    return adjusted


def make_label(term_id: str, term_name: Optional[str]) -> str:
    if term_name and term_name != term_id:
        return f"{term_id} {term_name}"
    return term_id


def run_enrichment(
    term_to_ids: Dict[str, Set[str]],
    study_ids: Set[str],
    universe_ids: Set[str],
    category: str,
    term_names: Optional[Dict[str, str]] = None,
    term_namespaces: Optional[Dict[str, str]] = None,
    min_term_size: int = 1,
    max_term_size: Optional[int] = None,
) -> pd.DataFrame:
    """Run one-sided hypergeometric enrichment for a term collection."""
    term_names = term_names or {}
    term_namespaces = term_namespaces or {}

    N = len(universe_ids)
    study_present = set(study_ids) & universe_ids
    n = len(study_present)
    if n == 0:
        raise ValueError(f"None of the query IDs were found in the {category} universe.")

    records: List[Dict[str, object]] = []
    tested_pvalues: List[float] = []
    tested_indices: List[int] = []

    for term_id, ids in term_to_ids.items():
        term_universe_ids = set(ids) & universe_ids
        K = len(term_universe_ids)

        if K < min_term_size:
            continue
        if max_term_size is not None and K > max_term_size:
            continue

        hits = term_universe_ids & study_present
        k = len(hits)

        # P(X >= k) for overlap between study genes and term genes.
        pvalue = float(hypergeom.sf(k - 1, N, K, n)) if k > 0 else 1.0
        expected = (n * K / N) if N else float("nan")
        fold_enrichment = (k / expected) if expected > 0 else float("nan")

        term_name = term_names.get(term_id, term_id if category in {"KEYWORD", "SPECIES"} else "")
        namespace = term_namespaces.get(term_id, "unknown") if category == "GO" else category

        records.append({
            "category": category,
            "namespace": namespace,
            "term_id": term_id,
            "term_name": term_name,
            "label": make_label(term_id, term_name),
            "query_count": k,
            "query_size": n,
            "universe_count": K,
            "universe_size": N,
            "expected_count": expected,
            "fold_enrichment": fold_enrichment,
            "pvalue": pvalue,
            "fdr_bh": 1.0,
            "query_gene_ids": ";".join(sorted(hits)),
        })
        tested_indices.append(len(records) - 1)
        tested_pvalues.append(pvalue)

    adjusted = benjamini_hochberg(tested_pvalues)
    for idx, qval in zip(tested_indices, adjusted):
        records[idx]["fdr_bh"] = qval

    columns = [
        "category", "namespace", "term_id", "term_name", "label",
        "query_count", "query_size", "universe_count", "universe_size",
        "expected_count", "fold_enrichment", "pvalue", "fdr_bh",
        "query_gene_ids",
    ]

    df = pd.DataFrame.from_records(records)
    if df.empty:
        return pd.DataFrame(columns=columns)

    # Keep only terms observed in the query list for reporting, after FDR correction.
    df = df[df["query_count"] > 0].copy()
    if df.empty:
        return pd.DataFrame(columns=columns)

    df = df.sort_values(["fdr_bh", "pvalue", "fold_enrichment"], ascending=[True, True, False])
    return df[columns]


def add_namespace_specific_fdr(go_df: pd.DataFrame) -> pd.DataFrame:
    """Add FDR corrected within each GO namespace."""
    if go_df.empty or "namespace" not in go_df.columns:
        go_df = go_df.copy()
        go_df["fdr_bh_namespace"] = []
        return go_df

    go_df = go_df.copy()
    go_df["fdr_bh_namespace"] = 1.0

    for _, idx in go_df.groupby("namespace").groups.items():
        pvals = go_df.loc[idx, "pvalue"].astype(float).tolist()
        go_df.loc[idx, "fdr_bh_namespace"] = benjamini_hochberg(pvals)

    cols = list(go_df.columns)
    cols.remove("fdr_bh_namespace")
    insert_at = cols.index("fdr_bh") + 1 if "fdr_bh" in cols else len(cols)
    cols.insert(insert_at, "fdr_bh_namespace")
    return go_df[cols]


def build_genus_abundance_table(
    query_ids: Set[str],
    taxa_universe_ids: Set[str],
    genus_to_ids: Dict[str, Set[str]],
    gene_to_genus_hit_counts: Dict[str, Dict[str, int]],
) -> pd.DataFrame:
    """Create a genus abundance table for query genes, without enrichment statistics."""
    query_present = set(query_ids) & set(taxa_universe_ids)
    total_query_hits = sum(
        sum(counts.values())
        for gene_id, counts in gene_to_genus_hit_counts.items()
        if gene_id in query_present
    )
    total_universe_hits = sum(sum(counts.values()) for counts in gene_to_genus_hit_counts.values())

    records: List[Dict[str, object]] = []

    for genus, ids in genus_to_ids.items():
        universe_gene_ids = set(ids) & set(taxa_universe_ids)
        query_gene_ids = universe_gene_ids & query_present
        if not query_gene_ids:
            continue

        query_hit_count = sum(gene_to_genus_hit_counts.get(gene_id, {}).get(genus, 0) for gene_id in query_gene_ids)
        universe_hit_count = sum(gene_to_genus_hit_counts.get(gene_id, {}).get(genus, 0) for gene_id in universe_gene_ids)

        records.append({
            "genus": genus,
            "query_gene_count": len(query_gene_ids),
            "query_hit_count": query_hit_count,
            "query_gene_fraction": (len(query_gene_ids) / len(query_present)) if query_present else float("nan"),
            "query_hit_fraction": (query_hit_count / total_query_hits) if total_query_hits > 0 else float("nan"),
            "universe_gene_count": len(universe_gene_ids),
            "universe_hit_count": universe_hit_count,
            "universe_gene_fraction": (len(universe_gene_ids) / len(taxa_universe_ids)) if taxa_universe_ids else float("nan"),
            "universe_hit_fraction": (universe_hit_count / total_universe_hits) if total_universe_hits > 0 else float("nan"),
            "query_gene_ids": ";".join(sorted(query_gene_ids)),
        })

    columns = [
        "genus", "query_gene_count", "query_hit_count", "query_gene_fraction",
        "query_hit_fraction", "universe_gene_count", "universe_hit_count",
        "universe_gene_fraction", "universe_hit_fraction", "query_gene_ids",
    ]

    df = pd.DataFrame.from_records(records)
    if df.empty:
        return pd.DataFrame(columns=columns)

    df = df.sort_values(
        ["query_gene_count", "query_hit_count", "universe_gene_count", "genus"],
        ascending=[False, False, False, True],
    )
    return df[columns]


def build_gene_genus_richness_table(
    query_ids: Set[str],
    taxa_universe_ids: Set[str],
    gene_to_genus_hit_counts: Dict[str, Dict[str, int]],
    gene_name_map: Optional[Dict[str, str]] = None,
) -> pd.DataFrame:
    """Create a gene-level table showing how many distinct genera hit each query gene."""
    gene_name_map = gene_name_map or {}
    query_present = set(query_ids) & set(taxa_universe_ids)

    records: List[Dict[str, object]] = []

    for gene_id in sorted(query_present):
        genus_counts = gene_to_genus_hit_counts.get(gene_id, {})
        if not genus_counts:
            continue

        ordered = sorted(genus_counts.items(), key=lambda x: (-x[1], x[0]))
        total_hits = sum(genus_counts.values())

        gene_name = gene_name_map.get(gene_id, "")
        display_label = f"{gene_name} ({gene_id})" if gene_name and gene_name != gene_id else gene_id

        records.append({
            "gene_id": gene_id,
            "gene_name": gene_name,
            "display_label": display_label,
            "distinct_genus_count": len(genus_counts),
            "total_hit_count": total_hits,
            "top_genus": ordered[0][0],
            "genus_hit_counts": ";".join(f"{genus}({count})" for genus, count in ordered),
        })

    columns = [
        "gene_id", "gene_name", "display_label", "distinct_genus_count",
        "total_hit_count", "top_genus", "genus_hit_counts",
    ]

    df = pd.DataFrame.from_records(records)
    if df.empty:
        return pd.DataFrame(columns=columns)

    df = df.sort_values(["distinct_genus_count", "total_hit_count", "gene_id"], ascending=[False, False, True])
    return df[columns]


def truncate_label(label: str, max_len: int = 90) -> str:
    label = str(label)
    return label if len(label) <= max_len else label[: max_len - 1] + "â€¦"


def safe_neglog10(values: Iterable[float]) -> List[float]:
    out: List[float] = []
    for v in values:
        try:
            v = float(v)
        except Exception:
            v = 1.0
        if not math.isfinite(v) or v <= 0:
            v = 1e-300
        out.append(-math.log10(v))
    return out


def bubble_size(count: float) -> float:
    """Convert a query gene count into a bubble area for dotplots."""
    try:
        count = float(count)
    except Exception:
        count = 1.0
    if not math.isfinite(count) or count <= 0:
        count = 1.0
    # Square-root scaling keeps large counts visible without overwhelming the plot.
    return 35.0 + 38.0 * math.sqrt(count)


def choose_fdr_column(df: pd.DataFrame, prefer_namespace_fdr: bool = False) -> str:
    if prefer_namespace_fdr and "fdr_bh_namespace" in df.columns:
        return "fdr_bh_namespace"
    return "fdr_bh"


def write_empty_plot(message: str, out_base: str, formats: Sequence[str]) -> None:
    for fmt in formats:
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=11)
        ax.axis("off")
        fig.tight_layout()
        fig.savefig(f"{out_base}.{fmt}", dpi=300, bbox_inches="tight")
        plt.close(fig)


def make_dotplot(
    df: pd.DataFrame,
    title: str,
    out_base: str,
    top_n: int = 20,
    alpha: float = 0.05,
    formats: Sequence[str] = ("png", "pdf"),
    prefer_namespace_fdr: bool = False,
) -> None:
    if df.empty:
        write_empty_plot(f"No enriched terms to plot: {title}", out_base, formats)
        return

    fdr_col = choose_fdr_column(df, prefer_namespace_fdr=prefer_namespace_fdr)
    plot_df = df.copy()
    sig = plot_df[plot_df[fdr_col] <= alpha]
    if not sig.empty:
        plot_df = sig

    plot_df = plot_df.sort_values([fdr_col, "pvalue", "fold_enrichment"], ascending=[True, True, False]).head(top_n)
    plot_df = plot_df.iloc[::-1].copy()

    labels = [truncate_label(x) for x in plot_df["label"]]
    x = plot_df["fold_enrichment"].astype(float).tolist()
    y = list(range(len(plot_df)))
    sizes = [bubble_size(v) for v in plot_df["query_count"]]
    color_values = safe_neglog10(plot_df[fdr_col])

    height = max(4.0, 0.35 * len(plot_df) + 1.8)
    width = 11.0

    for fmt in formats:
        fig, ax = plt.subplots(figsize=(width, height))
        sc = ax.scatter(x, y, s=sizes, c=color_values, alpha=0.85, edgecolors="black", linewidths=0.4)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=9)
        ax.set_xlabel("Fold enrichment")
        ax.set_title(title)
        ax.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.5)

        cbar = fig.colorbar(sc, ax=ax)
        if fdr_col == "fdr_bh_namespace":
            cbar.set_label("Significance: -log10(namespace FDR)")
        else:
            cbar.set_label("Significance: -log10(FDR) (higher = more significant; 1.30 â‰ˆ FDR 0.05)")

        handles = [
            ax.scatter([], [], s=bubble_size(c), c="none", edgecolors="black", label=str(c))
            for c in FIXED_BUBBLE_LEGEND_COUNTS
        ]
        ax.legend(handles=handles, title="Query genes", frameon=True, loc="lower right")

        fig.tight_layout()
        fig.savefig(f"{out_base}.{fmt}", dpi=300, bbox_inches="tight")
        plt.close(fig)


def make_barplot(
    df: pd.DataFrame,
    title: str,
    out_base: str,
    top_n: int = 20,
    alpha: float = 0.05,
    formats: Sequence[str] = ("png", "pdf"),
    prefer_namespace_fdr: bool = False,
) -> None:
    if df.empty:
        write_empty_plot(f"No enriched terms to plot: {title}", out_base, formats)
        return

    fdr_col = choose_fdr_column(df, prefer_namespace_fdr=prefer_namespace_fdr)
    plot_df = df.copy()
    sig = plot_df[plot_df[fdr_col] <= alpha]
    if not sig.empty:
        plot_df = sig

    plot_df = plot_df.sort_values([fdr_col, "pvalue", "query_count"], ascending=[True, True, False]).head(top_n)
    plot_df = plot_df.iloc[::-1].copy()

    labels = [truncate_label(x) for x in plot_df["label"]]
    values = safe_neglog10(plot_df[fdr_col])
    y = list(range(len(plot_df)))

    height = max(4.0, 0.35 * len(plot_df) + 1.8)
    width = 11.0

    for fmt in formats:
        fig, ax = plt.subplots(figsize=(width, height))
        ax.barh(y, values)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=9)
        ax.set_xlabel("-log10(namespace FDR)" if fdr_col == "fdr_bh_namespace" else "-log10(FDR)")
        ax.set_title(title)
        ax.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.5)

        if alpha > 0 and alpha < 1:
            ax.axvline(-math.log10(alpha), linestyle="--", linewidth=1.0, alpha=0.8)

        for yi, value, count in zip(y, values, plot_df["query_count"]):
            ax.text(value, yi, f"  n={int(count)}", va="center", fontsize=8)

        fig.tight_layout()
        fig.savefig(f"{out_base}.{fmt}", dpi=300, bbox_inches="tight")
        plt.close(fig)


def make_treemap(
    df: pd.DataFrame,
    title: str,
    out_base: str,
    top_n: int = 25,
    alpha: float = 0.05,
    formats: Sequence[str] = ("png",),
    size_by: str = "query_count",
) -> None:
    if not HAS_SQUARIFY:
        write_empty_plot("Treemap requires squarify. Install it to enable this plot.", out_base, formats)
        return
    if df.empty:
        write_empty_plot(f"No enriched terms to plot: {title}", out_base, formats)
        return

    plot_df = df.copy()
    sig = plot_df[plot_df["fdr_bh"] <= alpha]
    if not sig.empty:
        plot_df = sig

    plot_df = plot_df.sort_values(["fdr_bh", "pvalue", "fold_enrichment"], ascending=[True, True, False]).head(top_n)
    if plot_df.empty:
        write_empty_plot(f"No enriched terms to plot: {title}", out_base, formats)
        return

    sizes = plot_df[size_by].astype(float).tolist()
    labels: List[str] = []
    for _, row in plot_df.iterrows():
        labels.append(
            f"{truncate_label(row['label'], 42)}\n"
            f"n={int(row['query_count'])}\n"
            f"FE={float(row['fold_enrichment']):.2f}\n"
            f"FDR={float(row['fdr_bh']):.2e}"
        )

    color_values = safe_neglog10(plot_df["fdr_bh"])
    min_color = min(color_values)
    max_color = max(color_values)
    if max_color <= min_color:
        max_color = min_color + 1.0

    for fmt in formats:
        fig, ax = plt.subplots(figsize=(12, 8))
        norm = plt.Normalize(min_color, max_color)
        colors = plt.cm.viridis(norm(color_values))
        squarify.plot(
            sizes=sizes,
            label=labels,
            color=colors,
            alpha=0.9,
            pad=True,
            ax=ax,
            text_kwargs={"fontsize": 8},
        )
        ax.axis("off")
        ax.set_title(title)

        sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.01)
        cbar.set_label("-log10(FDR)")

        fig.tight_layout()
        fig.savefig(f"{out_base}.{fmt}", dpi=300, bbox_inches="tight")
        plt.close(fig)


def make_genus_abundance_barplot(
    df: pd.DataFrame,
    title: str,
    out_base: str,
    top_n: int = 20,
    formats: Sequence[str] = ("png", "pdf"),
) -> None:
    if df.empty:
        write_empty_plot(f"No genera found among query genes: {title}", out_base, formats)
        return

    plot_df = df.head(top_n).iloc[::-1].copy()
    labels = [truncate_label(x, 50) for x in plot_df["genus"]]
    values = plot_df["query_gene_count"].astype(int).tolist()
    hit_counts = plot_df["query_hit_count"].astype(int).tolist()
    y = list(range(len(plot_df)))

    height = max(4.0, 0.35 * len(plot_df) + 1.8)
    width = 11.0

    for fmt in formats:
        fig, ax = plt.subplots(figsize=(width, height))
        ax.barh(y, values)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=9)
        ax.set_xlabel("Query genes with hits to genus")
        ax.set_title(title)
        ax.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.5)

        for yi, value, hits in zip(y, values, hit_counts):
            ax.text(value, yi, f"  hit_rows={hits}", va="center", fontsize=8)

        fig.tight_layout()
        fig.savefig(f"{out_base}.{fmt}", dpi=300, bbox_inches="tight")
        plt.close(fig)


def make_gene_genus_richness_barplot(
    df: pd.DataFrame,
    title: str,
    out_base: str,
    top_n: int = 20,
    formats: Sequence[str] = ("png", "pdf"),
) -> None:
    if df.empty:
        write_empty_plot(f"No gene/genus hit diversity data to plot: {title}", out_base, formats)
        return

    plot_df = df.head(top_n).iloc[::-1].copy()
    label_column = "display_label" if "display_label" in plot_df.columns else "gene_id"
    labels = [truncate_label(x, 75) for x in plot_df[label_column]]
    values = plot_df["distinct_genus_count"].astype(int).tolist()
    total_hits = plot_df["total_hit_count"].astype(int).tolist()
    y = list(range(len(plot_df)))

    # Wider figure because labels can be long: GENENAME (ID).
    height = max(4.5, 0.38 * len(plot_df) + 1.8)
    width = 13.5

    for fmt in formats:
        fig, ax = plt.subplots(figsize=(width, height))
        ax.barh(y, values)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=9)
        ax.set_xlabel("Distinct genera per query gene")
        ax.set_title(title)
        ax.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.5)

        for yi, value, hits in zip(y, values, total_hits):
            ax.text(value, yi, f"  total_hit_rows={hits}", va="center", fontsize=8)

        fig.tight_layout()
        fig.savefig(f"{out_base}.{fmt}", dpi=300, bbox_inches="tight")
        plt.close(fig)


def write_summary(
    outpath: str,
    annotations: str,
    query_ids_path: str,
    universe_ids: Set[str],
    query_ids: Set[str],
    go_df: pd.DataFrame,
    kw_df: pd.DataFrame,
    alpha: float,
    min_term_size: int,
    max_term_size: Optional[int],
    go_names_loaded: int,
    go_namespaces_loaded: int,
    gene_names_loaded: int,
    genus_abundance_df: Optional[pd.DataFrame] = None,
    gene_genus_richness_df: Optional[pd.DataFrame] = None,
    taxa_universe_ids: Optional[Set[str]] = None,
    taxa_hits_path: Optional[str] = None,
) -> None:
    query_present = query_ids & universe_ids
    missing = query_ids - universe_ids

    with open(outpath, "wt", encoding="utf-8") as out:
        out.write("Functional enrichment summary\n")
        out.write("=============================\n")
        out.write(f"Annotation universe: {annotations}\n")
        out.write(f"Query IDs: {query_ids_path}\n")
        out.write(f"Universe size: {len(universe_ids)}\n")
        out.write(f"Query IDs supplied: {len(query_ids)}\n")
        out.write(f"Query IDs found in annotation universe: {len(query_present)}\n")
        out.write(f"Query IDs missing from annotation universe: {len(missing)}\n")
        out.write(f"Alpha/FDR cutoff used for plots: {alpha}\n")
        out.write(f"Minimum term size: {min_term_size}\n")
        out.write(f"Maximum term size: {max_term_size if max_term_size is not None else 'None'}\n")
        out.write(f"GO names loaded: {go_names_loaded}\n")
        out.write(f"GO namespaces loaded: {go_namespaces_loaded}\n")
        out.write(f"Gene names loaded from annotation TSV: {gene_names_loaded}\n\n")

        out.write(f"GO terms with >=1 query hit: {len(go_df)}\n")
        out.write(
            f"GO terms significant at global FDR <= {alpha}: "
            f"{int((go_df['fdr_bh'] <= alpha).sum()) if not go_df.empty else 0}\n"
        )

        if not go_df.empty and "namespace" in go_df.columns:
            out.write("\nGO namespace summary\n")
            for namespace in sorted(go_df["namespace"].dropna().unique()):
                sub = go_df[go_df["namespace"] == namespace]
                n_sig_global = int((sub["fdr_bh"] <= alpha).sum()) if "fdr_bh" in sub.columns else 0
                n_sig_ns = int((sub["fdr_bh_namespace"] <= alpha).sum()) if "fdr_bh_namespace" in sub.columns else 0
                out.write(
                    f"  {namespace}: terms_with_hits={len(sub)}, "
                    f"significant_global_FDR={n_sig_global}, "
                    f"significant_namespace_FDR={n_sig_ns}\n"
                )

        out.write("\n")
        out.write(f"KEYWORD terms with >=1 query hit: {len(kw_df)}\n")
        out.write(
            f"KEYWORD terms significant at FDR <= {alpha}: "
            f"{int((kw_df['fdr_bh'] <= alpha).sum()) if not kw_df.empty else 0}\n"
        )

        if taxa_hits_path is not None:
            taxa_universe_ids = taxa_universe_ids or set()
            taxa_query_present = query_ids & taxa_universe_ids
            out.write("\nTaxa/genus abundance summary\n")
            out.write(f"Taxa hits file: {taxa_hits_path}\n")
            out.write(f"Taxa universe size: {len(taxa_universe_ids)}\n")
            out.write(f"Query IDs found in taxa universe: {len(taxa_query_present)}\n")
            if genus_abundance_df is not None:
                out.write(f"Genera observed among query genes: {len(genus_abundance_df)}\n")
            if gene_genus_richness_df is not None:
                out.write(f"Query genes with at least one genus hit: {len(gene_genus_richness_df)}\n")

        if missing:
            out.write("\nMissing query IDs in annotation universe, first 100:\n")
            for gene_id in sorted(missing)[:100]:
                out.write(f"{gene_id}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="GO / KEYWORD enrichment with optional GO namespaces, keyword treemaps, and genus abundance summaries."
    )
    parser.add_argument("-a", "--annotations", required=True,
                        help="Annotation TSV used as the GO/KEYWORD universe.")
    parser.add_argument("-i", "--ids", required=True,
                        help="Text file with query IDs, one per line. If there are several columns, the first one is used.")
    parser.add_argument("-o", "--out-prefix", default="functional_enrichment",
                        help="Output file prefix inside the output directory. Default: functional_enrichment")
    parser.add_argument("-d", "--outdir", default="enrichment",
                        help="Output directory where all tables, plots and summary files will be written. Default: enrichment")

    parser.add_argument("--id-column", default="#ID",
                        help="ID column in annotation TSV. Default: #ID")
    parser.add_argument("--genename-column", default="GENENAME",
                        help="Gene name column in annotation TSV. Used for labels like GENENAME (ID). Default: GENENAME")
    parser.add_argument("--go-column", default="GO",
                        help="GO column in annotation TSV. Default: GO")
    parser.add_argument("--keyword-column", default="KEYWORD",
                        help="KEYWORD column in annotation TSV. Default: KEYWORD")
    parser.add_argument("--sep", default="\t",
                        help="Input separator for annotation and taxa TSV files. Default: tab")
    parser.add_argument("--term-delimiter-regex", default=r";",
                        help="Regex used to split GO/KEYWORD values. Default: ';'")

    parser.add_argument("--go-obo", default=None,
                        help="Optional go-basic.obo file to map GO IDs to names and namespaces.")
    parser.add_argument("--go-names-tsv", default=None,
                        help="Optional TSV/CSV mapping GO IDs to names. Example: GO:0005525<TAB>GTP binding")
    parser.add_argument("--no-split-go-namespaces", action="store_true",
                        help="Do not write separate GO namespace tables/plots, even if --go-obo is provided.")
    parser.add_argument("--use-namespace-fdr-for-namespace-plots", action="store_true",
                        help="For BP/MF/CC plots, use FDR corrected only within that namespace.")

    parser.add_argument("--taxa-hits", default=None,
                        help="Optional TSV with columns like ID ... HIT_SPECIES ... to summarize genus abundance and gene-level genus diversity.")
    parser.add_argument("--taxa-id-column", default="ID",
                        help="ID column in --taxa-hits. Default: ID")
    parser.add_argument("--taxa-species-column", default="HIT_SPECIES",
                        help="Species column in --taxa-hits. Default: HIT_SPECIES")
    parser.add_argument("--write-species-table", action="store_true",
                        help="Also write a SPECIES enrichment table from --taxa-hits. No default species plots are produced.")

    parser.add_argument("--make-keyword-treemap", action="store_true",
                        help="Write a treemap for KEYWORD enrichment.")

    parser.add_argument("--min-term-size", type=int, default=1,
                        help="Ignore terms annotated to fewer than this many universe IDs. Default: 1")
    parser.add_argument("--max-term-size", type=int, default=None,
                        help="Ignore terms annotated to more than this many universe IDs. Default: no limit")
    parser.add_argument("--alpha", type=float, default=0.05,
                        help="FDR cutoff used to prioritize terms in plots. Default: 0.05")
    parser.add_argument("--top", type=int, default=20,
                        help="Maximum number of terms shown in each plot. Default: 20")
    parser.add_argument("--plot-formats", default="png,pdf",
                        help="Comma-separated plot formats. Default: png,pdf")

    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if args.min_term_size < 1:
        raise ValueError("--min-term-size must be >= 1")
    if args.max_term_size is not None and args.max_term_size < args.min_term_size:
        raise ValueError("--max-term-size must be >= --min-term-size")

    plot_formats = [x.strip().lower() for x in args.plot_formats.split(",") if x.strip()]
    if not plot_formats:
        plot_formats = ["png"]

    os.makedirs(args.outdir, exist_ok=True)
    out_prefix = os.path.join(args.outdir, args.out_prefix)

    eprint("[1/6] Reading query IDs...")
    query_ids = read_query_ids(args.ids)
    if not query_ids:
        raise ValueError("No query IDs were found in the ID list.")

    eprint("[2/6] Reading annotation universe and building GO/KEYWORD mappings...")
    universe_ids, go_term_to_ids, kw_term_to_ids, gene_name_map = read_annotation_data(
        annotations_tsv=args.annotations,
        id_col=args.id_column,
        gene_name_col=args.genename_column,
        go_col=args.go_column,
        keyword_col=args.keyword_column,
        sep=args.sep,
        term_delimiter_regex=args.term_delimiter_regex,
    )

    query_present = query_ids & universe_ids
    if not query_present:
        raise ValueError("None of the query IDs are present in the annotation universe.")
    if len(query_present) < len(query_ids):
        eprint(f"WARNING: {len(query_ids) - len(query_present)} query IDs were not found in the annotation universe.")

    eprint("[3/6] Loading optional GO names and namespaces...")
    go_names: Dict[str, str] = {}
    go_namespaces: Dict[str, str] = {}

    if args.go_obo:
        obo_names, obo_namespaces = parse_go_obo_metadata(args.go_obo)
        go_names.update(obo_names)
        go_namespaces.update(obo_namespaces)

    if args.go_names_tsv:
        go_names.update(read_go_names_tsv(args.go_names_tsv))

    eprint("[4/6] Running GO/KEYWORD enrichment...")
    go_df = run_enrichment(
        term_to_ids=go_term_to_ids,
        study_ids=query_ids,
        universe_ids=universe_ids,
        category="GO",
        term_names=go_names,
        term_namespaces=go_namespaces,
        min_term_size=args.min_term_size,
        max_term_size=args.max_term_size,
    )
    go_df = add_namespace_specific_fdr(go_df)

    kw_df = run_enrichment(
        term_to_ids=kw_term_to_ids,
        study_ids=query_ids,
        universe_ids=universe_ids,
        category="KEYWORD",
        term_names=None,
        min_term_size=args.min_term_size,
        max_term_size=args.max_term_size,
    )

    genus_abundance_df: Optional[pd.DataFrame] = None
    gene_genus_richness_df: Optional[pd.DataFrame] = None
    species_df: Optional[pd.DataFrame] = None
    taxa_universe_ids: Set[str] = set()

    if args.taxa_hits:
        eprint("[5/6] Reading optional taxa file and summarizing genus abundance...")
        taxa_universe_ids, genus_to_ids, species_to_ids, gene_to_genus_hit_counts = read_taxa_terms(
            taxa_tsv=args.taxa_hits,
            id_col=args.taxa_id_column,
            species_col=args.taxa_species_column,
            sep=args.sep,
        )

        taxa_query_present = query_ids & taxa_universe_ids
        if not taxa_query_present:
            eprint("WARNING: None of the query IDs are present in the taxa universe. Genus abundance outputs will be skipped.")
        else:
            genus_abundance_df = build_genus_abundance_table(
                query_ids=query_ids,
                taxa_universe_ids=taxa_universe_ids,
                genus_to_ids=genus_to_ids,
                gene_to_genus_hit_counts=gene_to_genus_hit_counts,
            )
            gene_genus_richness_df = build_gene_genus_richness_table(
                query_ids=query_ids,
                taxa_universe_ids=taxa_universe_ids,
                gene_to_genus_hit_counts=gene_to_genus_hit_counts,
                gene_name_map=gene_name_map,
            )
            if args.write_species_table:
                species_df = run_enrichment(
                    term_to_ids=species_to_ids,
                    study_ids=query_ids,
                    universe_ids=taxa_universe_ids,
                    category="SPECIES",
                    term_names=None,
                    min_term_size=args.min_term_size,
                    max_term_size=args.max_term_size,
                )
    else:
        eprint("[5/6] No taxa file provided; skipping genus abundance summaries.")

    eprint("[6/6] Writing tables, plots and summary...")

    go_out = f"{out_prefix}.GO.enrichment.tsv"
    kw_out = f"{out_prefix}.KEYWORD.enrichment.tsv"
    go_df.to_csv(go_out, sep="\t", index=False)
    kw_df.to_csv(kw_out, sep="\t", index=False)

    if genus_abundance_df is not None:
        genus_abundance_df.to_csv(f"{out_prefix}.GENUS.abundance.tsv", sep="\t", index=False)

    if gene_genus_richness_df is not None:
        gene_genus_richness_df.to_csv(f"{out_prefix}.GENE.genus_richness.tsv", sep="\t", index=False)

    if species_df is not None:
        species_df.to_csv(f"{out_prefix}.SPECIES.enrichment.tsv", sep="\t", index=False)

    namespace_tables: List[Tuple[str, pd.DataFrame]] = []
    if not args.no_split_go_namespaces:
        if args.go_obo:
            namespaces_to_write = GO_NAMESPACES + sorted(
                ns for ns in go_df["namespace"].dropna().unique().tolist()
                if ns not in GO_NAMESPACES
            ) if not go_df.empty else GO_NAMESPACES

            for namespace in namespaces_to_write:
                sub = go_df[go_df["namespace"] == namespace].copy() if not go_df.empty else go_df.copy()
                namespace_tables.append((namespace, sub))
                safe_ns = namespace.replace(" ", "_")
                sub.to_csv(f"{out_prefix}.GO.{safe_ns}.enrichment.tsv", sep="\t", index=False)
        else:
            eprint("WARNING: GO namespaces require --go-obo. Namespace-specific GO outputs were not created.")

    make_dotplot(go_df, "GO term enrichment", f"{out_prefix}.GO.dotplot", args.top, args.alpha, plot_formats)
    make_barplot(go_df, "GO term enrichment", f"{out_prefix}.GO.barplot", args.top, args.alpha, plot_formats)

    for namespace, sub in namespace_tables:
        safe_ns = namespace.replace(" ", "_")
        title = GO_NAMESPACE_LABELS.get(namespace, f"GO {namespace}")
        make_dotplot(
            sub,
            f"{title} enrichment",
            f"{out_prefix}.GO.{safe_ns}.dotplot",
            args.top,
            args.alpha,
            plot_formats,
            prefer_namespace_fdr=args.use_namespace_fdr_for_namespace_plots,
        )
        make_barplot(
            sub,
            f"{title} enrichment",
            f"{out_prefix}.GO.{safe_ns}.barplot",
            args.top,
            args.alpha,
            plot_formats,
            prefer_namespace_fdr=args.use_namespace_fdr_for_namespace_plots,
        )

    make_dotplot(kw_df, "KEYWORD enrichment", f"{out_prefix}.KEYWORD.dotplot", args.top, args.alpha, plot_formats)
    make_barplot(kw_df, "KEYWORD enrichment", f"{out_prefix}.KEYWORD.barplot", args.top, args.alpha, plot_formats)

    if args.make_keyword_treemap:
        treemap_formats = [f for f in plot_formats if f in {"png", "pdf", "svg"}] or ["png"]
        make_treemap(
            kw_df,
            "KEYWORD enrichment treemap",
            f"{out_prefix}.KEYWORD.treemap",
            top_n=max(args.top, 25),
            alpha=args.alpha,
            formats=treemap_formats,
        )

    if genus_abundance_df is not None:
        make_genus_abundance_barplot(
            genus_abundance_df,
            "Genus abundance among query genes",
            f"{out_prefix}.GENUS.abundance",
            top_n=args.top,
            formats=plot_formats,
        )

    if gene_genus_richness_df is not None:
        make_gene_genus_richness_barplot(
            gene_genus_richness_df,
            "Query genes with the highest number of distinct genus hits",
            f"{out_prefix}.GENE.genus_richness",
            top_n=args.top,
            formats=plot_formats,
        )

    write_summary(
        outpath=f"{out_prefix}.summary.txt",
        annotations=args.annotations,
        query_ids_path=args.ids,
        universe_ids=universe_ids,
        query_ids=query_ids,
        go_df=go_df,
        kw_df=kw_df,
        alpha=args.alpha,
        min_term_size=args.min_term_size,
        max_term_size=args.max_term_size,
        go_names_loaded=len(go_names),
        go_namespaces_loaded=len(go_namespaces),
        gene_names_loaded=len(gene_name_map),
        genus_abundance_df=genus_abundance_df,
        gene_genus_richness_df=gene_genus_richness_df,
        taxa_universe_ids=taxa_universe_ids,
        taxa_hits_path=args.taxa_hits,
    )

    eprint("Done.")
    eprint(f"Output directory: {args.outdir}")
    eprint(f"GO table: {go_out}")
    eprint(f"KEYWORD table: {kw_out}")
    if genus_abundance_df is not None:
        eprint(f"GENUS abundance table: {out_prefix}.GENUS.abundance.tsv")
    if gene_genus_richness_df is not None:
        eprint(f"Gene genus-richness table: {out_prefix}.GENE.genus_richness.tsv")
    if species_df is not None:
        eprint(f"SPECIES table: {out_prefix}.SPECIES.enrichment.tsv")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        raise SystemExit(1)
    except Exception as exc:
        eprint(f"ERROR: {exc}")
        raise SystemExit(1)

