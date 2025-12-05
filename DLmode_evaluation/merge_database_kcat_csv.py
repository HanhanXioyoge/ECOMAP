#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
merge_database_kcat_csv.py

Merge rows in a CSV according to these rules:

1) Rows are mergeable only if ALL of the following A-keys are exactly equal:
      A_KEYS = ["Organism", "uniprot", "ec", "sequence"]

2) Within each A-group, rows that share ANY exact match among the B-keys
   belong to the same connected component (union-find):
      B_KEYS = ["MetaNetXID", "substrate", "substrate_smiles", "InChIKey"]

   This supports "chain connectivity": if row1 shares 'substrate' with row2,
   and row2 shares 'InChIKey' with row3, then {row1,row2,row3} collapses into
   one merged row, even if row1 and row3 have no direct B-key match.

3) For each connected component:
   - Keep all non-specified columns (including the B-keys) from the earliest
     row (smallest original index) in that component.
   - Aggregate the numeric column 'value' using either arithmetic mean or
     median (user-selectable).

4) Output:
   - Preserve the original column order of the input.
   - Preserve row order by the earliest original row index per merged component.

Usage:
    python merge_database_kcat_csv.py --input database.csv --output out.csv --agg mean
    python merge_database_kcat_csv.py -i database.csv -o out.csv -a median

Dependencies:
    pip install pandas numpy
"""

from __future__ import annotations
import argparse
import contextlib
import os
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

# ---- Configuration keys ----
A_KEYS = ["Organism", "uniprot", "ec", "sequence"]
B_KEYS = ["MetaNetXID", "substrate", "substrate_smiles", "InChIKey"]
VALUE_COL = "value"
AUX_IDX = "__orig_idx__"  # internal helper column for stable ordering


# ---------------- Union-Find (Disjoint Set Union) ----------------
class DSU:
    """Simple Union-Find structure to build connected components."""
    def __init__(self, n: int):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x: int) -> int:
        # Path compression
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        # Union by rank
        if self.rank[ra] < self.rank[rb]:
            self.parent[ra] = rb
        elif self.rank[ra] > self.rank[rb]:
            self.parent[rb] = ra
        else:
            self.parent[rb] = ra
            self.rank[ra] += 1


# ---------------- Utility helpers ----------------
def _normalize_str_col(s: pd.Series) -> pd.Series:
    """
    Normalize a string-like Series for exact matching:
    - cast to string
    - fill NaN with empty string
    - strip whitespace
    """
    return s.astype(str).fillna("").str.strip()


def _build_components_with_any_b_match(sub: pd.DataFrame) -> List[List[int]]:
    """
    Given a subgroup 'sub' where all A-keys are identical,
    connect rows into components if ANY of the B-keys match exactly.

    Returns:
        A list of components, each being a list of GLOBAL row indices.
    """
    m = len(sub)
    if m <= 1:
        return [sub.index.tolist()]  # single row → single component

    # Map local index <-> global index
    local_to_global = list(sub.index)  # global row indices of this subgroup
    dsu = DSU(m)

    # For each B-key, group identical (non-empty) values and union them
    for col in B_KEYS:
        if col not in sub.columns:
            continue

        series = _normalize_str_col(sub[col])
        val_map: Dict[str, List[int]] = {}
        for j, v in enumerate(series.tolist()):
            # Ignore empty strings so they do not connect everything
            if not v:
                continue
            val_map.setdefault(v, []).append(j)

        for local_idxs in val_map.values():
            if len(local_idxs) > 1:
                base = local_idxs[0]
                for other in local_idxs[1:]:
                    dsu.union(base, other)

    # Aggregate members by root
    comp_map: Dict[int, List[int]] = {}
    for j in range(m):
        root = dsu.find(j)
        comp_map.setdefault(root, []).append(j)

    # Convert local indices to GLOBAL row indices and sort members
    components: List[List[int]] = []
    for members in comp_map.values():
        components.append(sorted(local_to_global[idx] for idx in members))

    return components


# ---------------- Core merging logic ----------------
def merge_dataframe(df: pd.DataFrame, agg: str = "mean") -> pd.DataFrame:
    """
    Merge rows in df by:
      - all A_KEYS equal, and
      - ANY of B_KEYS equal (connected components).
    Aggregate 'value' by mean/median. Preserve original column order and
    use the earliest row in each component as representative for non-aggregated fields.

    Args:
        df  : input DataFrame
        agg : "mean" or "median" for the 'value' column

    Returns:
        Merged DataFrame.
    """
    if VALUE_COL not in df.columns:
        raise ValueError(f"Missing required numeric column '{VALUE_COL}' in input CSV.")

    original_cols = list(df.columns)

    # Add a stable original index used for ordering representatives.
    if AUX_IDX in df.columns:
        raise ValueError(f"Input already contains reserved column name '{AUX_IDX}'.")
    df = df.copy()
    df[AUX_IDX] = np.arange(len(df), dtype=int)

    # Normalize string columns used for matching
    for col in A_KEYS + B_KEYS:
        if col in df.columns:
            df[col] = _normalize_str_col(df[col])

    # Coerce value to numeric
    df[VALUE_COL] = pd.to_numeric(df[VALUE_COL], errors="coerce")

    # Choose aggregation
    agg = agg.lower().strip()
    if agg not in {"mean", "median"}:
        raise ValueError("agg must be 'mean' or 'median'.")
    agg_fn = np.nanmean if agg == "mean" else np.nanmedian

    # Ensure all A-keys exist to make groupby deterministic
    for col in A_KEYS:
        if col not in df.columns:
            df[col] = ""

    merged_rows = []

    # Group by A-keys (sort=False to preserve group encounter order)
    for _, sub in df.groupby(A_KEYS, sort=False, dropna=False):
        # Build connected components via ANY B-key match
        components = _build_components_with_any_b_match(sub)

        for comp_global_indices in components:
            comp_df = df.loc[comp_global_indices]

            # Representative = earliest original row in this component
            # (we compare by AUX_IDX value, but idxmin returns the label)
            rep_label = comp_df[AUX_IDX].idxmin()
            rep_row = df.loc[rep_label].copy()

            # Aggregate the 'value' over the component (ignoring NaN)
            agg_val = agg_fn(comp_df[VALUE_COL].values)
            if not np.isnan(agg_val):
                rep_row[VALUE_COL] = float(agg_val)
            # If all NaN, keep representative's existing value (which could be NaN)

            # Keep representative row with AUX_IDX for sorting; drop later
            merged_rows.append(rep_row)

    # Build the final DataFrame
    out_df = pd.DataFrame(merged_rows)

    # Guarantee we return a valid DataFrame with the original columns
    if out_df.empty:
        out_df = pd.DataFrame(columns=original_cols)
        return out_df.reset_index(drop=True)

    # Sort by original order of representatives, then drop helper column
    if AUX_IDX in out_df.columns:
        try:
            out_df = out_df.sort_values(by=AUX_IDX, kind="stable")
        except TypeError:
            # Older pandas may not support kind="stable"
            out_df = out_df.sort_values(by=AUX_IDX)
        out_df = out_df.drop(columns=[AUX_IDX])

    # Restore original column order
    out_df = out_df[original_cols].reset_index(drop=True)
    return out_df


def merge_csv(input_csv: str, output_csv: str, agg: str = "mean") -> pd.DataFrame:
    """
    Convenience wrapper: read CSV → merge → write CSV with Windows-safe replace.

    Args:
        input_csv  : input CSV path
        output_csv : output CSV path
        agg        : 'mean' or 'median'

    Returns:
        The merged DataFrame (also written to disk).
    """
    # Read as strings to avoid unintended type inference; we'll coerce 'value' later.
    df = pd.read_csv(input_csv, dtype=str)

    merged = merge_dataframe(df, agg=agg)

    # Ensure output directory exists
    out_path = Path(output_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Write via temporary file, then atomic replace (best effort on Windows)
    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
    try:
        merged.to_csv(tmp_path, index=False, encoding="utf-8")
        os.replace(tmp_path, out_path)  # may raise PermissionError if target is open in Excel
        print(f"[ok] saved to {out_path}")
    except PermissionError:
        # Fall back to a new name if the target is locked
        fallback = out_path.with_name(out_path.stem + "_new" + out_path.suffix)
        merged.to_csv(fallback, index=False, encoding="utf-8")
        print(f"[warn] PermissionError writing '{out_path}'. "
              f"Is it open in Excel?\n→ Wrote to '{fallback}' instead.")
    finally:
        with contextlib.suppress(FileNotFoundError):
            os.remove(tmp_path)

    return merged


# ---------------- CLI ----------------
def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Merge CSV rows where A-keys match and ANY B-key matches.")
    p.add_argument("--input", "-i", required=True, help="Input CSV path")
    p.add_argument("--output", "-o", required=True, help="Output CSV path")
    p.add_argument("--agg", "-a", choices=["mean", "median"], default="mean",
                   help="Aggregation method for 'value' (default: mean)")
    return p


if __name__ == "__main__":
    args = _build_argparser().parse_args()
    merge_csv(args.input, args.output, agg=args.agg)
