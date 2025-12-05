#!/usr/bin/env python3
"""
fill_pdbpath.py
Populate/replace a 'pdbpath' column in a CSV based on first occurrence order of 'ProteinID'.

Rules:
- The first unique ProteinID seen → "seq1.pdb", the next new ProteinID → "seq2.pdb", etc.
- Repeated ProteinID rows reuse the same assigned seqX.pdb as their first occurrence.
- Rows with missing/blank ProteinID get an empty pdbpath.
- If a 'pdbpath' column already exists, it will be overwritten unless --only-fill-empty is used.

Usage:
  python fill_pdbpath.py input.csv               # writes input_with_pdbpath.csv
  python fill_pdbpath.py input.csv -o output.csv # writes to output.csv
  python fill_pdbpath.py input.csv --in-place    # edits the file in place
  python fill_pdbpath.py input.csv --only-fill-empty  # keep non-empty pdbpath cells intact
"""

import argparse
import os
import sys
import pandas as pd
import numpy as np


def fill_pdbpath(
    df: pd.DataFrame,
    id_col: str = "ProteinID",
    pdb_col: str = "pdbpath",
    prefix: str = "seq",
    ext: str = ".pdb",
    only_fill_empty: bool = False,
) -> pd.DataFrame:
    if id_col not in df.columns:
        raise KeyError(f"Column '{id_col}' not found in CSV.")

    # Ensure pdb_col exists
    if pdb_col not in df.columns:
        df[pdb_col] = ""

    # Normalize ProteinID to pandas' nullable string for robust NA handling
    s = df[id_col].astype("string").str.strip()

    # Factorize in first-occurrence order; NA → code -1
    codes, uniques = pd.factorize(s, sort=False)

    # codes >= 0 → assign seq{code+1}.pdb; else empty
    assigned_idx = pd.Series(codes + 1, index=df.index)
    assigned_idx = assigned_idx.where(codes >= 0, pd.NA)

    # Build the target series
    target = assigned_idx.map(lambda x: f"{prefix}{int(x)}{ext}" if pd.notna(x) else "")

    if only_fill_empty:
        # Only write where the existing pdb_col is empty/NA
        is_empty = df[pdb_col].astype("string").fillna("").str.strip().eq("")
        df.loc[is_empty, pdb_col] = target[is_empty]
    else:
        df[pdb_col] = target

    return df


def main():
    ap = argparse.ArgumentParser(description="Fill 'pdbpath' based on first occurrence of ProteinID.")
    ap.add_argument("input_csv", help="Path to input CSV")
    ap.add_argument("-o", "--output", dest="output_csv", help="Path to output CSV (default: <input>_with_pdbpath.csv)")
    ap.add_argument("--in-place", action="store_true", help="Edit the input file in place")
    ap.add_argument("--id-col", default="ProteinID", help="Name of ProteinID column (default: ProteinID)")
    ap.add_argument("--pdb-col", default="pdbpath", help="Name of pdbpath column (default: pdbpath)")
    ap.add_argument("--prefix", default="seq", help="Prefix for generated file names (default: seq)")
    ap.add_argument("--ext", default=".pdb", help="Extension for generated file names (default: .pdb)")
    ap.add_argument("--only-fill-empty", action="store_true", help="Only fill empty pdbpath cells; keep existing non-empty values")
    ap.add_argument("--encoding", default=None, help="Optional file encoding for reading/writing (e.g., utf-8, gbk)")
    args = ap.parse_args()

    if args.in_place and args.output_csv:
        print("Error: --in-place and --output cannot be used together.", file=sys.stderr)
        sys.exit(2)

    # Decide output path
    if args.in_place:
        output_path = args.input_csv  # overwrite
    else:
        if args.output_csv:
            output_path = args.output_csv
        else:
            root, ext = os.path.splitext(args.input_csv)
            output_path = f"{root}_with_pdbpath.csv"

    # Read
    try:
        df = pd.read_csv(args.input_csv, dtype=str, encoding=args.encoding)
    except Exception as e:
        print(f"Failed to read CSV: {e}", file=sys.stderr)
        sys.exit(1)

    # Process
    try:
        df = fill_pdbpath(
            df,
            id_col=args.id_col,
            pdb_col=args.pdb_col,
            prefix=args.prefix,
            ext=args.ext,
            only_fill_empty=args.only_fill_empty,
        )
    except Exception as e:
        print(f"Error while filling pdbpath: {e}", file=sys.stderr)
        sys.exit(1)

    # Write
    try:
        df.to_csv(output_path, index=False, encoding=args.encoding if args.encoding else "utf-8")
    except Exception as e:
        print(f"Failed to write CSV: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Done. Wrote: {output_path}")


if __name__ == "__main__":
    main()
