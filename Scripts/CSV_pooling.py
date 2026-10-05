#!/usr/bin/env python3
"""
Pool molecule-level ensemble descriptors into one CSV table.

Utility
-------
properties.py writes one stats.csv file per molecule directory. Each stats.csv
contains conformer-level rows and one final row with:

    conf_id == "ensemble"

This script is intended to be run in a parent directory containing its
subdirectories for individual molecules. It scans those subdirectories, reads
each stats.csv file, extracts the ensemble row, adds the subdirectory name as
the molecule name, and writes:

    ensembles_pooled.csv

Prerequisites
-------------
Python 3 with pandas installed.

Expected directory layout
-------------------------
Run from a parent directory containing:

    molecule_1/stats.csv
    molecule_2/stats.csv
    molecule_3/stats.csv
    ...

Only immediate subdirectories are scanned.
"""

import os
import csv
import pandas as pd


# =========================
# Configuration
# =========================

INPUT_CSV_NAME = "stats.csv"
OUTPUT_CSV_NAME = "ensembles_pooled.csv"


# =========================
# Main logic
# =========================

def main():
    """Pool ensemble rows from all immediate subdirectories."""
    base_dir = os.getcwd()

    pooled_rows = []
    column_order = None  # taken from the first stats.csv that is processed

    # Loop over immediate subdirectories only.
    for folder_name in sorted(os.listdir(base_dir)):
        folder_path = os.path.join(base_dir, folder_name)
        if not os.path.isdir(folder_path):
            continue

        csv_path = os.path.join(folder_path, INPUT_CSV_NAME)
        if not os.path.isfile(csv_path):
            continue

        df = pd.read_csv(csv_path)

        if "conf_id" not in df.columns:
            print(f"[WARN] {csv_path}: no 'conf_id' column, skipping.")
            continue

        ensemble_rows = df[df["conf_id"] == "ensemble"]
        if ensemble_rows.empty:
            print(f"[WARN] {csv_path}: no ensemble row, skipping.")
            continue

        ensemble = ensemble_rows.iloc[0].to_dict()

        # Add folder name as first field. The folder name is used as the
        # molecule/calculation identifier in the pooled table.
        ensemble_with_name = {"name": folder_name}
        ensemble_with_name.update(ensemble)

        pooled_rows.append(ensemble_with_name)

        # Preserve the column order from the first processed file.
        if column_order is None:
            column_order = ["name"] + list(df.columns)

    if not pooled_rows:
        print("No ensemble rows found in any subfolder. Nothing written.")
        return

    out_path = os.path.join(base_dir, OUTPUT_CSV_NAME)
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=column_order)
        writer.writeheader()
        for row in pooled_rows:
            writer.writerow(row)

    print(f"Wrote {len(pooled_rows)} ensemble rows to {out_path}")


if __name__ == "__main__":
    main()