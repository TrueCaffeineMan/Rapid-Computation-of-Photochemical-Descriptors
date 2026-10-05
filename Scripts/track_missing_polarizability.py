#!/usr/bin/env python3
"""
Track missing conformer polarizabilities in stats.csv files from properties.py.

Utility
-------
properties.py stores polarizability_au for each conformer when the value can be
parsed from the xTB output. In the current ensemble averaging step, a missing
polarizability value is treated as 0.0. This is useful for avoiding a crash,
but it can bias the ensemble polarizability downward if one or more conformers
with missing values have non-negligible Boltzmann weights.

This script does not rerun xTB and does not edit stats.csv. It scans either the
current molecule directory or the immediate subdirectories of the current
working directory, reads stats.csv, and writes one CSV file per molecule:

    polarizability_missing_tracker.csv

Each per-molecule file contains conformer-level rows plus a final __SUMMARY__
row. The script also writes a pooled summary in the starting directory:

    polarizability_missing_pooled.csv

The summary reports the Boltzmann weight coverage of conformers with valid
polarizability values, the weight assigned to missing values, the zero-filled
weighted value, and the renormalized present-only average.

Prerequisites
-------------
Python 3 standard library only.

Expected directory layout
-------------------------
Run either:

1. inside one molecule directory containing stats.csv

or

2. in a parent directory containing immediate molecule subdirectories, each
   containing stats.csv.
"""

import csv
import math
import os
from typing import Dict, List, Optional, Tuple


# =========================
# Configuration
# =========================

STATS_CSV_NAME = "stats.csv"
OUTPUT_CSV_NAME = "polarizability_missing_tracker.csv"
POOLED_OUTPUT_CSV_NAME = "polarizability_missing_pooled.csv"


# =========================
# Helpers
# =========================

def parse_float(value: object) -> Optional[float]:
    """Parse a finite float from a CSV field; return None for blank/invalid values."""
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    if not math.isfinite(number):
        return None
    return number


def format_float(value: Optional[float], digits: int = 10) -> str:
    """Format floats compactly for CSV output; keep missing values blank."""
    if value is None:
        return ""
    return f"{value:.{digits}g}"


def load_stats(path: str) -> List[Dict[str, str]]:
    """Load stats.csv rows as dictionaries."""
    with open(path, "r", newline="") as handle:
        return list(csv.DictReader(handle))


def discover_molecule_dirs(start_dir: str) -> List[str]:
    """
    Find molecule calculation directories.

    If the current directory itself contains stats.csv, audit only it.
    Otherwise, audit immediate subdirectories that contain stats.csv.
    """
    if os.path.isfile(os.path.join(start_dir, STATS_CSV_NAME)):
        return [start_dir]

    molecule_dirs: List[str] = []
    for name in sorted(os.listdir(start_dir)):
        path = os.path.join(start_dir, name)
        if not os.path.isdir(path):
            continue
        if os.path.isfile(os.path.join(path, STATS_CSV_NAME)):
            molecule_dirs.append(path)
    return molecule_dirs


def classify_polarizability(value: Optional[float]) -> str:
    """Return a compact status label for polarizability_au."""
    if value is None:
        return "missing"
    return "present"


def summarize_one_molecule(molecule_dir: str) -> Tuple[List[Dict[str, str]], Dict[str, str]]:
    """Create conformer-level rows and one summary row for a molecule directory."""
    molecule_name = os.path.basename(os.path.abspath(molecule_dir))
    stats_path = os.path.join(molecule_dir, STATS_CSV_NAME)
    rows = load_stats(stats_path)

    conformer_rows = [
        row for row in rows
        if row.get("conf_id", "") and row.get("conf_id", "") != "ensemble"
    ]
    ensemble_rows = [row for row in rows if row.get("conf_id", "") == "ensemble"]
    ensemble = ensemble_rows[0] if ensemble_rows else {}

    detail_rows: List[Dict[str, str]] = []

    n_present = 0
    n_missing = 0
    weight_present = 0.0
    weight_missing = 0.0
    zero_filled_weighted_sum = 0.0
    present_weighted_sum = 0.0
    missing_conf_ids: List[str] = []

    for row in conformer_rows:
        conf_id = row.get("conf_id", "")
        polarizability = parse_float(row.get("polarizability_au"))
        weight = parse_float(row.get("boltz_weight"))
        if weight is None:
            weight = 0.0

        status = classify_polarizability(polarizability)
        if polarizability is None:
            n_missing += 1
            weight_missing += weight
            missing_conf_ids.append(conf_id)
            weighted_contribution = 0.0
        else:
            n_present += 1
            weight_present += weight
            weighted_contribution = polarizability * weight
            present_weighted_sum += weighted_contribution
            zero_filled_weighted_sum += weighted_contribution

        detail_rows.append(
            {
                "row_type": "conformer",
                "molecule": molecule_name,
                "molecule_path": os.path.abspath(molecule_dir),
                "conf_id": conf_id,
                "polarizability_status": status,
                "polarizability_au": format_float(polarizability),
                "boltz_weight": format_float(weight),
                "weighted_contribution_zero_filled": format_float(weighted_contribution),
                "E_xtb_Eh": row.get("E_xtb_Eh", ""),
                "relE_kJmol": row.get("relE_kJmol", ""),
                "ensemble_reported_polarizability_au": "",
                "n_conformers": "",
                "n_present": "",
                "n_missing": "",
                "weight_present": "",
                "weight_missing": "",
                "zero_filled_weighted_polarizability_au": "",
                "renormalized_present_only_polarizability_au": "",
                "bias_zero_filled_minus_renormalized_au": "",
                "relative_bias_percent": "",
                "missing_conf_ids": "",
            }
        )

    if weight_present > 0.0:
        renormalized_present_only = present_weighted_sum / weight_present
        bias = zero_filled_weighted_sum - renormalized_present_only
        relative_bias_percent = (
            100.0 * bias / renormalized_present_only
            if renormalized_present_only != 0.0
            else None
        )
    else:
        renormalized_present_only = None
        bias = None
        relative_bias_percent = None

    ensemble_reported = parse_float(ensemble.get("polarizability_au"))

    summary = {
        "row_type": "__SUMMARY__",
        "molecule": molecule_name,
        "molecule_path": os.path.abspath(molecule_dir),
        "conf_id": "__SUMMARY__",
        "polarizability_status": "missing_detected" if n_missing else "complete",
        "polarizability_au": "",
        "boltz_weight": "",
        "weighted_contribution_zero_filled": "",
        "E_xtb_Eh": "",
        "relE_kJmol": "",
        "ensemble_reported_polarizability_au": format_float(ensemble_reported),
        "n_conformers": str(len(conformer_rows)),
        "n_present": str(n_present),
        "n_missing": str(n_missing),
        "weight_present": format_float(weight_present),
        "weight_missing": format_float(weight_missing),
        "zero_filled_weighted_polarizability_au": format_float(zero_filled_weighted_sum),
        "renormalized_present_only_polarizability_au": format_float(renormalized_present_only),
        "bias_zero_filled_minus_renormalized_au": format_float(bias),
        "relative_bias_percent": format_float(relative_bias_percent),
        "missing_conf_ids": ";".join(missing_conf_ids),
    }

    detail_rows.append(summary)
    return detail_rows, summary


def fieldnames() -> List[str]:
    """Stable column order for detail and summary CSV files."""
    return [
        "row_type",
        "molecule",
        "molecule_path",
        "conf_id",
        "polarizability_status",
        "polarizability_au",
        "boltz_weight",
        "weighted_contribution_zero_filled",
        "E_xtb_Eh",
        "relE_kJmol",
        "ensemble_reported_polarizability_au",
        "n_conformers",
        "n_present",
        "n_missing",
        "weight_present",
        "weight_missing",
        "zero_filled_weighted_polarizability_au",
        "renormalized_present_only_polarizability_au",
        "bias_zero_filled_minus_renormalized_au",
        "relative_bias_percent",
        "missing_conf_ids",
    ]


def write_rows(path: str, rows: List[Dict[str, str]]) -> None:
    """Write rows with a stable header."""
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames(), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


# =========================
# Main logic
# =========================

def main() -> None:
    """Run the missing-polarizability audit from the current working directory."""
    start_dir = os.getcwd()
    molecule_dirs = discover_molecule_dirs(start_dir)

    if not molecule_dirs:
        print("No stats.csv files found. Nothing written.")
        return

    pooled_summary_rows: List[Dict[str, str]] = []
    for molecule_dir in molecule_dirs:
        detail_rows, summary = summarize_one_molecule(molecule_dir)
        out_path = os.path.join(molecule_dir, OUTPUT_CSV_NAME)
        write_rows(out_path, detail_rows)
        pooled_summary_rows.append(summary)

        molecule_name = os.path.basename(os.path.abspath(molecule_dir))
        print(
            f"{molecule_name}: "
            f"missing {summary['n_missing']} of {summary['n_conformers']} conformer polarizabilities; "
            f"wrote {OUTPUT_CSV_NAME}"
        )

    pooled_path = os.path.join(start_dir, POOLED_OUTPUT_CSV_NAME)
    write_rows(pooled_path, pooled_summary_rows)
    print(f"Wrote pooled summary to {pooled_path}")


if __name__ == "__main__":
    main()