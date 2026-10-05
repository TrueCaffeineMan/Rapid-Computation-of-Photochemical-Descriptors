#!/usr/bin/env python3
"""
Check conformer-level sTDA-xTB outputs for lower-energy dark states.

Utility
-------
properties.py selects the lowest-energy bright state, using
f >= BRIGHT_F_THRESHOLD, and stores that state as E_S1_eV/lambda_S1_nm/f_S1.
This is a useful photochemical descriptor, but it is not necessarily the
exact S1 state. If one or more lower-energy states have oscillator strengths
below the same threshold, the conformer has dark states below the selected
bright state. Such conformers should be flagged before interpreting E_S1_eV
as a literal S1 energy.

NOTE: No dark S1 states were found throughout the study.

This script scans either the current molecule directory or the immediate
subdirectories of the current working directory, reads:

    conformers/conf_*/sTDA.out

and writes one CSV file per molecule directory:

    dark_state_flags.csv

A pooled file is also written in the starting directory:

    dark_state_flags_pooled.csv

By default, the per-molecule and pooled CSV files contain only conformers with
lower-energy dark states below the selected bright state. Set INCLUDE_UNFLAGGED
to True to write all parsed conformers with an explicit flag column.

Prerequisites
-------------
Python 3 standard library only.

Expected directory layout
-------------------------
Run either:

1. inside one molecule directory containing:
       conformers/conf_001/sTDA.out
       conformers/conf_002/sTDA.out
       ...
       stats.csv

or

2. in a parent directory containing immediate molecule subdirectories, each with
   the same conformers/ and stats.csv layout.
"""

import csv
import os
from typing import Dict, Iterable, List, Optional, Tuple


# =========================
# Configuration
# =========================

STATS_CSV_NAME = "stats.csv"
CONFORMER_ROOT = "conformers"
STDA_OUT_NAME = "sTDA.out"
OUTPUT_CSV_NAME = "dark_state_flags.csv"
POOLED_OUTPUT_CSV_NAME = "dark_state_flags_pooled.csv"

# The same oscillator-strength threshold used in properties.py for selecting
# the lowest-energy bright state.
BRIGHT_F_THRESHOLD = 0.05

# False: report only flagged conformers.
# True: report every conformer, with a flag column.
INCLUDE_UNFLAGGED = False


# =========================
# Helpers
# =========================

def format_float(value: Optional[float], digits: int = 8) -> str:
    """Format floats compactly for CSV output; keep missing values blank."""
    if value is None:
        return ""
    return f"{value:.{digits}g}"


def join_values(values: Iterable[object]) -> str:
    """Join list-like values into a semicolon-separated CSV-safe string."""
    return ";".join(str(v) for v in values)


def parse_stda_output(stda_file: str) -> List[Dict[str, float]]:
    """
    Parse sTDA excited states from the standard output block.

    Returns
    -------
    list of dict
        Each dictionary contains:
            index, E_eV, lambda_nm, f
    """
    with open(stda_file, "r", errors="replace") as handle:
        lines = handle.readlines()

    start = None
    for i, line in enumerate(lines):
        if "excitation energies, transition moments and TDA amplitudes" in line:
            start = i
            break

    if start is None:
        raise RuntimeError("no excitation-energy block found")

    states: List[Dict[str, float]] = []
    for line in lines[start + 2:]:
        stripped = line.strip()
        if not stripped:
            break
        if stripped.lower().startswith("alpha tensor"):
            break
        if not stripped[0].isdigit():
            continue

        tokens = line.split()
        if len(tokens) < 4:
            continue

        try:
            state_index = int(tokens[0])
            energy_ev = float(tokens[1])
            lambda_nm = float(tokens[2])
            oscillator_strength = float(tokens[3])
        except ValueError:
            continue

        states.append(
            {
                "index": state_index,
                "E_eV": energy_ev,
                "lambda_nm": lambda_nm,
                "f": oscillator_strength,
            }
        )

    if not states:
        raise RuntimeError("no states parsed from excitation-energy block")

    return states


def pick_selected_bright_state(
    states: List[Dict[str, float]],
    f_threshold: float,
) -> Tuple[Dict[str, float], str]:
    """
    Apply the same state-selection rule as properties.py.

    First choice:
        lowest-energy state with f >= threshold.

    Fallback:
        state with maximal f if no state reaches the threshold.
    """
    bright_states = [state for state in states if state["f"] >= f_threshold]
    if bright_states:
        return min(bright_states, key=lambda state: state["E_eV"]), "lowest_f_ge_threshold"
    return max(states, key=lambda state: state["f"]), "max_f_fallback"


def load_stats_rows(molecule_dir: str) -> Dict[str, Dict[str, str]]:
    """Load stats.csv rows keyed by conf_id, if the file exists."""
    stats_path = os.path.join(molecule_dir, STATS_CSV_NAME)
    if not os.path.isfile(stats_path):
        return {}

    with open(stats_path, "r", newline="") as handle:
        reader = csv.DictReader(handle)
        return {
            row.get("conf_id", ""): row
            for row in reader
            if row.get("conf_id", "")
        }


def discover_molecule_dirs(start_dir: str) -> List[str]:
    """
    Find molecule calculation directories.

    If the current directory itself contains conformers/, audit only it.
    Otherwise, audit immediate subdirectories that contain conformers/.
    """
    current_conformers = os.path.join(start_dir, CONFORMER_ROOT)
    if os.path.isdir(current_conformers):
        return [start_dir]

    molecule_dirs: List[str] = []
    for name in sorted(os.listdir(start_dir)):
        path = os.path.join(start_dir, name)
        if not os.path.isdir(path):
            continue
        if os.path.isdir(os.path.join(path, CONFORMER_ROOT)):
            molecule_dirs.append(path)
    return molecule_dirs


def conformer_dirs(molecule_dir: str) -> List[Tuple[str, str, str]]:
    """Return tuples of (conf_id, conf_dir, sTDA_out_path)."""
    root = os.path.join(molecule_dir, CONFORMER_ROOT)
    if not os.path.isdir(root):
        return []

    result: List[Tuple[str, str, str]] = []
    for conf_id in sorted(os.listdir(root)):
        conf_dir = os.path.join(root, conf_id)
        if not os.path.isdir(conf_dir):
            continue
        stda_path = os.path.join(conf_dir, STDA_OUT_NAME)
        if os.path.isfile(stda_path):
            result.append((conf_id, conf_dir, stda_path))
    return result


def audit_one_molecule(molecule_dir: str) -> List[Dict[str, str]]:
    """Audit one molecule directory and write dark_state_flags.csv there."""
    stats_by_conf = load_stats_rows(molecule_dir)
    molecule_name = os.path.basename(os.path.abspath(molecule_dir))
    rows: List[Dict[str, str]] = []

    for conf_id, conf_dir, stda_path in conformer_dirs(molecule_dir):
        stats_row = stats_by_conf.get(conf_id, {})

        try:
            states = parse_stda_output(stda_path)
            selected, selection_rule = pick_selected_bright_state(
                states,
                BRIGHT_F_THRESHOLD,
            )
        except RuntimeError as exc:
            rows.append(
                {
                    "molecule": molecule_name,
                    "molecule_path": os.path.abspath(molecule_dir),
                    "conf_id": conf_id,
                    "conf_dir": os.path.relpath(conf_dir, molecule_dir),
                    "parse_ok": "False",
                    "flag_lower_dark_below_selected": "",
                    "error": str(exc),
                    "n_states": "",
                    "true_s1_index": "",
                    "true_s1_E_eV": "",
                    "true_s1_f": "",
                    "selected_state_index": "",
                    "selected_state_E_eV": "",
                    "selected_state_nm": "",
                    "selected_state_f": "",
                    "selected_state_rule": "",
                    "n_lower_dark_states": "",
                    "lowest_dark_index": "",
                    "lowest_dark_E_eV": "",
                    "lowest_dark_nm": "",
                    "lowest_dark_f": "",
                    "dark_state_indices_below_selected": "",
                    "dark_state_E_eV_below_selected": "",
                    "dark_state_f_below_selected": "",
                    "gap_selected_minus_lowest_dark_eV": "",
                    "reported_E_S1_eV": stats_row.get("E_S1_eV", ""),
                    "E_xtb_Eh": stats_row.get("E_xtb_Eh", ""),
                    "relE_kJmol": stats_row.get("relE_kJmol", ""),
                    "boltz_weight": stats_row.get("boltz_weight", ""),
                }
            )
            continue

        true_s1 = min(states, key=lambda state: state["E_eV"])
        lower_dark_states = [
            state
            for state in states
            if state["E_eV"] < selected["E_eV"] and state["f"] < BRIGHT_F_THRESHOLD
        ]
        flag = bool(lower_dark_states)

        if not flag and not INCLUDE_UNFLAGGED:
            continue

        lowest_dark = min(lower_dark_states, key=lambda state: state["E_eV"]) if flag else None
        gap = selected["E_eV"] - lowest_dark["E_eV"] if lowest_dark is not None else None

        rows.append(
            {
                "molecule": molecule_name,
                "molecule_path": os.path.abspath(molecule_dir),
                "conf_id": conf_id,
                "conf_dir": os.path.relpath(conf_dir, molecule_dir),
                "parse_ok": "True",
                "flag_lower_dark_below_selected": str(flag),
                "error": "",
                "n_states": str(len(states)),
                "true_s1_index": str(int(true_s1["index"])),
                "true_s1_E_eV": format_float(true_s1["E_eV"]),
                "true_s1_f": format_float(true_s1["f"]),
                "selected_state_index": str(int(selected["index"])),
                "selected_state_E_eV": format_float(selected["E_eV"]),
                "selected_state_nm": format_float(selected["lambda_nm"]),
                "selected_state_f": format_float(selected["f"]),
                "selected_state_rule": selection_rule,
                "n_lower_dark_states": str(len(lower_dark_states)),
                "lowest_dark_index": str(int(lowest_dark["index"])) if lowest_dark else "",
                "lowest_dark_E_eV": format_float(lowest_dark["E_eV"]) if lowest_dark else "",
                "lowest_dark_nm": format_float(lowest_dark["lambda_nm"]) if lowest_dark else "",
                "lowest_dark_f": format_float(lowest_dark["f"]) if lowest_dark else "",
                "dark_state_indices_below_selected": join_values(
                    int(state["index"]) for state in lower_dark_states
                ),
                "dark_state_E_eV_below_selected": join_values(
                    format_float(state["E_eV"]) for state in lower_dark_states
                ),
                "dark_state_f_below_selected": join_values(
                    format_float(state["f"]) for state in lower_dark_states
                ),
                "gap_selected_minus_lowest_dark_eV": format_float(gap),
                "reported_E_S1_eV": stats_row.get("E_S1_eV", ""),
                "E_xtb_Eh": stats_row.get("E_xtb_Eh", ""),
                "relE_kJmol": stats_row.get("relE_kJmol", ""),
                "boltz_weight": stats_row.get("boltz_weight", ""),
            }
        )

    out_path = os.path.join(molecule_dir, OUTPUT_CSV_NAME)
    write_rows(out_path, rows)
    return rows


def fieldnames() -> List[str]:
    """Stable column order for output CSV files."""
    return [
        "molecule",
        "molecule_path",
        "conf_id",
        "conf_dir",
        "parse_ok",
        "flag_lower_dark_below_selected",
        "error",
        "n_states",
        "true_s1_index",
        "true_s1_E_eV",
        "true_s1_f",
        "selected_state_index",
        "selected_state_E_eV",
        "selected_state_nm",
        "selected_state_f",
        "selected_state_rule",
        "n_lower_dark_states",
        "lowest_dark_index",
        "lowest_dark_E_eV",
        "lowest_dark_nm",
        "lowest_dark_f",
        "dark_state_indices_below_selected",
        "dark_state_E_eV_below_selected",
        "dark_state_f_below_selected",
        "gap_selected_minus_lowest_dark_eV",
        "reported_E_S1_eV",
        "E_xtb_Eh",
        "relE_kJmol",
        "boltz_weight",
    ]


def write_rows(path: str, rows: List[Dict[str, str]]) -> None:
    """Write rows with a stable header, even if the row list is empty."""
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames(), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


# =========================
# Main logic
# =========================

def main() -> None:
    """Run the dark-state audit from the current working directory."""
    start_dir = os.getcwd()
    molecule_dirs = discover_molecule_dirs(start_dir)

    if not molecule_dirs:
        print("No molecule directories with conformers/ found. Nothing written.")
        return

    pooled_rows: List[Dict[str, str]] = []
    for molecule_dir in molecule_dirs:
        rows = audit_one_molecule(molecule_dir)
        pooled_rows.extend(rows)
        print(
            f"{os.path.basename(os.path.abspath(molecule_dir))}: "
            f"wrote {len(rows)} rows to {OUTPUT_CSV_NAME}"
        )

    pooled_path = os.path.join(start_dir, POOLED_OUTPUT_CSV_NAME)
    write_rows(pooled_path, pooled_rows)
    print(f"Wrote {len(pooled_rows)} pooled rows to {pooled_path}")


if __name__ == "__main__":
    main()