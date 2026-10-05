#!/usr/bin/env python3
"""
Compute conformer-level and ensemble-averaged photochemical descriptors.

Utility
-------
This is the main workflow script used after CREST conformational sampling.
It expects to be run inside one molecule directory containing:

    crest_conformers.xyz

The script performs the following operations:

1. Splits crest_conformers.xyz into individual conformer directories:
       conformers/conf_001/conf.xyz
       conformers/conf_002/conf.xyz
       ...

2. For each conformer, runs GFN2-xTB geometry optimization with ALPB(DMSO),
   tight optimization criteria, and Molden output generation.
IMPORTANT: Molden files take much disc volume!

3. Parses:
       - total xTB energy
       - molecular polarizability, when available

4. Runs xtb4stda on the optimized xtb geometry.

5. Runs sTDA-xTB and parses vertical excitation energies, excitation
   wavelengths, and oscillator strengths.

6. Selects the lowest-energy bright state using BRIGHT_F_THRESHOLD. If no state
   reaches the brightness threshold, the brightest state is selected as a
   fallback.

7. Computes conformer-level descriptors:
       E_S1_eV
       lambda_S1_nm
       f_S1
       dE_S1_405_eV
       dE_S1_420_eV
       dE_S1_450_eV
       f_sum_window

8. Computes Boltzmann weights from optimized GFN2-xTB energies at 298.15 K and
   writes ensemble-averaged descriptors.

9. For the lowest-energy conformer, runs an ORCA/otherm thermochemistry step
   and estimates a solution phase Gibbs energy by adding an xTB ALPB(DMSO)
   single-point solvation correction to the gas-phase QRRHO Gibbs energy.

Output
------
stats.csv

The output contains one row per conformer and a final row:

    conf_id == "ensemble"

Prerequisites
-------------
External programs required:

- xtb
- xtb4stda
- sTDA-xTB / stda
- ORCA
- otherm (by the Duarte group) https://github.com/duartegroup/otherm

Python standard library modules are used only; no additional Python packages
are required by this script.

BEFORE USE
----------
Edit the external executable paths in the "External program paths" section
below.
"""

import os
import csv
import math
import subprocess
import re
import shutil


# =========================
# External program paths
# =========================
# Replace these dummy paths with the correct paths on the system where the
# workflow is run. These are the only machine-specific paths in the script.

XTB_PATH = "/path/to/xtb"
XTB4STDA_PATH = "/path/to/xtb4stda"
STDA_PATH = "/path/to/stda"
ORCA_PATH = "/path/to/orca"
OTHERM_PATH = "/path/to/otherm.py"


# =========================
# Configuration
# =========================

CREST_FILE = "crest_conformers.xyz"
CONFORMER_ROOT = "conformers"

# Excitation energy limit for sTDA in eV, tailored for the present study.
STDA_E_MAX = 4.5

# Lamp wavelengths in nm.
LAMP_WAVELENGTHS = [405.0, 420.0, 450.0]

# Energy window for f_sum_window in eV.
WINDOW_MIN = 2.7
WINDOW_MAX = 3.5

# Bright-state threshold for oscillator strength.
BRIGHT_F_THRESHOLD = 0.05

# Boltzmann parameters.
TEMPERATURE = 298.15  # K
R_KJ_MOL_K = 8.314462618e-3  # kJ mol-1 K-1
EH_TO_KJ_MOL = 2625.49962  # 1 Eh in kJ mol-1


# =========================
# Helpers
# =========================

def split_crest_xyz(filename, out_root):
    """
    Split a CREST multi-structure XYZ file into individual conformer folders.

    Parameters
    ----------
    filename : str
        Name/path of the CREST multi-structure XYZ file.
    out_root : str
        Output directory for conformer folders.

    Returns
    -------
    list of tuple
        Each tuple contains:
            conf_id, conf_dir, conf_xyz_path, crest_energy_au_or_None
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"{filename} not found")

    os.makedirs(out_root, exist_ok=True)
    conformers = []

    with open(filename, "r") as f:
        idx = 0
        while True:
            line = f.readline()
            if not line:
                break
            line = line.strip()
            if not line:
                continue
            try:
                natoms = int(line.split()[0])
            except ValueError:
                # Corrupted line or footer.
                break

            comment = f.readline()
            if not comment:
                break
            comment = comment.rstrip("\n")

            coords = []
            for _ in range(natoms):
                coord_line = f.readline()
                if not coord_line:
                    raise IOError("Unexpected EOF while reading coordinates")
                coords.append(coord_line.rstrip("\n"))

            idx += 1
            conf_id = f"conf_{idx:03d}"
            conf_dir = os.path.join(out_root, conf_id)
            os.makedirs(conf_dir, exist_ok=True)

            conf_xyz = os.path.join(conf_dir, "conf.xyz")
            with open(conf_xyz, "w") as out:
                out.write(f"{natoms}\n")
                out.write(comment + "\n")
                for c in coords:
                    out.write(c + "\n")

            # Try to parse energy from the CREST comment line.
            crest_energy = None
            tokens = comment.strip().split()
            for t in tokens:
                try:
                    crest_energy = float(t)
                    break
                except ValueError:
                    continue

            conformers.append((conf_id, conf_dir, conf_xyz, crest_energy))

    return conformers


def run_xtb(conf_dir, conf_xyz):
    """
    Run GFN2-xTB optimization with ALPB(DMSO) and Molden output.

    Returns
    -------
    tuple
        total energy in Eh, molecular polarizability in a.u. or None.
    """
    cmd = [
        XTB_PATH,
        os.path.basename(conf_xyz),
        "--gfn", "2",
        "--alpb", "DMSO",
        "--opt", "tight",
        "--molden",
    ]
    res = subprocess.run(
        cmd,
        cwd=conf_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    xtb_out_path = os.path.join(conf_dir, "xtb.out")
    with open(xtb_out_path, "w") as f:
        f.write(res.stdout)

    if res.returncode != 0:
        raise RuntimeError(f"xtb failed in {conf_dir}:\n{res.stdout}")

    text = res.stdout

    m = re.search(r"total\s+energy\s+(-?\d+\.\d+)", text, re.IGNORECASE)
    if not m:
        if os.path.exists(xtb_out_path):
            with open(xtb_out_path, "r") as f:
                text2 = f.read()
            m = re.search(r"total\s+energy\s+(-?\d+\.\d+)", text2, re.IGNORECASE)
    if not m:
        raise RuntimeError(f"Could not find total energy in xtb output for {conf_dir}")
    energy_au = float(m.group(1))

    # Parse molecular polarizability line:
    # "Mol. α(0) /au : ..."
    # A dot wildcard is used for the Greek alpha to avoid encoding issues.
    m_pol = re.search(r"Mol\.\s*.\(0\)\s*/au\s*:\s*([+-]?\d+\.\d+(?:[Ee][+-]?\d+)?)", text)
    if not m_pol and os.path.exists(xtb_out_path):
        with open(xtb_out_path, "r") as f:
            text2 = f.read()
        m_pol = re.search(r"Mol\.\s*.\(0\)\s*/au\s*:\s*([+-]?\d+\.\d+(?:[Ee][+-]?\d+)?)", text2)

    if m_pol:
        polarizability_au = float(m_pol.group(1))
    else:
        # Missing values are allowed and tracked separately by
        # track_missing_polarizability.py.
        polarizability_au = None

    return energy_au, polarizability_au


def run_xtb4stda(conf_dir, geom_file="xtbopt.xyz"):
    """Run xtb4stda on the optimized geometry."""
    geom_path = os.path.join(conf_dir, geom_file)
    if not os.path.exists(geom_path):
        raise FileNotFoundError(f"{geom_path} not found for xtb4stda")

    cmd = [XTB4STDA_PATH, os.path.basename(geom_path)]
    res = subprocess.run(
        cmd,
        cwd=conf_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if res.returncode != 0:
        raise RuntimeError(f"xtb4stda failed in {conf_dir}:\n{res.stdout}")


def run_stda(conf_dir, out_name="sTDA.out"):
    """Run sTDA and capture output to out_name."""
    out_path = os.path.join(conf_dir, out_name)
    with open(out_path, "w") as fout:
        res = subprocess.run(
            [STDA_PATH, "-xtb", "-e", str(STDA_E_MAX)],
            cwd=conf_dir,
            stdout=fout,
            stderr=subprocess.STDOUT,
            text=True,
        )
    if res.returncode != 0:
        raise RuntimeError(f"stda failed in {conf_dir}, see {out_name}")
    return out_path


def parse_stda_output(stda_file):
    """
    Parse excited states from sTDA output.

    Returns
    -------
    list of dict
        Each dictionary contains:
            index, E, nm, f
    """
    states = []
    with open(stda_file, "r") as f:
        lines = f.readlines()

    start = None
    for i, line in enumerate(lines):
        if "excitation energies, transition moments and TDA amplitudes" in line:
            start = i
            break

    if start is None:
        raise RuntimeError(f"No excitation block found in {stda_file}")

    for line in lines[start + 2:]:
        if line.strip() == "" or line.lstrip().startswith("alpha tensor"):
            break
        stripped = line.strip()
        if not stripped:
            break
        if not stripped[0].isdigit():
            continue
        tokens = line.split()
        if len(tokens) < 4:
            continue
        try:
            idx = int(tokens[0])
            E_eV = float(tokens[1])
            lam_nm = float(tokens[2])
            fL = float(tokens[3])
        except ValueError:
            continue
        states.append(
            {"index": idx, "E": E_eV, "nm": lam_nm, "f": fL}
        )

    if not states:
        raise RuntimeError(f"No states parsed from {stda_file}")

    return states


def pick_bright_lowest_state(states, f_threshold=BRIGHT_F_THRESHOLD):
    """
    Pick the lowest-energy bright state.

    The first selected state is the lowest-energy state with f >= f_threshold.
    If no state reaches the threshold, the state with the largest oscillator
    strength is selected.
    """
    bright = [s for s in states if s["f"] >= f_threshold]
    if bright:
        return min(bright, key=lambda s: s["E"])
    return max(states, key=lambda s: s["f"])


def compute_descriptors(states):
    """
    Compute conformer-level excitation descriptors from parsed sTDA states.
    """
    chosen = pick_bright_lowest_state(states)
    desc = {
        "E_S1_eV": chosen["E"],
        "lambda_S1_nm": chosen["nm"],
        "f_S1": chosen["f"],
    }

    for lam in LAMP_WAVELENGTHS:
        E_exp = 1240.0 / lam
        key = f"dE_S1_{int(lam)}_eV"
        desc[key] = chosen["E"] - E_exp

    f_sum = sum(
        s["f"]
        for s in states
        if WINDOW_MIN <= s["E"] <= WINDOW_MAX
    )
    desc["f_sum_window"] = f_sum

    return desc


def boltzmann_weights(energies_au):
    """Return Boltzmann weights from energies in Hartree."""
    E_min = min(energies_au)
    weights = []
    for E in energies_au:
        dE_kJ = (E - E_min) * EH_TO_KJ_MOL
        w = math.exp(-dE_kJ / (R_KJ_MOL_K * TEMPERATURE))
        weights.append(w)
    norm = sum(weights)
    return [w / norm for w in weights], E_min


def run_xtb_sp(workdir, xyz_name, use_alpb):
    """
    Run a single-point xTB calculation on xyz_name in workdir.

    Parameters
    ----------
    use_alpb : bool
        If True, use ALPB(DMSO). If False, run a gas-phase single point.

    Returns
    -------
    float
        Total energy in Eh.
    """
    cmd = [
        XTB_PATH,
        xyz_name,
        "--gfn", "2",
    ]
    if use_alpb:
        cmd += ["--alpb", "DMSO"]

    res = subprocess.run(
        cmd,
        cwd=workdir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if res.returncode != 0:
        raise RuntimeError(f"xTB SP failed in {workdir} (ALPB={use_alpb}):\n{res.stdout}")

    text = res.stdout
    m = re.search(r"total\s+energy\s+(-?\d+\.\d+)", text, re.IGNORECASE)
    if not m:
        raise RuntimeError(f"Could not find total energy in xTB SP output for {workdir}")
    return float(m.group(1))


# =========================
# Main logic
# =========================

def main():
    """Run the complete conformer-level and ensemble descriptor workflow."""
    base_dir = os.getcwd()
    crest_path = os.path.join(base_dir, CREST_FILE)
    if not os.path.exists(crest_path):
        raise FileNotFoundError(
            f"{CREST_FILE} not found in {base_dir}. Run this script where crest_conformers.xyz is."
        )

    print(f"Splitting {CREST_FILE} into individual conformers...")
    conformers = split_crest_xyz(crest_path, CONFORMER_ROOT)
    print(f"Found {len(conformers)} conformers.")

    conf_data = []
    energies_au = []

    lowest_E_xtb = None
    lowest_conf_dir_xtb = None

    for conf_id, conf_dir, conf_xyz, crest_E in conformers:
        print(f"Processing {conf_id} ...")

        # 1) xTB optimization:
        #    writes xtbopt.xyz, molden.input, xtb.out; returns energy and polarizability.
        E_xtb, polarizability_au = run_xtb(conf_dir, conf_xyz)

        if lowest_E_xtb is None or E_xtb < lowest_E_xtb:
            lowest_E_xtb = E_xtb
            lowest_conf_dir_xtb = conf_dir

        # 2) Prepare sTDA-xTB input on optimized geometry.
        run_xtb4stda(conf_dir, "xtbopt.xyz")

        # 3) Run sTDA-xTB.
        stda_out = run_stda(conf_dir, "sTDA.out")

        # 4) Parse sTDA output. If no excitation block is present, skip this
        #    conformer from ensemble averaging.
        try:
            states = parse_stda_output(stda_out)
        except RuntimeError as e:
            print(f"[WARN] {conf_id}: {e} -> skipping this conformer from ensemble.")
            continue

        desc = compute_descriptors(states)

        rec = {
            "conf_id": conf_id,
            "conf_dir": os.path.relpath(conf_dir, base_dir),
            "E_xtb_Eh": E_xtb,
            "E_crest_Eh": crest_E if crest_E is not None else "",
            "polarizability_au": polarizability_au,
        }
        rec.update(desc)
        conf_data.append(rec)
        energies_au.append(E_xtb)

    if not conf_data:
        print("No conformers with excitation data; nothing to average.")
        return

    # 5) Boltzmann weights from xTB energies.
    print("Computing Boltzmann weights at T = {:.2f} K".format(TEMPERATURE))
    w_list, E_min = boltzmann_weights(energies_au)

    for rec, w, E in zip(conf_data, w_list, energies_au):
        dE_kJ = (E - E_min) * EH_TO_KJ_MOL
        rec["relE_kJmol"] = dE_kJ
        rec["boltz_weight"] = w

    # 6) ORCA + otherm QRRHO on lowest-energy conformer.
    print("Running ORCA + otherm QRRHO on lowest-energy conformer...")
    if lowest_conf_dir_xtb is None:
        raise RuntimeError("No lowest-energy conformer recorded for ORCA/otherm step.")

    orca_dir = os.path.join(base_dir, "orca")
    os.makedirs(orca_dir, exist_ok=True)

    xtbopt_path = os.path.join(lowest_conf_dir_xtb, "xtbopt.xyz")
    if not os.path.exists(xtbopt_path):
        raise FileNotFoundError(f"{xtbopt_path} not found for ORCA input")
    init_xyz_path = os.path.join(orca_dir, "init.xyz")
    shutil.copyfile(xtbopt_path, init_xyz_path)

    opt_inp_path = os.path.join(orca_dir, "opt.inp")
    with open(opt_inp_path, "w") as f:
        f.write("""# orca_run
! TightOpt NUMFREQ XTB2

* xyzfile 0 1 init.xyz

%xtb
ACCURACY 0.01
end

%scf
MaxIter 800
end

%pal
nprocs 1
end

%geom
MaxIter 500
end
""")

    opt_out_path = os.path.join(orca_dir, "opt.out")
    with open(opt_out_path, "w") as f_out:
        res_orca = subprocess.run(
            [ORCA_PATH, "opt.inp"],
            cwd=orca_dir,
            stdout=f_out,
            stderr=subprocess.STDOUT,
            text=True,
        )
    if res_orca.returncode != 0:
        raise RuntimeError(f"ORCA failed in {orca_dir}, see opt.out")

    thermo_out_path = os.path.join(orca_dir, "thermo.out")
    res_otherm = subprocess.run(
        [
            OTHERM_PATH,
            "-t", f"{TEMPERATURE:.2f}",
            "-ss", "1M",
            "-m", "grimme",
            "-sn", "1",
            "-r", "opt.out",
        ],
        cwd=orca_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    with open(thermo_out_path, "w") as ft:
        ft.write(res_otherm.stdout)
    if res_otherm.returncode != 0:
        raise RuntimeError(f"otherm failed in {orca_dir}, see thermo.out")

    # Parse G from "For convenience E, H, G in Hartrees" block.
    G_QRRHO_Eh = None
    lines = res_otherm.stdout.splitlines()
    for i, line in enumerate(lines):
        if "For convenience E, H, G in Hartrees" in line:
            if i + 1 < len(lines):
                nums_line = lines[i + 1].strip()
                parts = [p.strip() for p in nums_line.split(",")]
                if len(parts) >= 3:
                    try:
                        G_QRRHO_Eh = float(parts[2])
                    except ValueError:
                        G_QRRHO_Eh = None
            break
    if G_QRRHO_Eh is None:
        raise RuntimeError(f"Could not parse G from {thermo_out_path}")

    # Estimate xTB ALPB(DMSO) solvation correction:
    # G_solution ≈ G_QRRHO_gas + [E_xTB_ALPB(DMSO) - E_xTB_gas]
    opt_xyz_for_xtb = "opt.xyz"
    opt_xyz_path = os.path.join(orca_dir, opt_xyz_for_xtb)
    if not os.path.exists(opt_xyz_path):
        raise FileNotFoundError(f"{opt_xyz_path} not found for xTB SP solvation step")

    Egas_Eh = run_xtb_sp(orca_dir, opt_xyz_for_xtb, use_alpb=False)
    Esol_Eh = run_xtb_sp(orca_dir, opt_xyz_for_xtb, use_alpb=True)

    G_QRRHO_solv_Eh = G_QRRHO_Eh + (Esol_Eh - Egas_Eh)

    # 7) Ensemble-averaged descriptors.
    ensemble = {"conf_id": "ensemble", "conf_dir": ""}

    E_ens = sum(r["E_xtb_Eh"] * r["boltz_weight"] for r in conf_data)
    ensemble["E_xtb_Eh"] = E_ens
    ensemble["E_crest_Eh"] = ""
    ensemble["relE_kJmol"] = 0.0
    ensemble["boltz_weight"] = 1.0

    ensemble["G_QRRHO_Eh"] = G_QRRHO_Eh
    ensemble["G_QRRHO_solv_Eh"] = G_QRRHO_solv_Eh

    keys_to_avg = [
        "E_S1_eV",
        "f_S1",
        "f_sum_window",
        "polarizability_au",
    ]
    for lam in LAMP_WAVELENGTHS:
        keys_to_avg.append(f"dE_S1_{int(lam)}_eV")

    for key in keys_to_avg:
        val = sum(
            (r[key] if r[key] is not None else 0.0) * r["boltz_weight"]
            for r in conf_data
        )
        ensemble[key] = val

    if ensemble.get("E_S1_eV", None):
        ensemble["lambda_S1_nm"] = 1240.0 / ensemble["E_S1_eV"]
    else:
        ensemble["lambda_S1_nm"] = ""

    conf_data.append(ensemble)

    # 8) Write CSV.
    out_csv = os.path.join(base_dir, "stats.csv")
    print(f"Writing {out_csv}")

    fieldnames = sorted({k for rec in conf_data for k in rec.keys()})

    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for rec in conf_data:
            writer.writerow(rec)

    print("Done.")


if __name__ == "__main__":
    main()