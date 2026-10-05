# A Workflow for Rapid Computation of Photochemical Descriptors

This repository is the version-of-record code deposit accompanying the RSC Advances article **“Substituent-Controlled Visible Light Photoisomerization and Cyclodimerization of Styryl Pyrimidinones.”** The deposited Python code comprises the main computational workflow and three utility scripts. The workflow was developed for rapid screening of absorption-related photochemical descriptors from isomeric SMILES input.

## Version of Record

- **Article:** *Substituent-Controlled Visible Light Photoisomerization and Cyclodimerization of Styryl Pyrimidinones*
- **Authors:** Mikhail V. Polynski, Zaruhi A. Hovasyan, Arthur A. Harutyunyan, Astghik A. Shahkhatuni, and Aleksan G. Shahkhatuni
- **Journal:** RSC Advances
- **DOI:** `[TO BE ADDED AFTER PUBLICATION]`
- **Archived release and raw data:** `https://doi.org/10.5281/zenodo.23167218`

Please cite the version of record when using the workflow or data.

## Article Abstract

Visible light-driven *E/Z* photoisomerization is an attractive strategy for controlling reactivity and biological activity; however, targeted design of heterocyclic bioactive compounds requires a solid understanding of how substituent effects govern productive isomerization and competing photochemical pathways. Here, a series of 6-methyl-2-styryl-pyrimidine-4(3*H*)-one derivatives with previously reported antifungal and antimicrobial activity of *E*-isomers were studied by LED-NMR spectroscopy under visible light irradiation at 405–450 nm. Pronounced substituent-dependent differences were observed in the efficiency of *E/Z* photoisomerization, the stability of the corresponding *Z*-isomers, and the formation of [2 + 2] cycloaddition products. Derivatives bearing electron-donating groups showed particularly efficient photoresponse, with up to 95% conversion to the *Z*-isomer within 1–3 h at moderate irradiation powers below 1 W. The resulting *Z*-isomers showed appreciable thermal stability, with half-lives ranging from several hours to several days. Quantum chemical calculations were used to rationalize the observed substituent effects and photoreaction outcomes. Among the computed descriptors, the oscillator strength *f*<sub>S1</sub>(*E*) was found to be the most informative. This descriptor represented the absorption-related contribution to photoconversion, and favorable transition intensities were calculated for derivatives bearing OMe-type donor-substituted phenyl groups attached to the vinylene fragment. These results establish substituted 6-methyl-2-styryl-pyrimidine-4(3*H*)-ones as a visible light-responsive scaffold in which photoconversion, *Z*-isomer stability, and competing cycloaddition can be tuned by substituent design.

## Additional Info

The main workflow was used to generate conformer ensembles and calculate conformer-level and Boltzmann-averaged descriptors. The calculated quantities include vertical excitation energies, excitation wavelengths, oscillator strengths, molecular polarizabilities, and Gibbs free energies of isomerization.

## Repository Contents

| File | Utility |
| --- | --- |
| `properties.py` | Main workflow for conformer processing, descriptor calculation, Boltzmann averaging, and thermochemistry. |
| `CSV_pooling.py` | Utility for pooling Boltzmann-averaged descriptors |
| `Check_dark_states.py` | Utility for checking if lower-energy dark states present |
| `track_missing_polarizability.py` | Utility for locating missing conformer polarizabilities and quantifying their effect on the ensemble average. |

## Workflow Summary

The calculations were performed separately for each stereoisomer.

1. Encode the required stereochemistry in an isomeric SMILES string.
2. Generate a three-dimensional structure with Open Babel and pre-optimize it with UFF.
3. Generate a conformer ensemble with CREST using GFN2-xTB and ALPB(DMSO).
4. Run `properties.py` in the directory containing `crest_conformers.xyz`.
5. Check for dark states and missing polarizabilities with the corresponding utility scripts (none were found).
6. Pool the ensemble rows from multiple molecular directories with `CSV_pooling.py`.

## Software Requirements

The calculations reported in the article involved external computational software:

- Open Babel 3.1.0
- CREST with GFN2-xTB and ALPB(DMSO)
- xtb 6.7.1
- `xtb4stda`
- sTDA-xTB (`stda` executable)
- ORCA 5.0.3
- [otherm](https://github.com/duartegroup/otherm)

## Selecting the Python Interpreter

The supplied scripts begin with a default hashbang:

```python
#!/usr/bin/env python3
```

For script execution, please provide the interpreter path corresponding to your (Conda - recommended) environment, for example:


```python
#!/Users/<USER>/anaconda3/envs/photochemistry/bin/python
```

On Windows, you can invoke the scripts using `python` or `conda run`.

## Configuring External Program Paths

Before running `properties.py`, replace the placeholder paths in the **External program paths** section:

```python
XTB_PATH = "/absolute/path/to/xtb"
XTB4STDA_PATH = "/absolute/path/to/xtb4stda"
STDA_PATH = "/absolute/path/to/stda"
ORCA_PATH = "/absolute/path/to/orca"
OTHERM_PATH = "/absolute/path/to/otherm.py"
```

## Preparing an Input

Use a separate, initially clean directory for every molecule or stereoisomer, for example:

```text
calculations/
└── molecule_E/
    └── init.smi
```

`init.smi` should contain the isomeric SMILES string for the required structure. Stereochemistry should be encoded before three-dimensional structure generation.

### 1. Generate and Pre-optimize a 3D Structure

Run the following command inside the molecular directory:

```bash
obabel -ismi init.smi -oxyz -O init.xyz --gen3d -h && \
obabel -ixyz init.xyz -oxyz -O init_min.xyz --minimize --ff UFF --steps 500 && \
mv init_min.xyz init.xyz
```

### 2. Generate a Conformer Ensemble

Run CREST with GFN2-xTB and ALPB(DMSO):

```bash
crest init.xyz --gfn2 --alpb dmso -ewin <energy_window> -T <num_cores> | tee crest.out
```

Replace `<energy_window>` with the required CREST energy window in kcal/mol and `<num_cores>` with the number of CPU cores. `tee` saves the CREST output in `crest.out` while displaying progress in the terminal. The subsequent modeling requires `crest_conformers.xyz`.

### 3. Running the Main Workflow

In the same directory, run:

```bash
python /path/to/repository/properties.py
```

The following operations are performed for every conformer:

- geometry optimization with GFN2-xTB, tight convergence, and ALPB(DMSO);
- calculation and parsing of the total energy and molecular polarizability;
- preprocessing with `xtb4stda`;
- sTDA-xTB calculation of excited states up to 4.5 eV;
- selection of the lowest-energy bright state with oscillator strength `f >= 0.05`;
- calculation of excitation energy, excitation wavelength, oscillator strength, etc.

If no state reaches `f >= 0.05`, the state with the largest oscillator strength is selected. A conformer for which no excitation block can be parsed is omitted from ensemble averaging and is reported through a warning.

For the lowest-energy conformer, ORCA and otherm are used for a QR-RHO thermochemistry calculation with a 1 M standard state. In the ORCA input generated by the present script, zero charge, singlet state, `XTB2`, `TightOpt`, `NUMFREQ`, and one process are specified.

## Main Output

The principal output is `stats.csv`. One row is provided for every conformer included in ensemble averaging, followed by a final row with:

```text
conf_id = ensemble
```

The `conformers/conf_*/` directories contain conformer geometries and program outputs. The `orca/` directory contains the thermochemistry calculation. Molden files generated during xTB optimization can require substantial disk space.


## Pooling Ensemble Descriptors

After all molecular calculations have been completed, arrange the directories as follows:

```text
calculations/
├── molecule_1/
│   └── stats.csv
├── molecule_2/
│   └── stats.csv
└── molecule_3/
    └── stats.csv
```

And run `CSV_pooling.py`. Only the direct subdirectories are scanned. The first row with `conf_id == "ensemble"` is extracted from each `stats.csv`, and the directory name is added as the molecular identifier. The resulting table is written to:

```text
ensembles_pooled.csv
```
