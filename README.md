# OstravaMag

A linear spin wave theory (LSWT) code for computing magnon properties from first-principles exchange interaction constants.

## Overview

OstravaMag computes magnon band structures, density of states (DOS), and thermodynamic quantities from exchange interaction parameters *J* obtained via:

- **TB2J** (default)
- **OpenMX**
- **GPAW**
- **SPRKKR** (experimental)

The code solves the LSWT Hamiltonian on a *k*-point grid or along a specified *k*-path and integrates the resulting magnon spectrum to yield thermodynamic observables including internal energy, free energy, entropy, and heat capacity.

## Installation
To obtain working installation on your system:
```bash
pip install git+ssh://git@github.com:Andreas314/OstravaMag.git
```
To run the examples, we recomend clonning the repository and running the code inside a virtual environmet:
```bash
git clone git@github.com:Andreas314/OstravaMag.git
cd OstravaMag
python -m venv test-env
source test-env/bin/activate
pip install -e .
```
The whole installation can then be erased with the cloned repository.
## Usage
See all the available options with:
```bash
OstravaMag -h
```
## Command Line Options

| Option | Description | Default |
|--------|------------|---------|
| `-h, --help` | Show help message and exit | — |
| `--max_R MAX_R MAX_R MAX_R` | Maximal R to consider | `2 2 2` |
| `--directory_name DIRECTORY_NAME` | Name of the directory with the input file | `.` |
| `--root_name ROOT_NAME` | Root name of the output files from OpenMX/GPAW | `none` |
| `--type TYPE` | Type of the input file | `TB2J` |
| `--exclude_elements [EXCLUDE_ELEMENTS ...]` | Indices of elements to exclude from the calculation | `none` |
| `--compute_magnon_bands` | Compute magnon band structure | Off |
| `--compute_magnon_dos` | Compute magnon density of states | Off |
| `--compute_magnon_thermo` | Compute thermodynamical quantities | Off |
| `--compute_magnon_projected_dos (experimental)` | Compute DOS projected on included atoms | Off |
| `--kpath_name KPATH_NAME` | File with k-path for magnon band structure calculation | None |
| `--k_grid K_GRID K_GRID K_GRID` | Number of k-points in each BZ direction | `8 8 8` |
| `--smearing SMEARING` | Gaussian smearing size (meV) | `5 meV` |
| `--max_T MAX_T` | Maximum temperature for thermodynamical calculations (K) | `300 K` |
| `--num_T NUM_T` | Number of temperature points | `300` |
| `--n_points N_POINTS` | Number of points in DOS | `50` |
| `--compute_qha` | Compute thermal properties in quasi-harmonic approximation *(not working yet)* | Off |
| `--directory_name_qha DIRECTORY_NAME_QHA` | Directory for QHA calculations | `.` |

## Input file formats

### TB2J
Place the TB2J output directory (containing `TB2J.pickle`) in `--directory_name`. Optionally set `--type TB2J`.

### OpenMX
Set `--type OpenMX`. The input file in the directory `--directory_name` are expected to be the result of the scf calculation `"--root_name".out` and the result of the `jx` calculation `jx.out`.

### GPAW
Set `--type GPAW`. The input file in the directory `--directory_name` are expected to be the result of the scf calculation `"--root_name".txt`, the numpy file containing the exchange coupling constants `j_"--root_name".npy` and the numpy file containing the mesh used for the calculations `q_"--root_name".npy`. Note that the code expects the calculations to be carried out in the upper half of the BZ, where $g_z\geq0$.

### SPRKKR
Set `--type KKR`. The input file in the directory `--directory_name` are expected to be the potential file with structural information and magnetic moments `"--root_name".pot` and the result of the exchange coupling constant calculations `"--root_name"_JXC_XCPLTEN_Jij.dat`.

### *k*-path file (`--kpath_name`)
A plain-text file specifying high-symmetry points and labels, one per line, the first line contains number of points along each path:

```
20
0 0 0 $\Gamma$
0.5 0 0 M

0.5 0 0 M
0.33 0.33 0 K

0.33 0.33 0 K
0 0 0 $\Gamma$
```

## Output files

All output files are written to `OstravaMag_results_` with the suffix `--type` in the directory of the input files `--directory_name`. All the results are contained in a data file `.dat` with a plot in a `.png` file.

| File | Contents |
|---|---|
| `magnon_bands.dat/.png` | Magnon band energies along the *k*-path |
| `magnon_tdos.dat/.png` | Total magnon DOS as a function of energy |
| `magnon_pdos.dat/.png` | Magnon DOS of each branch as a function of energy|
| `magnon_projected_dos.dat/.png` | Magnon DOS projected on each atom as a function of energy|
| `thermo.dat/.png` | U, F, S, Cv as a function of temperature |
| `energy_grid.dat` | Magnon energies in each point of the *k*-grid|
| `Jr_GPAW.dat` | J constants in real space used in GPAW calculations|
## Examples
### DOS

```bash
OstravaMag --directory_name examples/GPAW/Fe_bcc/ --type GPAW --root_name Fe_bcc --compute_magnon_dos --max_R 4 4 4 --n_points 1000 --smearing 6
```

### Magnon band structure along a *k*-path

```bash
OstravaMag --directory_name examples/GPAW/Fe_bcc/ --type GPAW --root_name Fe_bcc --compute_magnon_bands --kpath_name ./examples/GPAW/Fe_bcc/kpath.txt --max_R 4 4 4
```

