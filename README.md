# OstravaMag

A linear spin wave theory (LSWT) code for computing magnon properties from first-principles exchange interaction constants.

## Overview

OstravaMag computes magnon band structures, density of states (DOS), and thermodynamic quantities from exchange interaction parameters *J* obtained via:

- **TB2J** (default)
- **OpenMX**
- **GPAW**
- **SPRKKR**

The code solves the LSWT Hamiltonian on a *k*-point grid or along a specified *k*-path and integrates the resulting magnon spectrum to yield thermodynamic observables including internal energy, free energy, entropy, and heat capacity. Also allows to calculate the critical temperature of the system and magnon free energy to phonopy's input for Quasi harmonic approximation calculations.

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
For all the given functionalities, one has to specify maximal lattice vector to consider via `--max_R`, type of the input files via `--type`, working directory via `--directory_name` (or `--directory_name_qha` for the interafce with phonopy) and the given switch keywords specifying what kind of calculations are to be performed. Further options influencing the calculations are specified below.

The code is able to perform the following calculations:

- Bandstructure calculations on a given k-path (`--kpath_name` has to be specified)
- Density of states calculations (`--k_grid`, `--smearing` and `--n_points` can be specified, the interval of energies is chosen automatically, GPAW uses k-grid from an input file-see below)
- Dos projected on the individual magnetic atoms (dos calculations have to be turned on, same specifications apply as for dos calculations)
- Internal energy, free energy, entropy and heat capacity calculations (dos calculations have to be turned on; `--max_T` and `--num_T` can be specified)
- Critical temperature calculations in MFA and RPA (`--k_grid`, `--epsilon`, `--tol` and `--max_it` can be specified)
  
The code also automatically exclude atoms with sufficiently small magnetic moment. One can exclude more atoms via `--exclude_elements`. The indexing starts from 0.

It also allows to add the magnon contributions of free energy into the fe-v.dat input file of the phonopy's qha calculations. For this, one has to put folders containing the input files for multiple volumes into a directory specified by `--directory_name_qha` along the fe-v.dat file. The volumes in the fe-v.dat file have to be subset of the volumes for which the exchange constants were calculated. Otherwise, the exceeding volumes are cut off. Dos and thermodynamic properties are automatically calculated for each volume. Other options can be used alongside this option as well, but the temperature range is taken from the fe-v.dat file.
## Input file formats
The code expects all the input files to be the results of collinear calculations. 
### TB2J
Place the TB2J output directory (containing `TB2J.pickle`) in `--directory_name`. Optionally set `--type TB2J`.

### OpenMX
Set `--type OpenMX`. The expected input files in the directory `--directory_name` are the result of the scf calculation `"--root_name".out` and the result of the `jx` calculation `jx.out`.

### GPAW
Set `--type GPAW`. The expected input file in the directory `--directory_name` are the result of the scf calculation `"--root_name".txt`, the numpy file containing the exchange coupling constants `j_"--root_name".npy` and the numpy file containing the mesh used for the calculations `q_"--root_name".npy`. Note that the code expects the calculations to be carried out in the upper half of the BZ, where $g_z\geq0$.

### SPRKKR
Set `--type KKR`. The expected input file in the directory `--directory_name` are the potential file with structural information and magnetic moments `"--root_name".pot` and the result of the exchange coupling constant calculations `"--root_name"_JXC_XCPLTEN_Jij.dat`.

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

All output files are written to `OstravaMag_results_` with the suffix `--type` in the directory of the input files `--directory_name`. All the results are contained in a data file `.dat` possibly with a plot in a `.png` file.

| File | Contents |
|---|---|
| `magnon_bands.dat/.png` | Magnon band energies along the *k*-path |
| `magnon_tdos.dat/.png` | Total magnon DOS as a function of energy |
| `magnon_pdos.dat/.png` | Magnon DOS of each branch as a function of energy|
| `magnon_projected_dos.dat/.png` | Magnon DOS projected on each atom as a function of energy|
| `thermo.dat/.png` | U, F, S, Cv as a function of temperature |
| `energy_grid.dat` | Magnon energies in each point of the *k*-grid|
| `Jr.dat` | J constants in real space (for OpenMX/GPAW input)|
| `fe-v-OstravaMag.dat`|Modified input fe-v.dat file with magnon free energy added|
## Examples
### DOS

```bash
OstravaMag --directory_name examples/GPAW/Fe_bcc/ --type GPAW --root_name Fe_bcc --compute_magnon_dos --max_R 4 4 4 --n_points 1000 --smearing 6
```

### Magnon band structure along a *k*-path

```bash
OstravaMag --directory_name examples/GPAW/Fe_bcc/ --type GPAW --root_name Fe_bcc --compute_magnon_bands --kpath_name ./examples/GPAW/Fe_bcc/kpath.txt --max_R 4 4 4
```

