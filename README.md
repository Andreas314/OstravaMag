# OstravaMag

A linear spin wave theory (LSWT) code for computing magnon properties from first-principles exchange interaction constants.

## Overview

OstravaMag computes magnon band structures, density of states (DOS), and thermodynamic quantities from exchange interaction parameters *J* obtained via:

- **TB2J** (default)
- **OpenMX**
- **GPAW**

The code solves the LSWT Hamiltonian on a *k*-point grid or along a specified *k*-path and integrates the resulting magnon spectrum to yield thermodynamic observables including internal energy, free energy, entropy, and heat capacity.

## Installation

```bash
*ADD PIP INSTALLATION LATER!*
pip install -r requirements.txt
```

## Usage

```bash
OstravaMag [options]
```

### Basic example

```bash
OstravaMag --directory_name ./examples/Fe2Ti/FM/TB2J --type TB2J \
           --k_grid 20 20 20 --compute_magnon_dos --compute_magnon_thermo
```

### Magnon band structure along a *k*-path

```bash
OstravaMag --directory_name ./Fe2Ti --root_name Fe2Ti \
           --compute_magnon_bands --kpath_name kpath.txt \
           --smearing 3
```

## Options

| Option | Default | Description |
|---|---|---|
| `--directory_name` | `.` | Directory containing the input file |
| `--root_name` | `none` | Root name for output files and the OpenMX/GPAW input file |
| `--type` | `TB2J` | Input file type: `TB2J`, `OpenMX`, or `GPAW` |
| `--max_R` | `2 2 2` | Maximum real-space range of exchange interactions to include (in unit cells) |
| `--exclude_elements` | none | Indices of atomic sites to exclude from the calculation |
| `--compute_magnon_bands` | off | Compute magnon band structure along `--kpath_name` |
| `--compute_magnon_dos` | off | Compute magnon density of states on the `--k_grid` |
| `--compute_magnon_thermo` | off | Compute thermodynamic quantities (requires DOS) |
| `--kpath_name` | none | File defining the *k*-path for band structure (required with `--compute_magnon_bands`) |
| `--k_grid` | `8 8 8` | Number of *k*-points in each BZ direction for DOS/thermo |
| `--smearing` | `5` | Gaussian smearing width in meV for DOS broadening |
| `--max_T` | `300` | Maximum temperature in K for thermodynamic quantities |
| `--num_T` | `300` | Number of temperature points in the interval [0, `--max_T`] |

## Input file formats

### TB2J
Place the TB2J output directory (containing `TB2J.pickle`) in `--directory_name`. Optionally set `--type TB2J`.

### OpenMX / GPAW
Provide the root name of the output via `--root_name` and set `--type OpenMX` or `--type GPAW` accordingly.

### *k*-path file (`--kpath_name`)
A plain-text file specifying high-symmetry points and labels, one per line, the first line contains number of points along each path:

```
20
0 0 0 G
0.5 0 0 M

0.5 0 0 M
0.33 0.33 0 K

0.33 0.33 0 K
0 0 0 G
```

## Output files

All output files are written to `OstravaMag_results_` with the suffix `--type`.

| File | Contents |
|---|---|
| `magnon_bands.dat` | Magnon band energies along the *k*-path |
| `magnon_dos.dat` | Magnon DOS as a function of energy |
| `thermo.dat` | U, F, S, Cv as a function of temperature |
| `magnons_on_grid.dat` | Magnon energies in each point of the *k*-grid|

## Thermodynamic quantities

Computed within the harmonic (LSWT) approximation via integration over the magnon DOS *g*(E):

| Quantity | Formula |
|---|---|| Internal energy | $U = \int g(E)\, E\left(\frac{1}{2} + n(E,T)\right) dE$ |
| Free energy | $F = \int g(E)\left[\frac{E}{2} + k_BT \ln\!\left(1 - e^{-E/k_BT}\right)\right] dE$ |
| Entropy | $S = k_B \int g(E)\left[\frac{x}{e^x - 1} - \ln\!\left(1 - e^{-x}\right)\right] dE$ |
| Heat capacity | $C_v = k_B \int g(E)\, x^2\, n(E,T)\left[n(E,T)+1\right] dE$ |

where $x = E / k_BT$ and $n(E,T) = 1/(e^x - 1)$ is the Bose-Einstein occupation.
