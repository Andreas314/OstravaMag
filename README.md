# OstravaMag

A linear spin wave theory (LSWT) code for computing magnon properties from first-principles exchange interaction constants.  
## Overview

OstravaMag computes magnon band structures, density of states (DOS), and thermodynamic quantities from exchange interaction parameters *J* obtained via:

- **TB2J** (default)
- **OpenMX**
- **GPAW**
- **SPRKKR**

The code solves the LSWT Hamiltonian on a *k*-point grid or along a specified *k*-path and integrates the resulting magnon spectrum to yield thermodynamic observables including internal energy, free energy, entropy, and heat capacity. Also allows to calculate the critical temperature of the system and magnon free energy to phonopy's input for Quasi harmonic approximation calculations.  
Download the code with all the examples and install it inside a virtual environment with the following set of commands
```
git clone git@github.com:Andreas314/OstravaMag.git
cd OstravaMag
python -m venv test-env
source test-env/bin/activate
pip install -e .
```
The virtual environment must then be turned on upon each usage with
```
source test-env/bin/activate
```
The whole installation can be easily erased by erasing the directory.
## Documentation
Detailed description can be found in the documentation, which can be found in :
<https://ostravamag.readthedocs.io/en/latest/>  
or compiled to a pdf manual by running 
```
./make_docs.sh
```
in the docs folder. It uses pandoc and latex extensions, which can be installed using
```
sudo apt install -y pandoc \
  texlive-xetex texlive-latex-recommended texlive-latex-extra \
  texlive-fonts-recommended texlive-fonts-extra texlive-science \
  texlive-plain-generic lmodern fonts-dejavu fonts-libertinus \
  imagemagick
```
