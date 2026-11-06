# OstravaMag: A Python package for calculating the properties of magnons

OstravaMag is lightweight open-source a Python package aimed to provide a simple way to calculate magnon-based properties, exchange constants in an unified convention, magnon contributions to quasi-harmonic approximation (QHA) along with a phonopy interface, and critical temperatures of magnetic bulk crystals. Main emphasis is put on collinear calculations, but the code can handle non-collinear calculations as well, to some extend.  
&emsp;The code accepts input from calculations done by other packages (TB2J [1], OpenMX's jx [2], GPAW [3], GPAW [4]) and can be easily extended to be compatible with other codes as well. Additionally, a script is provided, which turns results of experimental measurements or total-energy mappings into an input that can be automatically read.  
&emsp;The magnon dispersion is calculated from the the bosonic Bogoliubov-de Gennes. The critical temperatures are calculated in the Mean-Field Approximation (MFA) and Random Phase Approximation (RPA). All the data outputs are supported by preliminary plots for a fast evaluation of the results. Moreover, a file containing the exchange constants based on the distance of sites in a unifying format is written, since each code uses a different formalism and may not provide a human readable way.  
&emsp;Note: Most of the calculations that has been performed to test this code were done on an average office laptop and did not take more than few minutes.


Source code can be found in:  
 <https://github.com/Andreas314/OstravaMag>

The installation can be carries out with
```text
pip install git+ssh://git@github.com:Andreas314/OstravaMag.git
```

   