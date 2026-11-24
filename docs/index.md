# OstravaMag: A Python package for calculating the properties of magnons
  

OstravaMag is a lightweight open-source Python package aimed at providing a simple way to calculate magnon-based properties, exchange constants in a unified convention, magnon contributions to the quasi-harmonic approximation (QHA) along with a phonopy interface, and critical temperatures of magnetic bulk crystals. The main emphasis is on collinear calculations, but the code can handle non-collinear calculations as well, to some extent.  
&emsp;The code accepts input from calculations done by other packages (TB2J [1], OpenMX's jx [2], GPAW [3], KKR [4]) and can be easily extended to be compatible with other codes as well. Additionally, a script is provided that turns the results of experimental measurements or total-energy mappings into an input that can be read automatically.  
&emsp;The magnon dispersion is calculated from the bosonic Bogoliubov-de Gennes. The critical temperatures are calculated in the Mean-Field Approximation (MFA) and Random Phase Approximation (RPA). All the data outputs are accompanied by preliminary plots for a fast evaluation of the results. Moreover, a file containing the exchange constants based on the distance between sites in a unifying format is written, since each code uses a different formalism and may not provide a human-readable way.  
&emsp;Note: Most of the calculations that have been performed to test this code were done on an average office laptop and did not take more than a few minutes.
  
  

Source code can be found in:  
<https://github.com/Andreas314/OstravaMag>
  

The installation can be carried out with
```text
pip install git+ssh://git@github.com/Andreas314/OstravaMag.git
```

```{toctree}
:maxdepth: 1
:caption: Theory

theory
```

```{toctree}
:maxdepth: 1
:caption: Settings

settings/inputs
settings/keywords
settings/outputs
```

```{toctree}
:maxdepth: 1
:caption: Tutorials

tutorials/basic
tutorials/minimp
tutorials/qha
```

```{toctree}
:maxdepth: 1
:caption: References

references
```
