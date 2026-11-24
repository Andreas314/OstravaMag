# Tutorial #1-bcc Fe with jx
In this section, we will show how to calculate the magnon band structure, density of states, thermodynamic properties and Curie temperature of bcc Fe.
## Downloading the examples
The setup for all the examples used in the reference article is as follows
```
git clone git@github.com:Andreas314/OstravaMag.git
cd OstravaMag
python -m venv test-env
source test-env/bin/activate
pip install -e .
```
These commands clone the repository, set up a virtual environment and install OstravaMag with all of its dependencies inside that environment. The whole folder can then be deleted with
```
rm -rf OstravaMag
```
which erases the files, along with the environment and the installation.
## Preparation of the input files
Iron is a classic example for calculations of the exchange constants. To obtain the structure and $J_{ij}$, follow the tutorial provided by OpenMX (we are interested in the **bulk** part of the tutorial):  
<https://www.openmx-square.org/openmx_man3.9/node125.html>   
To capture the output of _*jx*_, please run the calculations of the exchange constants on $n$ processors as follows

```
mpirun -np n ./jx Fe_Bulk_jx.scfout Fe_Bulk_jx.config > jx.out
```
  

Once those calculations have been run, we need the _*Fe_Bulk_jx.out*_ and _*jx.out*_ files. These files can be found in the examples directory, _*OstravaMag/examples/OpenMX/Fe*_. Create a directory called _*Fe_tutorial*_ and copy the results of the calculations, or the example inputs, into this directory.
## Band structure
We begin by visualizing the band structure along the $\Gamma\rightarrow N\rightarrow P \rightarrow\Gamma\rightarrow H$ path. Create a file called _*kpath.txt*_ inside _*Fe_tutorial*_ and put the following text defining the path into it
```txt
50
0 0 0 $\Gamma$
0 0 0.5 N
  

0 0 0.5 N
0.25 0.25 0.25 P
  

0.25 0.25 0.25 P
0 0 0 $\Gamma$
  

0 0 0 $\Gamma$
0.5 -0.5 0.5 H
```
We can now run the calculation with the following command:
```
OstravaMag --directory_name Fe_tutorial --type OpenMX --root_name Fe_Bulk_jx \
--compute_magnon_bands --kpath_name Fe_tutorial/kpath.txt --max_R 2 2 2
```
This creates a directory called _*OstravaMag_results_OpenMX*_, which contains the results. If everything went well, the band structure should be visualized in $magnon_bands.png$ and should look like this
  

![Resulting band structure of bcc Fe from the OpenMX input](./figures/fe_bands.png)

One can cross-check the results with those obtained by the TB2J magnon code:  
<https://tb2j.readthedocs.io/en/latest/src/magnon_band.html>  
If you use the example data, you might also try running the calculation with the "--max_R 4 4 4" option.
## Density of states
Now that we see that the dispersion seems to be correct, we can proceed to calculate the dos and the related quantities. Keep the files in the original folder and run
```
OstravaMag --directory_name Fe_tutorial --type OpenMX --root_name Fe_Bulk_jx \
--compute_magnon_dos --k_grid 20 20 20
```
This writes additional files into the results folder. The total density looks like this (normalized so that the integral is equal to 1; the normalization is printed to standard output and should be equal to the number of magnetic atoms in the calculation)
  

![Magnon DOS of bcc Fe for too few points](./figures/fe_dos_bad.png)

We see that it is very spiky, so the default number of energy points and the default smearing are not sufficient. Run the calculation again, raising both quantities and using a denser grid
```
OstravaMag --directory_name Fe_tutorial --type OpenMX --root_name Fe_Bulk_jx \
--compute_magnon_dos --k_grid 30 30 30 --n_points 1000 --smearing 7
```

![Magnon DOS of bcc Fe with enough points and larger smearing](./figures/fe_dos_good.png)

This looks better; the quality can be improved by running the calculations with more k-points, more energy points and a smaller smearing. In this case, the calculation takes a few seconds, but it could become more demanding for more complex structures.
## Thermodynamic quantities
The dos now looks reasonable, and we would like to see the free energy, internal energy, entropy and constant-volume heat capacity. Run the following command, utilizing the correct setup for the dos
```
OstravaMag --directory_name Fe_tutorial --type OpenMX --root_name Fe_Bulk_jx \
--compute_magnon_dos --k_grid 30 30 30 --compute_magnon_thermo --n_points 1000 \
--smearing 7 --max_T 5000 --num_T 5001
```
This computes the quantities up to 5000 K with a 1 K resolution (5001 points; the 0 K point must be counted in). The result should look like this

![Thermodynamic quantities from the calculated DOS](./figures/fe_thermo.png)

This result should be quite robust with respect to the setup of the dos calculation. We see that the high-temperature limit of $C_V$ appears at around 4000 K and that $F$, $E$ and $S$ satisfy the thermodynamic relationship.