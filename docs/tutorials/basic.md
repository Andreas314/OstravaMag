# Tutorial #1-bcc Fe with jx
In this section, we will show, how to calculate magnon band structure, density of states, thermodynamic properties and Curie temperature of bcc Fe.
## Downloading the examples
All the examples that were used in the reference article can be run in the following manner
```
git clone git@github.com:Andreas314/OstravaMag.git
cd OstravaMag
python -m venv test-env
source test-env/bin/activate
pip install -e .
```
those commands clone the the repository, set up a virtual environment and install OstravaMag with all the dependencies inside the environment. The whole folder then can be deleted with
```
rm -rf OstravaMag
```
which erases the files along with the environment and the installation.
## Preparation of the input files
Iron is a classical example for calculations of the exchange constants. To obtain the structure and $J_{ij}$, follow the tutorial provided by OpenMX (we are interested in the **bulk** part of the tutorial):  
<https://www.openmx-square.org/openmx_man3.9/node125.html>   
To capture the output of *jx*, please run the calculations of the exchange constants on $n$ processors as follows
```
mpirun -np n ./jx Fe_Bulk_jx.scfout Fe_Bulk_jx.config > jx.out
```

Once those calculation are run, we need the *Fe_Bulk_jx.out* and *jx.out* files. Those files can be found in the directory with examples *OstravaMag/examples/OpenMX/Fe*. Create a directory called *Fe_tutorial* and copy the results of the calculations or the example inputs into this directory.
## Band structure
We begin by visualizing the band structure along the $\Gamma\rightarrow N\rightarrow P \rightarrow\Gamma\rightarrow H$ path. Create a file called *kpath.txt* inside the *Fe_tutorial* and put the following text defining the path in there
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
OstravaMag --directory_name Fe_tutorial --type OpenMX --root_name Fe_Bulk_jx --compute_magnon_bands --kpath_name Fe_tutorial/kpath.txt --max_R 2 2 2
```
This creates a directory called *OstravaMag_results_OpenMX* which contains the results. If everything went right, the band structure should be visualized in $magnon_bands.png$ and look like this

![Fe_bands](./figures/fe_bands.png)
One can crosscheck the results with tose obtained by the TB2J magnon code:  
<https://tb2j.readthedocs.io/en/latest/src/magnon_band.html>  
If you use the example data, you might also try running the calculation with "--max_R 4 4 4" option
## Density of states
Now that we see that the dispersion seem to be correct, we can proceed to calculate the dos and related quantities. Keep the files in the original folder and run
```
OstravaMag --directory_name Fe_tutorial --type OpenMX --root_name Fe_Bulk_jx --compute_magnon_dos --k_grid 20 20 20
```
This writes additional files into the folder with results. The total density looks like this (normalized to have integral equal to 1, the normalization is printed to standard output and should be equal to number of magnetic atoms in the calculation)

![Fe_dos_bad](./figures/fe_dos_bad.png)
We see that it is very spiky, so the default number of energy points and smearing is not sufficient. Run the calculation again raising both quantities and using denser grid
```
OstravaMag --directory_name Fe_tutorial --type OpenMX --root_name Fe_Bulk_jx --compute_magnon_dos --k_grid 30 30 30 --n_points 1000 --smearing 7
```
![Fe_dos_good](./figures/fe_dos_good.png)
This looks better, the quality can be improved by running the calculations with more k-points, energy points and smaller smearing. In this case, the calculation takes a few second, but could get more demanding on more complex structures.  
## Thermodynamic quantities
The dos now looks reasonable, we would wish to see the free energy, internal energy, entropy and constant-volume heat capacity. Run the following command utilizing the correct setup for dos
```
OstravaMag --directory_name Fe_tutorial --type OpenMX --root_name Fe_Bulk_jx --compute_magnon_dos --k_grid 30 30 30 --compute_magnon_thermo --n_points 1000 --smearing 7 --max_T 5000 --num_T 5001
```
This computed the quantities up to 5000 K with 1 K resolution (5001 points, 0 K points must be counted in). The result should look like this  
![Fe_thermo](./figures/fe_thermo.png)
This result should be quite robust with the setup of the dos calculation. We see that the high-temperature limit of $C_V$ appears at around 4000 K and $F$, $E$ and $S$ meet the thermodynamic relationship.