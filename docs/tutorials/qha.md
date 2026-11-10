# Tutorial #3-FeRh (AFM phase) qha calculation
In this section, we will show, how to calculate qha with magnon contributions using OstravaMag. As an example, we will use AFM phase of FeRh with a fcc supercell (two times larger than the unit cell) with 2 magnetic atoms.  
First, make a directory called *FeRh_tutorial* and copy the example  OpenMX calculations into it 
```
cp -r examples/QHA/FeRh/AFM/* FeRh_tutorial/
```
This contains calculations of exchange constants on multiple volumes as well as a file called "OstravaMaf_qha.txt" containing the total energies from DFT calculations in the ascending order with respect to volume. With the right naming convention of the folders, this file can be simply obtained from the OpenMX output files with the following command
```
grep "Utot\." V_*/FeRh.out| tr -s " " | cut -d " " -f 3 | /
awk '{$1=$1*27.211386245; printf "%.7f\n",$1}' > OstravaMag_qha.txt
```
Important here is that the energies must be in eV and OpenMX outputs them in Hartree. The naming of the folders does not matter, but the total energies must be sorted in the right order. Once we have that, the qha calculations can be rn with the following minimal command
```
OstravaMag --compute_qha --directory_name_qha FeRh_tutorial/ --type OpenMX /
--root_name FeRh --n_cells 2 --max_T 1000 --num_T 1001
```
It is important to specify here the number of cells in the supercell. The results are stored the *OstravaMag_results_qha* directory. We can check the quality of the fits of the free energy for a few temperatures in the *qha_fits.png* file

![FeRh_fits](./figures/ferh_fits.png)

We see that even with the default dos setting, we obtained a reasonable polynomial fit. The constant-pressure heat capacity should look like this

![FeRh_fits](./figures/ferh_cp.png)
