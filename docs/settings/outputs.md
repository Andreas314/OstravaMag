# Output files
In this section, we provide a short overview of the output files and their content. All the data are written into data files into a *OstravaMag_results* folder created in the input directory. Additionally, band structure, density of states and thermodynamic properties are plotted for a quick analysis of the outputs. The following data files are saved:

- magnon_bands.dat  
    Contains the results of the band structure calculations.
- magnon_bands_tics.dat  
    Contains the position of the labels of the high-symmetry points on the k-path. 
- magnon_tdos.dat  
    Contains the total density of states. The data are normalized to have the integral equal to 1.
- magnon_pdos.dat  
    Contains the partial density of states (dos of each magnon branch).
- magnon_projected_dos.dat  
    Contains the density of states after it was projected on the magnetic atoms. The order of the magnetic sites is the same as the order in the calculation. 
- thermo.dat  
    Contains free energy, entropy, constant-volume heat capacity and internal energy. The quantities are output per magnetic atom.
- energy_grid.dat  
    Contains magnon energies in each k-point.
- Jr.dat  
    Contains the lattice vectors, the positions of the atoms in the unit cells and the exchange constants in the TB2J conventions sorted according to the distance between two sites. 
- fe-v-OstravaMag.dat  
    Contains data from the input fe-v.dat phonopy file to which free energy contribution of magnons was added.
- qha_results.dat
    Contains results of the qha calculations ($V_0$, $C_p$ and $G$). The quantities are output per formula unit similarly to phonopy output.

The results of the calculations of the critical temperature are written into the standard input along with the information about the convergence of the RPA algorithm. 