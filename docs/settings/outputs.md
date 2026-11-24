# Output files
In this section, we provide a short overview of the output files and their contents. All the data are written into data files in an OstravaMag_results folder (with a suffix of the interface used) created in the input directory. Additionally, the band structure, density of states and thermodynamic properties are plotted for a quick analysis of the outputs. The following data files are saved:
  

- magnon_bands.dat  
    Contains the results of the band structure calculations.
- magnon_bands_ticks.dat  
    Contains the positions of the labels of the high-symmetry points on the k-path.
- magnon_tdos.dat  
    Contains the total density of states. The data are normalized so that the integral is equal to 1.
- magnon_pdos.dat  
    Contains the partial density of states (dos of each magnon branch).
- magnon_projected_dos.dat  
    Contains the density of states after it has been projected onto the magnetic atoms. The order of the magnetic sites is the same as the order used in the calculation.
- thermo.dat  
    Contains the free energy, entropy, constant-volume heat capacity and internal energy. The quantities are output per magnetic atom.
- energy_grid.dat  
    Contains the magnon energies at each k-point.
- Jr.dat  
    Contains the lattice vectors, the positions of the atoms in the unit cell and the exchange constants in the TB2J conventions, sorted according to the distance between two sites.
- fe-v-OstravaMag.dat  
    Contains the data from the input fe-v.dat phonopy file to which the free energy contribution of magnons has been added.
- qha_results.dat  
    Contains the results of the qha calculations ($V_0$, $C_p$ and $G$). The quantities are output per formula unit, similarly to the phonopy output.
  

The results of the calculations of the critical temperature are written to standard output, along with information about the convergence of the RPA algorithm.