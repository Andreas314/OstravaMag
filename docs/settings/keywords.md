# Keywords
The code is controlled through arguments passed on the command line using the argparse module of Python. To list the available arguments, run
```
OstravaMag -h
```
Here, we describe the options in a bit more detail.
## Managing input
- --directory_name  
    Sets the name of the directory that contains the input files.
- --root_name  
    Sets the prefix of the input files; refer to the section about input for more information (not used for the TB2J input).
- --type  
    Sets the kind of input files that are read. Can be TB2J, OpenMX, KKR or GPAW.
- --max_R  
    Takes three positive integers as input. Sets the maximum index of the displacement vector to take into account during the calculation. Even if not all the exchange constants are included in the input within the given range, the calculation is still run.
- --exclude_elements  
    Takes multiple integers as input. Sets which sites (magnetic and non-magnetic) should be excluded from the calculation. The numbering starts at 0 and is the same as in the file containing the structure.
- --kpath_name  
    Sets the name (with path) of the file that contains the high-symmetry points in the Brillouin zone.
  

- --directory_name_qha  
    Sets the name of the directory that contains the input files for qha. If the qha calculation is turned on, the --directory_name option is disregarded.
## Flags
The following options turn on a certain type of calculation; they require no input.
- --compute_magnon_bands  
    Turns on the calculation of the magnon band structure.
  

- --compute_magnon_dos  
    Turns on the calculation of the magnon density of states.
- --compute_magnon_thermo  
    Turns on the calculation of the magnon-related thermodynamic properties (dos calculations must be turned on).
- --compute_magnon_projected_dos  
    Turns on the calculation of the magnon density of states projected onto the magnetic atoms.
- --compute_qha  
    Turns on the qha calculations; can be combined with other flags, but automatically turns the dos and thermodynamic calculations on.
- --compute_tc  
    Turns on the calculation of the critical temperature.
## Managing calculations
- --k_grid  
    Takes three integers as input. Sets the number of points of the k-grid in each direction (used for the dos and RPA calculations).
- --smearing  
    Size of the Gaussian smearing used to calculate the dos (in meV).
- --max_T  
    Maximum temperature to use in the thermodynamic calculations (always starts at 0 K).
- --num_T  
    Number of temperature points in the given range.
- --n_points  
    Number of energy points in the dos calculations; the upper and lower limits of the energies are chosen automatically.
- --epsilon  
    Epsilon parameter used in the RPA algorithm; refer to the reference for more information.
- --tol  
    Tolerance of the RPA algorithm, which is used as a stopping criterion.
- --max_it  
    Maximum number of iterations of the RPA algorithm.
- --n_cells  
    Number of unit cells in the supercell used. It is necessary to specify this during the **qha** calculations to obtain the correct scaling of the results.
- --mixing  
    Mixing value in the RPA algorithm