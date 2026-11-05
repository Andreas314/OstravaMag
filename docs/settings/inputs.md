# Input files
This section offers a manual on how to prepare the input files before running the code. OstravaMag has four interfaces-TB2J, jx, SPR-KKR and GPAW-each of them requires different kind of files to be able to load the necessary parameters. Additionally, one can create a TB2J style input from the first few known exchange constants through a simple input text file. Finally, the k-path must be specified in an external text file.  

**Important Note:** When "input directory" is referred, the directory  given by the "--directory_name" argument is meant. Similarly, when files are labeled by a name containing "*", the name given by the "--root_name" argument should be substituted into the filename. 
## TB2J interface
Aside from human readable files and plots, TB2J also outputs a Python readable file *TB2J.pickle*. This file contains all the necessary data to run the magnon calculation. It cannot be renamed and must be placed in the input directory.  
&emsp;Although not thoroughly tested (yet). When TB2J option is used, the magnetic moments are read as vectors, the calculation for non-collinear structures should be possible and preliminary results seem to be working well.  

For more information about **TB2J**, please refer to:  
<https://tb2j.readthedocs.io/en/latest/index.html>
## OpenMX's jx interface
To run the calculation based on OpenMX's DFT calculations and jx's calculations of the exchange constants, two input files must be placed in the input directory. The structural information is read from the file, which contains the results of the DFT calculations and should be called *.out. The jx code writes the exchange constants into the standard output. This should be captured in a file called *jx.out*, which is then read by the interface scripts.

**Note:** This interface yet employs no symmetry for OpenMX, all the used pairs must therefore be calculated using jx or a custom script must be used.

Only version 3.9 of OpenMX was tested. For more information, please refer to:  
<https://www.openmx-square.org/openmx_man3.9/index.html>  
Specifically for the jx section:  
<https://www.openmx-square.org/openmx_man3.9/node119.html>
## SPR-KKR interface
The structural information for SPR-KKR is read from the *.pot potential file. The exchange constants are written in the *_JXC_XCPLTEN_Jij.dat file, which contains all the exchanges.

For more information, please refer to the **SPR-KKR** manual:  
<https://www.sprkkr.org/manual/>

## GPAW interface
In the case of GPAW, $J_{ij}$ are output directly in reciprocal space. The exchange constants and a list of points on which they are calculated should be saved in a file called j_*.npy and q_*.npy. The structure is read from the output, that GPAW writes to the standard input. It should be contained in a file called *.txt.  
&emsp; This interface takes advantage of the fact that $J_{ij}(\boldsymbol{k})=J_{ij}^*(\boldsymbol{-k})$. Therefore, the exchange constants should only be calculated for k-points with a non-negative $z$ coordinate (upper part of the Brillouin zone).

**Note:** Here, is is hard to identify, for which sites in the unit cells the exchange constants were calculated. So in case there are atoms which were not included in the calculations of the exchange constants, it is necessary to exclude them explicitly using the "--exclude_elements" keyword.

For more information on the calculations of the exchange constants using **GPAW**, please visit:
<https://gpaw.readthedocs.io/tutorialsexercises/magnetic/mft/mft.html>
## Minimum input
To perform the calculation, it is enough to know the lattice vectors, positions of the magnetic atoms, their magnetic moments and $n$ nearest interactions. Those information can be put in a text file as follows
```txt
lattice
    4.180022      2.090011      2.090011
    2.090011      4.180022      2.090011
    2.090011      2.090011      4.180022

atoms fractional                      # symbol  x y z   mx my mz  (mu_B)
    Ni  0.00  0.00  0.00     0  0   2
    Ni  0.50  0.50  0.50     0  0   -2

exchange                              # distance [A]   J [eV]
    2.9557   0.7
    4.1800   -9.5
```
This file is composed of three segments, each labeled by a keyword. The
segment labeled by the keyword "lattice" contains three lines, with one lattice
vector per line. The next keyword "atoms" must be followed by "fractional" or
"cartesian", depending on what type of coordinates we use in this segment. Each
line then contains a label, followed by three coordinates and three components
of the magnetic moment. The number of lines matches the number of atoms in
the calculation. The last segment starts with "exchange" and contains the distance
between the sites followed by the exchange constant (in meV).  
&emsp; This file serves as an input for the *OstravaMag_PrepareInput* script, which transforms it into the *TB2J.pickle* file, which can be read by *OstravaMag* as described above. How much the distances between the sites in the input file can differ from the ones found by the code is specified by the "--tol" keyword. If it happens that a distance is off, possible distances found are output to the standard output. In case the distances are omitted and just the exchange constants are written to the input, they are assigned to the $n$ lowest shells found by the code.  

**Note:** One limitation of this script is when there are two $J_{ij}$ with different values but the same distance. In this case, one might try breaking this symmetry by varying the lattice vectors/atomic positions slightly.

## kpath file
To plot the band structure, a definition of the kpath must be given. For OstravaMag, it is written in an external file, which position and name is given by the "--kpath_name" argument. An example of this file could look like this
```txt
20
0 0.5 0.5 X
0 0 0 $\Gamma$

0.5 0.5 0.5 L
0 0 0 $\Gamma$
```
This file defines a $X\rightarrow\Gamma|L\rightarrow \Gamma$ path for an fcc lattice with 20 points on each segment. The first line contains number of points on each segment. Each subsequent pair of lines contains the starting and finishing point of the segment. Each line contains relative coordinates of a high-symmetry point and label at the end of the line. Blank lines are allowed anywhere. To use greek letters in the plots, use labels as indicated in the example. 

## QHA
### Magnon contributions alone
First, let us describe how to obtain the magnon contributions to the Gibbs free energy and constants-pressure heat capacity.