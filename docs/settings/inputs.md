# Input files
This section offers a manual on how to prepare the input files before running the code. OstravaMag has four interfaces: TB2J, jx, SPR-KKR and GPAW; each of which requires different kinds of files in order to load the necessary parameters. Additionally, one can create a TB2J-style input from the first few known exchange constants through a simple input text file. Finally, the k-path must be specified in an external text file.

**Important Note:** When the "input directory" is referred to, the directory given by the "--directory_name" argument is meant. Similarly, when files are labeled with a name containing "\*", the name given by the "--root_name" argument should be substituted into the filename.

## TB2J interface
Aside from human-readable files and plots, TB2J also outputs a Python-readable file *TB2J.pickle*. This file contains all the necessary data to run the magnon calculation. It cannot be renamed and must be placed in the input directory.  
&emsp;Although not thoroughly tested (yet), when the TB2J option is used, the magnetic moments are read as vectors; calculations for non-collinear structures should therefore be possible, and preliminary results seem to be working well.

For more information about **TB2J**, please refer to:  
<https://tb2j.readthedocs.io/en/latest/index.html>

## OpenMX's jx interface
To run the calculation based on OpenMX's DFT calculations and jx's calculations of the exchange constants, two input files must be placed in the input directory. The structural information is read from a file that contains the results of the DFT calculations and should be called  *\*.out*. The jx code writes the exchange constants to standard output. These should be captured in a file called *jx.out*, which is then read by the interface scripts.

**Note:** This interface does not yet employ any symmetry for OpenMX; therefore, all the pairs used must be calculated using jx, or a custom script must be used.

Only version 3.9 of OpenMX has been tested. For more information, please refer to:  
<https://www.openmx-square.org/openmx_man3.9/index.html>  
Specifically for the jx section:  
<https://www.openmx-square.org/openmx_man3.9/node119.html>

## SPR-KKR interface
The structural information for SPR-KKR is read from the *\*.pot* potential file. The exchange constants are written to the *\*_JXC_XCPLTEN_Jij.dat* file, which contains all the exchange constants.

For more information, please refer to the **SPR-KKR** manual:  
<https://www.sprkkr.org/manual/>

## GPAW interface
In the case of GPAW, $J_{ij}$ are output directly in reciprocal space. The exchange constants and a list of points on which they are calculated should be saved in files called *j_\*.npy* and *q_\*.npy*. The structure is read from the output that GPAW writes to the standard output. It should be contained in a file called *\*.txt*.  
&emsp;This interface takes advantage of the fact that $J_{ij}(\boldsymbol{k})=J_{ij}^*(\boldsymbol{-k})$. Therefore, the exchange constants should only be calculated for k-points with a non-negative $z$ coordinate (upper part of the Brillouin zone).

**Note:** Here, it is hard to identify for which sites in the unit cell the exchange constants were calculated. So, in case there are atoms that were not included in the calculation of the exchange constants, it is necessary to exclude them explicitly using the "--exclude_elements" keyword.

For more information on the calculation of the exchange constants using **GPAW**, please visit:  
<https://gpaw.readthedocs.io/tutorialsexercises/magnetic/mft/mft.html>

## Minimum input
To perform the calculation, it is enough to know the lattice vectors, the positions of the magnetic atoms, their magnetic moments and $n$ nearest interactions. That information can be put into a text file as follows:
```txt
lattice
    4.180022      2.090011      2.090011
    2.090011      4.180022      2.090011
    2.090011      2.090011      4.180022

atoms fractional                      # symbol  x y z   mx my mz  (mu_B)
    Ni  0.00  0.00  0.00     0  0   2
    Ni  0.50  0.50  0.50     0  0   -2

exchange                              # distance [A]   J [meV]
    2.9557   0.7
    4.1800   -9.5
```
This file is composed of three segments, each labeled by a keyword. The segment labeled by the keyword "lattice" contains three lines, with one lattice vector per line. The next keyword "atoms" must be followed by "fractional" or "cartesian", depending on which type of coordinates is used in this segment. Each line then contains a label, followed by three coordinates and three components of the magnetic moment. The number of lines matches the number of atoms in the calculation. The last segment starts with "exchange" and contains the distance between the sites followed by the exchange constant (in meV).  
&emsp;This file serves as an input for the *OstravaMag_PrepareInput* script, which transforms it into the *TB2J.pickle* file that can be read by *OstravaMag* as described above. How much the distances between the sites in the input file may differ from the ones found by the code is specified by the "--tol" keyword. If a distance is off, the possible distances found are written to standard output. In case the distances are omitted and only the exchange constants are written to the input, they are assigned to the $n$ lowest shells found by the code.

**Note:** One limitation of this script arises when there are two $J_{ij}$ with different values but the same distance. In this case, one might try breaking this symmetry by varying the lattice vectors/atomic positions slightly.

## Kpath file
To plot the band structure, a definition of the kpath must be given. For OstravaMag, it is written in an external file, whose position and name are given by the "--kpath_name" argument. An example of this file could look like this:
```txt
20
0 0.5 0.5 X
0 0 0 $\Gamma$

0.5 0.5 0.5 L
0 0 0 $\Gamma$
```
This file defines an $X\rightarrow\Gamma|L\rightarrow \Gamma$ path for an fcc lattice with 20 points on each segment. The first line contains the number of points on each segment. Each subsequent pair of lines contains the starting and finishing points of the segment. Each line contains the relative coordinates of a high-symmetry point and a label at the end of the line. Blank lines are allowed anywhere. To use Greek letters in the plots, use labels as indicated in the example.

## QHA  
In this section, we will show what files are needed to run the qha calculations.  
&emsp;For all qha calculations, exchange constants for multiple volumes are needed. Those should be put in multiple folders inside a separate directory, which is set using the "--directory_name_qha" keyword. No other folders may be in this directory. Moreover, a file called *e-v.dat* needs to be included here, which contains a list of volumes and total energies (the same format as in phonopy's qha). The volumes do not have to coincide with the ones used in the calculation of the exchange constants, but only interpolation can be performed; this means that any volumes in the *e-v.dat* file that lie outside the range are cut off.  
&emsp;To be able to add the magnon contributions to the electronic free energy and DFT total energies in the *fe-v.dat* files of phonopy, just put the *fe-v.dat* file into the qha directory. It must be named that way so that OstravaMag can read it. All the files needed for the magnon calculations alone are needed as well.  
&emsp;One important thing to note here is that the temperature grid used in *fe-v.dat* is also used by the code when the phonopy interface is used. Furthermore, these contributions to the free energies are used in OstravaMag's qha. In this way, one may run qha with just the electronic and magnon contributions.

For more information on qha, please refer to phonopy's documentation:  
<https://phonopy.github.io/phonopy/qha.html>