## Tutorial #2-NiO with minimum input 
In this tutorial, we would like to show, how to properly set au a calculation with minimum input from experimental data. We will demonstrate it on an example of NiO. We would like to see magnon spectrum based on exchange constants from Table II of Ref [7]. Additionally, we wanna see, what its Néel temperature based on this data. 
## Band structure calculations
&emsp;The lattice constant used in the text is $a_0\approx 4.17 $ A , spin moment of the Ni ions in NiO is $S=1$ and it is in a collinear anti-ferromagnetic state. We can write the supercell geometry the input file for *OstravaMag_PrepareInput* as follows (we use magnetic moments, not spin moments)
```
lattice
    4.170000      2.085000      2.085000
    2.085000      4.170000      2.085000
    2.085000      2.085000      4.170000

atoms fractional
    Ni  0.00  0.00  0.00     0  0   2
    Ni  0.50  0.50  0.50     0  0   -2

exchange
    2.95   0.7
    4.1700   -9.5
```
Create a directory called NiO_tutorial and save this into a file called *input.txt* inside the directory. We want to plot band structure on the following k-path
```
100
0.5 0.5 0.5 X
0 0 0 G

0 0 0 G
0 0.5 0.5 M
```
Save this into *kpath.txt* in the *NiO_tutorial* directory. This gives us all the necessary input to run the calculation. To create the TB2J input from the plain text file and move it into the input directory, run the following command
```
OstravaMag_PrepareInput NiO_tutorial/input.txt --tol 0.01 
mv TB2J.pickle NiO_tutorial/
```
The "--tol" option says how much the shell distances in the input file may differ from the ones calculated by the program. If we instead set it to 0.001, we get the following warning
```
No neighbour shell at distance 2.95 A (tol 0.001).
 Shell   distances found: [2.9486, 4.17,
  5.1072, 5.8973, 6.5933, 7.2227] ...
```
We see that the first two shell distances are in agreement to the ones we want to use for $J_1$ and $J_2$, so writing them to the input file is obsolete. This however might always **not** be the case.  
&emsp;We now have the actual input file for *OstravaMag*, so the band structure can be calculated as (the TB2J format is set as default)
```
OstravaMag --directory_name NiO_tutorial/ --kpath_name NiO_tutorial/kpath.txt /
 --compute_magnon_bands
```
If everything went right, the result should look something like this

![NiO_bands](./figures/nio_bands.png) 
Which can is comparable to Figure 2 from Ref [8].
## Néel temperature calculations
The critical temperature can be simply calculated running the following command
```
OstravaMag --directory_name NiO_tutorial/ --k_grid 50 50 50 --compute_tc /
--max_it 20 --tol 1e-15 
```
which prints the following output to the standard output
```
Results without the (1+S)/S correction:
MFA: Tc= 440.97 K
RPA: Tc= 296.35 K
RPA finished with 2 iterations and error of 7.632783294297951e-17
Results with the (1+S)/S correction:
MFA: Tc= 881.94 K
RPA: Tc= 592.69 K
RPA finished with 3 iterations and error of 2.42861286636753e-16
```
which of those two results if relevant for us, depends on the kinds of compound we are dealing with. In case of itinerant magnetism (e.g. bcc Fe) the classical result without the correction is relevant. On the other case, if the system has localized spin moments, the second result is preferred, which is the case here.  
&emsp; We can compare the results with the results in Table II of Ref. [9].