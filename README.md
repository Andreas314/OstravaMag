All inside the src folder:\
For options run ./J2mag -h\
Examples:\
./J2mag.py --name "../examples/Fe2Ti/FM/TB2J_results/TB2J.pickle" --max_R 2 2 2 --compute_magnon_bands --kpath_name ../examples/Fe2Ti/kpath.txt\
./J2mag.py --name "../examples/Fe2Ti/AFM/TB2J_results/TB2J.pickle" --max_R 2 2 2 --compute_magnon_bands --kpath_name ../examples/Fe2Ti/kpath.txt\
With DOS:
./J2mag.py --name "../examples/Fe2Ti/FM/TB2J_results/TB2J.pickle" --max_R 2 2 2 --compute_magnon_bands --kpath_name ../examples/Fe2Ti/kpath.txt --k_grid 15 15 15 --n_points 100 --compute_magnon_dos --min_E -30 --max_E 350 --smearing 2\
./J2mag.py --name "../examples/Fe2Ti/AFM/TB2J_results/TB2J.pickle" --max_R 2 2 2 --compute_magnon_bands --kpath_name ../examples/Fe2Ti/kpath.txt --k_grid 15 15 15 --n_points 100 --compute_magnon_dos --min_E -30 --max_E 350 --smearing 2
