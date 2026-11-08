import argparse

def arg_parse():
    parser = argparse.ArgumentParser(prog='OstravaMag', usage='%(prog)s [options]')
    parser.add_argument('--max_R', nargs=3, type = int,  help='Maximal R to consider (default: 2,2,2)', default = [2, 2, 2])
    parser.add_argument('--directory_name', type = str,  help='Name of the directory with the input files (default: .)', default = ".")
    parser.add_argument('--root_name', type = str,  help='Root name of the input files (default: none)', default = "")
    parser.add_argument('--type', type = str,  help='The type of the input file (default: TB2J )', default = "TB2J")
    parser.add_argument('--exclude_elements',nargs = '*', type = int,  help='Indices of elements to exclude from the calculation (default: none)', default = [])
    parser.add_argument('--compute_magnon_bands', action = "store_true",  help='Turns on the computation of the magnon band structure' )
    parser.add_argument('--compute_magnon_dos', action = "store_true",  help='Turns on the computation of the magnon density of states' )
    parser.add_argument('--compute_magnon_thermo', action = "store_true",  help='Turns on the computation of the related thermodynamical quantities' )
    parser.add_argument('--compute_magnon_projected_dos', action = "store_true",  help='Turns on the computation of the dos projected on the included atoms' )
    parser.add_argument('--kpath_name', type = str,  help='Name of a file with a kpath for the calculation of the magnon bandstructure (no default)', default = "__UNDEFINED__" )
    parser.add_argument('--k_grid', nargs=3, type = int,  help='Number of kpoints in each direction of the BZ (default: 8,8,8)', default = [8, 8, 8])
    parser.add_argument('--smearing', type = float,  help='Size of gaussian smearing in meV (default: 5 meV)', default = 5)
    parser.add_argument('--max_T',  type = float,  help='Maximal temperature to calculatet the thermodynamical quantities (default: 300 K)', default = 300)
    parser.add_argument('--num_T',  type = int,  help='Number of points in the temperature interval (default: 300)', default = 300)
    parser.add_argument('--n_points', type = int,  help='Number of points in DOS (default: 50)', default = 50)
    #qha
    parser.add_argument('--compute_qha', action = "store_true",  help='Turns on the quasi-harmonic approximation calculations' )
    parser.add_argument('--directory_name_qha', type = str,  help='Name of the directory with the directories for the qha calculations (default: .)', default = ".")
    parser.add_argument('--phonopy', action = "store_true", help = 'When turned on, reads phonopy\'s fe-v.dat files and adds the magnon contributions to free energy (qha must be turned on)')
    parser.add_argument('--n_cells', type = int,  help='Number of unit cells used in the magnon calculation, used to normalise data in qha (default: 1)', default = 1)
    #critical temperature
    parser.add_argument('--compute_tc', action = "store_true",  help='Turns on the computation of the critical temperature' )
    parser.add_argument('--epsilon', type = float,  help='The epsilon parameter inside the RPA calculation (default: 0.1 K)', default = 0.1)
    parser.add_argument('--tol', type = float,  help='The minimal error tolerance of the RPA calculation (default: 1e-8)', default = 1e-8)
    parser.add_argument('--max_it', type = int,  help='The maximum number of iterations of the RPA calculation (default: 100)', default = 100)
    arguments = parser.parse_args()
    return arguments

