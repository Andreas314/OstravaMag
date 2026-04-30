import argparse

def arg_parse():
    parser = argparse.ArgumentParser(prog='J2mag', usage='%(prog)s [options]')
    parser.add_argument('--max_R', nargs=3, type = int,  help='Maximal R to consider (default: 2,2,2)', default = [2, 2, 2])
    parser.add_argument('--name', type = str,  help='Relative path to the input file (default: ./TB2J.pickle)', default = "./TB2J.pickle")
    parser.add_argument('--type', type = str,  help='The type of the input file (default: TB2J )', default = "TB2J")
    parser.add_argument('--compute_magnon_bands', action = "store_true",  help='Turns on computation of magnon band structure' )
    parser.add_argument('--compute_magnon_dos', action = "store_true",  help='Turns on computation of magnon density of states' )
    parser.add_argument('--kpath_name', type = str,  help='Name of a file with a kpath for the calculation of the magnon bandstructure (no default)', default = "__UNDEFINED__" )
    parser.add_argument('--k_grid', nargs=3, type = int,  help='Number of kpoints in each direcction of the BZ (default: 8,8,8)', default = [8, 8, 8])
    parser.add_argument('--smearing', type = float,  help='Size of gaussian smearing in meV (default: 5 meV)', default = 5)
    parser.add_argument('--max_E',  type = float,  help='Maximal energy of DOS to calculate (default: 600 meV)', default = 600)
    parser.add_argument('--min_E',  type = float,  help='Minimal energy of DOS to calculate (default: 0 meV)', default = 0)
    parser.add_argument('--n_points', type = int,  help='Number of points in DOS (default: 50 meV)', default = 50)
    arguments = parser.parse_args()
    return arguments

