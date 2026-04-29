import argparse

def arg_parse():
    parser = argparse.ArgumentParser(prog='J2mag', usage='%(prog)s [options]')
    parser.add_argument('--max_R', nargs=3, type = int,  help='Maximal R to consider (default: 2,2,2)', default = [2, 2, 2])
    parser.add_argument('--name', type = str,  help='Relative path to the input file (default: ./TB2J.pickle)', default = "./TB2J.pickle")
    parser.add_argument('--type', type = str,  help='The type of the input file (default: TB2J )', default = "TB2J")
    parser.add_argument('--compute_magnon_bands', action = "store_true",  help='Turns on computation of magnon band structure' )
    parser.add_argument('--kpath_name', type = str,  help='Name of a file with a kpath for the calculation of the magnon bandstructure (no default)', default = "__UNDEFINED__" )
    arguments = parser.parse_args()
    return arguments

