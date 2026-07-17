import pickle
import numpy as np
import os
from gpaw.response.heisenberg import get_q0_index
from interfaces.interface_TB2J import read_TB2J
from interfaces.interface_OpenMX import read_OpenMX
from interfaces.interface_GPAW import read_GPAW

def read_kpath(filename):
    with open(filename, "r") as f:
        # Read the first line (single number)
        first_line = f.readline().strip()
        num_points = int(first_line)
        letters = []
        # Read the remaining lines (three integers per line)
        ii = 0
        for line_number, line in enumerate(f, start=2):
            parts = line.strip().split()
            if len(parts) == 0:
                continue
            if len(parts) != 4:
                raise RuntimeError(f"Skipping line {line_number}: does not contain exactly 3 numbers and a letter")
                continue

            try:
                a, b, c = map(float, parts[0:3])
                letters.append(parts[3])
                tmp = np.array([a, b, c])
                if ii == 0:
                    kpath = tmp
                else:
                    kpath = np.vstack((kpath, tmp))
                ii += 1
            except ValueError:
                print(f"Skipping line {line_number}: contains non-integer values")
        if ( ii % 2 != 0):
            raise RuntimeError(f"Number of lines has to be even!")
        new_letters = [letters[0]]
        for ii in range(1, len(letters)- 1, 2):
            if (letters[ii] != letters[ii + 1]):
                new_letters.append(letters[ii] + "|" + letters[ii + 1])
            else:
                new_letters.append(letters[ii])
        new_letters.append(letters[-1])
    return num_points, kpath, new_letters

def read_data(data_type, data_name, max_R, structure_name = "", to_exclude = []):
    print("Data format: " + data_type)
    if (data_type == "TB2J"):
        result = read_TB2J(data_name, max_R, to_exclude)
    elif (data_type == "OpenMX"):
        result = read_OpenMX(data_name, max_R, structure_name, to_exclude)
    elif (data_type == "GPAW"):
        result = read_GPAW(data_name, structure_name, to_exclude)
    else:
        raise ValueError("Unknown data type " + data_type +"!")
    return result

