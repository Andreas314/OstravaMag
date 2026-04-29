import pickle
import numpy as np
import os
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

def read_data(data_type, data_name, max_R):
    if (data_type == "TB2J"):
        result = read_TB2J(data_name, max_R)
    else:
        raise ValueError("Unknown data type " + data_type +"!")
    return result

def read_TB2J(data_name, max_R):
    with open(data_name, 'rb') as file:
        data = pickle.load(file)
    R_size_x = max_R[0] * 2 + 1
    R_size_y = max_R[1] * 2 + 1
    R_size_z = max_R[2] * 2 + 1
    moments = data["spinat"]
    is_computed = data["index_spin"]
    cell = data["atoms"].cell[:]
    magnetic_atoms = find_non_mag(moments, is_computed)
    moments = moments[magnetic_atoms]
    J_size = len(magnetic_atoms)
    result = np.zeros(shape=(R_size_x, R_size_y, R_size_z, J_size, J_size), dtype=float, order='C')
    for R_x in range(-max_R[0], max_R[0] + 1):
        for R_y in range(-max_R[1], max_R[1] + 1):
            for R_z in range(-max_R[2], max_R[2] + 1):
                vec = (R_x, R_y, R_z)
                if not any(key[0] == vec for key in data['exchange_Jdict'].keys()):
                    raise IndexError("Couldn't find vector " + str(vec) + " (wrong TB2J cutoff?)!")
                for J_1 in range(J_size):
                    for J_2 in range(J_size):
                        J_x = magnetic_atoms[J_1]
                        J_y = magnetic_atoms[J_2]
                        key = ((R_x, R_y, R_z), J_x, J_y)
                        index = get_index(R_x, R_y, R_z, max_R)
                        if key in data['exchange_Jdict'].keys():
                            result[index[0], index[1], index[2], J_1, J_2] = 2 * data['exchange_Jdict'][key]
    return result, cell, moments

def get_index(R_x, R_y, R_z, max_R):
    index = np.zeros(shape=(3,1), dtype = int)
    index[0] = max_R[0] + R_x
    index[1] = max_R[1] + R_y
    index[2] = max_R[2] + R_z
    return index

def find_non_mag(moments, is_computed):
    indx = []
    ii = 0
    for line in moments:
        size = np.sqrt(sum(i * i for i in line))
        if size > 1e-3 and is_computed[ii] != -1:
            indx.append(ii)
        ii+=1
    print("Magnetic atoms found at indices: " + str(indx))
    return indx

