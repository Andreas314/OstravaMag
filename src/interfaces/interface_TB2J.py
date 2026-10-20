import pickle
import numpy as np
import os
from src.interfaces.common import find_non_mag

def read_TB2J(data_name, max_R, to_exclude):
    print("Reading data from: " + data_name + "/TB2J.pickle")
    with open(data_name + "/TB2J.pickle", 'rb') as file:
        data = pickle.load(file)
    R_size_x = max_R[0] * 2 + 1
    R_size_y = max_R[1] * 2 + 1
    R_size_z = max_R[2] * 2 + 1
    moments = data["spinat"]
    is_computed = data["index_spin"]
    cell = data["atoms"].cell[:]
    magnetic_atoms = find_non_mag(moments, is_computed, to_exclude)
    moments = moments[magnetic_atoms]
    J_size = len(magnetic_atoms)
    result = np.zeros(shape=(R_size_x, R_size_y, R_size_z, J_size, J_size), dtype=float, order='C')
    warning = False
    for R_x in range(-max_R[0], max_R[0] + 1):
        for R_y in range(-max_R[1], max_R[1] + 1):
            for R_z in range(-max_R[2], max_R[2] + 1):
                vec = (R_x, R_y, R_z)
                if not any(key[0] == vec for key in data['exchange_Jdict'].keys()):
                    warning = True
                for J_1 in range(J_size):
                    for J_2 in range(J_size):
                        J_x = magnetic_atoms[J_1]
                        J_y = magnetic_atoms[J_2]
                        key = ((R_x, R_y, R_z), J_x, J_y)
                        index = get_index(R_x, R_y, R_z, max_R)
                        if key in data['exchange_Jdict'].keys():
                            result[index[0], index[1], index[2], J_1, J_2] = 2 * data['exchange_Jdict'][key] * np.sign(moments[J_1, 2]) * np.sign(moments[J_2, 2])
    if warning:
        print("Warning: Some values were not found, those are replaced by 0.")
    return result, cell, moments, magnetic_atoms

def get_index(R_x, R_y, R_z, max_R):
    index = np.zeros(shape=(3,1), dtype = int)
    index[0] = max_R[0] + R_x
    index[1] = max_R[1] + R_y
    index[2] = max_R[2] + R_z
    return index
