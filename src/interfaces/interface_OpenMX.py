import numpy as np
from src.interfaces.common import find_non_mag

def read_OpenMX(data_name, max_R, structure_name, to_exclude):
    moments, mag_atoms = read_moments_OpenMX(structure_name, data_name, to_exclude)
    J_size = len(mag_atoms)
    cell = read_cell_OpenMX(data_name, structure_name)
    R_size_x = max_R[0] * 2 + 1
    R_size_y = max_R[1] * 2 + 1
    R_size_z = max_R[2] * 2 + 1
    max_elements_count = R_size_x * R_size_y * R_size_z * J_size * J_size
    result = np.zeros(shape=(R_size_x, R_size_y, R_size_z, J_size, J_size), dtype=float, order='C')
    seen_elements = []
    switch = False
    elements = 0
    with open(data_name + "/jx.out", "r") as f:
        for line in f:
            if not switch:
                if line == 97*"-" + "\n":
                     switch = True
                continue
            else:
                try:
                    if line.strip():
                        parts = line.strip().split()
                        if (parts[0] == "Elapsed"):
                            break
                        index = []
                        flag = False
                        for ii in range(2,5):
                            if (abs(int(parts[ii])) > max_R[ii - 2]):
                                flag = True
                        if flag:
                            continue
                        for ii in range(5):
                            if ii < 2:
                                index.append(int(parts[ii]))
                            else:
                                index.append(int(parts[ii]) + max_R[ii - 2])
                        index = tuple(index)
                        if (not (index[0] - 1 in mag_atoms)) or (not (index[1] - 1 in mag_atoms)):
                            continue
                        if not index in seen_elements:  
                            elements += 1
                            seen_elements.append(index)
                            J_1 = np.where(mag_atoms == index[0] - 1)[0]
                            J_2 = np.where(mag_atoms == index[1] - 1)[0]
                            result[index[2], index[3], index[4], J_1, J_2] = (2 * float(parts[5]) / 1000 
                                    *np.dot(moments[J_1][0], moments[J_2][0]) / np.linalg.norm(moments[J_1][0]) / np.linalg.norm(moments[J_2][0]))
                            if (index[2] == max_R[0] and index[3] == max_R[1] and index[4] == max_R[2]) and (J_1 == J_2):
                                result[index[2], index[3], index[4], J_1, J_2] = 0
                        else:
                            print("Warning, element repeats at indices " + str(index))
                except ValueError:
                    continue
       # if elements != max_elements_count:
       #     raise RuntimeError("Wrong number of elements in the input file!")
        return result, cell, moments, mag_atoms

def read_cell_OpenMX(data_name, structure_name):
    #Reads structural input from the *.out file
    switch = False
    cell = np.zeros(shape = (3,3), dtype = float)
    print("Reading structure from: " + data_name + "/" + structure_name + ".out")
    with open(data_name + "/" + structure_name + ".out", "r") as f:
        for line in f:
            ##Find the structural block
            if line == "<Atoms.UnitVectors\n":
                switch = True
                ii = 0
                continue
            if line == "Atoms.UnitVectors>\n":
                return cell
            ###
            ###Read the block
            if switch:
                parts = line.strip().split()
                for jj in range(3):
                    cell[ii, jj] = float(parts[jj])
                ii += 1
            ###
                    

def read_moments_OpenMX(structure_name, data_name, to_exclude):
    #Read the magnetic moments and filter the elements
    moments = np.array([0, 0, 0])
    switch = False
    with open(data_name + "/" + structure_name + ".out", "r") as f:
        for line in f:
            ###Find the moment block
            if not switch:
                if ("Total spin" in line):
                    switch = True
                    next(f, None)
                    next(f, None)
                continue
            else:
            ###
                if not line.strip():
                    moments = np.delete(moments, 0, axis = 0)
                    ###Find which indices are used in the J calculation
                    indices = find_indices_OpenMX(data_name, to_exclude)
                    ###Filter the atoms
                    magnetic_atoms = find_non_mag(moments, indices, to_exclude)
                    moments = moments[magnetic_atoms]
                    return moments, magnetic_atoms
            ###Read the block
                else:
                    parts = line.strip().split()
                    moments = np.vstack((moments, np.array([ 0, 0, float(parts[5])])))
            ###

def find_indices_OpenMX(data_name, to_exclude):
    switch = False
    left_indices = set([])
    right_indices = set([])
    atoms = 0
    print("Reading J_ij from: " + data_name + "/jx.out")
    with open(data_name + "/jx.out", "r") as f:
        for line in f:
            parts = line.strip().split()
            if not switch:
                if line.strip():
                    if parts[0] == "atomnum":
                        atoms = int(parts[2])
                if line == 97*"-"+"\n":
                     switch = True
                continue
            else:
                try:
                    if line.strip():
                        if not (int(parts[0]) - 1 in to_exclude or int(parts[1]) - 1 in to_exclude):
                            left_indices.add(int(parts[0]) - 1)
                            right_indices.add(int(parts[1]) - 1)
                    else:
                        break
                except ValueError:
                    continue
    if (left_indices == right_indices):
        indices = np.zeros(atoms, dtype = bool)
        for index in range(atoms):

            if index in right_indices:
                indices[index] = True
            else:
                indices[index] = False
    else:
        raise RuntimeError("Some elements are missing! (i indices are different then j indices in jx the file)")
    return indices
