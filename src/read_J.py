import pickle
import numpy as np
import os
from gpaw.response.heisenberg import get_q0_index
###PART FOR KPATH

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
###

###COMMON AREA

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

def find_non_mag(moments, is_computed, to_exclude):
    indx = []
    ii = 0
    for line in moments:
        size = np.sqrt(sum(i * i for i in line))
        if size > 1e-2 and is_computed[ii] != -1 and (not ii in to_exclude):
            indx.append(ii)
        ii+=1
    print("Magnetic atoms found at indices: " + str(indx))
    return np.array(indx, dtype = int)
###


###GPAW
def read_GPAW(data_name, structure_name, to_exclude):
    jq_name = data_name + "/" + "j_qab_" + structure_name + ".npy"
    q_name = data_name + "/" + "q_qc_" + structure_name + ".npy"
    weights_name = data_name + "/q_weights_" + structure_name + ".txt"
    str_name = data_name + "/" + structure_name + ".txt"
    print("Reading J_ab(q) from: " + jq_name)
    jq = np.load(jq_name)
    print("Reading q from: " + q_name)
    q = np.load(q_name)
    print("Reading structure from: " + str_name)
    cell = np.zeros(shape = (3, 3), dtype = float)
    null = np.array([0, 0, 0, 0], dtype = float)
    with open(str_name) as f:
        for i, line in enumerate(f):
            words = line.strip().split()
            if words[0] == "lattice":
                a, b, c = np.array(words[1:4], dtype = float)
                alpha, beta, gamma = np.radians(np.array(words[4:], dtype = float))
                cell[0] = [a, 0, 0]
                cell[1] = [b*np.cos(gamma), b*np.sin(gamma), 0]
                
                cx = c*np.cos(beta)
                cy = c*(np.cos(alpha) - np.cos(beta)*np.cos(gamma)) / np.sin(gamma)
                cz = np.sqrt(c**2 - cx**2 - cy**2)
                
                cell[2] = [cx, cy, cz]
            if words[0] == "magmoms":
                mag_len = len(words) - 1
                moments = np.zeros(shape = (mag_len, 3), dtype = float)
                for j, word in enumerate(words[1:]):
                    moments[j, -1] = float(word)
    mag_elements = find_non_mag(moments, np.ones(len(moments)), to_exclude)
    if (len(mag_elements) != len(jq[0])):
        raise RuntimeError("There is " + str(len(mag_elements)) + " magnetic elements and " + str(len(jq[0])) + " elements were used in the J_ij calculations!\n Please add/exclude the additional elemenets.\n Make sure the non-magnetic elements were not used in the calculations.")
    moments = moments[mag_elements]
    q_new = -q
    j_new = np.conj(jq)
    q0 = get_q0_index(q)
    q_new = np.delete(q_new, q0, axis = 0)
    j_new = np.delete(j_new, q0, axis = 0)
    to_delete = []
    for i in range(len(q_new)):
        if q_new[i][ 2] == 0:
            to_delete.append(i)
    q_new = np.delete(q_new, to_delete, axis = 0)
    j_new = np.delete(j_new, to_delete, axis = 0)
    jq = np.vstack((jq, j_new))
    q = np.vstack((q, q_new))
    return jq, q, mag_elements, cell, moments
###


###OPENMX

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
                if line.strip():
                    parts = line.strip().split()
                    if (parts[0] == "Elapsed"):
                        break
                    index = []
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
                        #match TB2J convention with signs and no self-interaction + 2*
                        result[index[2], index[3], index[4], J_1, J_2] = 2 * float(parts[5])* np.sign(moments[J_1, 2]) * np.sign(moments[J_2, 2]) / 1000
                        if (index[2] == max_R[0] and index[3] == max_R[1] and index[4] == max_R[2]) and (J_1 == J_2):
                            result[index[2], index[3], index[4], J_1, J_2] = 0
                    else:
                        print("Warning, element repeats at indices " + str(index))
        if elements != max_elements_count:
            raise RuntimeError("Wrong number of elements in the input file!")
        return result, cell, moments

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
                    indices = find_indices_OpenMX(data_name)
                    ###Filter the atoms
                    magnetic_atoms = find_non_mag(moments, indices, to_exclude)
                    moments = moments[magnetic_atoms]
                    return moments, magnetic_atoms
            ###Read the block
                else:
                    parts = line.strip().split()
                    moments = np.vstack((moments, np.array([ 0, 0, float(parts[5])])))
            ###

def find_indices_OpenMX(data_name):
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
                if line.strip():
                    left_indices.add(int(parts[0]) - 1)
                    right_indices.add(int(parts[1]) - 1)
                else:
                    break
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
###


###TB2J

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
###


