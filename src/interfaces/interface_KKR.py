import numpy as np
from src.interfaces.common import find_non_mag

def read_KKR(data_name, max_R, structure_name, to_exclude):
    cell, moments, mag_atoms = read_structure_KKR(structure_name, data_name, to_exclude)
    J_size = len(mag_atoms)
    R_size_x = max_R[0] * 2 + 1
    R_size_y = max_R[1] * 2 + 1
    R_size_z = max_R[2] * 2 + 1
    max_elements_count = R_size_x * R_size_y * R_size_z * J_size * J_size - J_size
    result = np.zeros(shape=(R_size_x, R_size_y, R_size_z, J_size, J_size), dtype=float, order='C')
    switch = False
    elements = 0
    indices = []
    with open(data_name + "/" + structure_name + "_JXC_XCPLTEN_Jij.dat", "r") as f:
        for line in f:
            if not switch:
                if line.split():
                    words = line.strip().split()
                    if words[0] == "IT":
                         switch = True
                    continue
            else:
                line = line.replace("-", " -")
                words = line.strip().split()
                R1 = int(words[4])
                R2 = int(words[5])
                R3 = int(words[6])
                if abs(R1) > max_R[0] or abs(R2) > max_R[1] or abs(R3) > max_R[2]:
                    continue
                R1 += max_R[0]
                R2 += max_R[1]
                R3 += max_R[2]
                i = int(words[1])
                j = int(words[3])
                if (i - 1 not in mag_atoms) or (j - 1 not in mag_atoms):
                    continue
                index = (i, j, R1, R2, R3)
                indices.append(index)
                J_1 = np.where(mag_atoms == index[0] - 1)[0]
                J_2 = np.where(mag_atoms == index[1] - 1)[0]
                Jxx = float(words[11])
                Jyy = float(words[12])
                J = (Jxx + Jyy) / 2
                result[index[2], index[3], index[4], J_1, J_2] = J / 1000 * 2
                elements += 1
        if elements != max_elements_count:
            raise RuntimeError("Wrong number of elements in the input file (wrong cut-off?)!")
        return result, cell, moments, mag_atoms


def read_structure_KKR(structure_name, data_name, to_exclude):
    #Reads structural input from the *.out file
    BOHR =  5.29177E-1
    cell = np.zeros(shape = (3,3), dtype = float)
    print("Reading structure from: " + data_name + "/" + structure_name + ".pot")
    moments = []
    types = []
    with open(data_name + "/" + structure_name + ".pot", "r") as f:
        for line in f:
            words = line.strip().split()
            
            if words[0] == "ALAT":
                const = float(words[1]) * BOHR
                for ii in range(3):
                    next_line = next(f)
                    next_words = next_line.strip().split()
                    for jj in range(1, 4):
                        cell[ii][jj - 1] = float(next_words[jj]) * const
            
            if words[0] == "OCCUPATION":
                next(f)
                while True:
                    next_line = next(f)
                    next_words = next_line.strip().split()
                    if "".join(set(next_words[0])) == "*":
                        break
                    types.append(int(next_words[1]))

            if words[0] == "MOMENTS":
                mom_dic = {}
                moments = []
                while True:
                    try:
                        next_words = next(f).strip().split()
                        mag_type = int(next_words[1])
                        next_line = next(f)
                        next_words = next_line.strip().split()
                        moment = float(next_words[2])
                        mom_dic[mag_type]= [0, 0, moment]
                        next(f)
                    except StopIteration:
                        break
                for t in types:
                    moments.append(mom_dic[t])
                moments = np.array(moments)
    indices = find_indices_KKR(structure_name, data_name, to_exclude)
    mag_atoms = find_non_mag(moments, indices, to_exclude)
    moments = moments[mag_atoms]
    return cell, moments, mag_atoms
    


def find_indices_KKR(structure_name, data_name, to_exclude):
    left_indices = []
    right_indices = []
    with open(data_name + "/" + structure_name + "_JXC_XCPLTEN_Jij.dat", "r") as f:
        for line in f:
            words = line.strip().split()
            try:
                i = int(words[1])
                j = int(words[3])
                if i not in left_indices:
                    left_indices.append(i)
                if j not in right_indices:
                    right_indices.append(j)
            except (ValueError, IndexError):
                continue
    if set(left_indices) != set(right_indices):
        raise RuntimeError("Some elements are missing! (i indices are different then j indices in jx the file)")
    return left_indices

                    
