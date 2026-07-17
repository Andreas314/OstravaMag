import numpy as np
from interfaces.common import find_non_mag
from gpaw.response.heisenberg import get_q0_index

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
