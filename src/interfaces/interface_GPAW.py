import numpy as np
from interfaces.common import find_non_mag
from gpaw.response.heisenberg import get_q0_index

def read_GPAW(data_name, structure_name, to_exclude):
    jq_name = data_name + "/" + "j_" + structure_name + ".npy"
    q_name = data_name + "/" + "q_" + structure_name + ".npy"
    str_name = data_name + "/" + structure_name + ".txt"
    print("Reading J_ab(q) from: " + jq_name)
    jq = np.load(jq_name)
    print("Reading q from: " + q_name)
    q = np.load(q_name)
    print("Reading structure from: " + str_name)
    cell = np.zeros(shape = (3, 3), dtype = float)
    null = np.array([0, 0, 0, 0], dtype = float)
    num_at = 0
    with open(str_name, "r+") as f:
        for line in f:
            words = line.strip().split()
            if line.startswith("Number of atoms:"):
                num_at = int(words[3])
            if line.strip() == "Unit cell:":
                next(f)
                for ii in range(3):
                    next_line = next(f)
                    next_words = next_line.strip().split()
                    for jj in range(3):
                        cell[ii][jj] = float(next_words[3 + jj])
            if line.strip() == "Local magnetic moments:":
                moments = np.zeros(shape = (num_at, 3), dtype = float)
                for jj in range(num_at):
                    next_line = next(f)
                    next_words = next_line.strip().split()
                    moments[jj, -1] = float(next_words[5][:-1])
    mag_elements = find_non_mag(moments, np.ones(len(moments)), to_exclude)
    if (len(mag_elements) != len(jq[0])):
        raise RuntimeError("There is " + str(len(mag_elements)) + " magnetic elements and " + str(len(jq[0])) + " elements were used in the J_ij calculations!\n Please add/exclude the additional elemenets.\n Make sure the non-magnetic elements were not used in the calculations.")
    moments = moments[mag_elements]
    q_new = -q
    j_new = np.conj(jq)
    to_delete = []
    for qq in range(len(q_new)): #make sure there are no duplicates
        if q_new[qq, 2] == 0 or q_new[qq, 2] == -0.5:
            to_delete.append(qq)
    q_new = np.delete(q_new, to_delete, axis = 0)
    j_new = np.delete(j_new, to_delete, axis = 0)
    jq = np.vstack((jq, j_new))
    q = np.vstack((q, q_new))
    return jq, q, mag_elements, cell, moments
