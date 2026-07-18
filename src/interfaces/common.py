import numpy as np

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

