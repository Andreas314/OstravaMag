#!/bin/env python3
import numpy as np
import sys

def write_J_OpenMX(location, structure, output):
    #read the positions of the atoms inside the cell + lattice vectors    
    open_block_0 = "<Atoms.SpeciesAndCoordinates"
    close_block_0 = "Atoms.SpeciesAndCoordinates>"
    counter_0 = 0
    
    open_block_1 = "<Atoms.UnitVectors"
    close_block_1 = "Atoms.UnitVectors>"
    counter_1 = 0
    
    word_atoms = "Atoms.Number"
    
    inside_block = 0
    block_type = 0
    lattice = np.zeros(shape=(3,3))
    with open(location + "/" + structure + ".out", "r") as f:
        for line in f:
            word = line.strip().split()
            if word == []:
                continue
            if not inside_block:
                if word[0] == open_block_0:
                    inside_block = 1
                    continue
                elif word[0] == open_block_1:
                    inside_block = 1
                    block_type = 1
                    continue
                elif word[0] == word_atoms:
                    bumber = int(word[1])
                    positions = np.zeros(shape = (bumber, 3))
                    labels = [""] * bumber
                else:
                    if block_type == 1:
                        break
                    continue
            else:
                if block_type == 0:
                    if word[0] == close_block_0:
                        inside_block = 0
                        continue
                    else:
                        labels[counter_0] = word[1] + word[0]
                        for i in range(3):
                            positions[counter_0][i] = float(word[i + 2])
                        counter_0 += 1
                elif block_type == 1:
                    if word[0] == close_block_1:
                        inside_block = 0
                        continue
                    else:
                        for i in range(3):
                            lattice[counter_1][i] = float(word[i])
                        counter_1 += 1
    for i in range(bumber):
            positions[i] = np.matmul(np.transpose(lattice), positions[i])
    js = []
    open_js = "-------------------------------------------------------------------------------------------------"
    inside_js = 0
    with open(location+"/jx.out", "r") as f:
        for line in f:
            word = line.strip().split()
            if word == []:
                continue
            if word[0] == "Elapsed":
                break
            if not inside_js:
                if word[0] == open_js:
                    inside_js = 1
            else:
                R = np.zeros(shape = 3)
                new_line = [0]*5
                try:
                    i = int(word[0]) - 1
                    j = int(word[1]) -1 
                except ValueError:
                    continue
                new_line[0] = labels[i]
                new_line[1] = labels[j]
                new_line[2] = float(word[5])
                for ii in range(3):
                    R[ii] = int(word[ii + 2])
                new_line[4] = (R[0], R[1], R[2])
                distance_cell = np.matmul(np.transpose(lattice), R)
                distance_atom = positions[j] - positions[i]
                distance = (distance_cell + distance_atom)
                new_line[3] = np.linalg.norm(distance)
                if new_line[3] != 0:
                    js.append(new_line)
    js.sort(key=lambda x: x[3])
    #go from relative to absolute positions
    with open(output, "w") as f:
        f.write("Lattice:"+"\n")
        f.write(str(lattice)+"\n")
        f.write("Positions relative:"+"\n")
        for i in range(bumber):
            f.write(str(labels[i]) + " " + str(positions[i])+"\n")
        f.write("Positions absolute:"+"\n")
        for i in range(bumber):
            f.write(str(labels[i]) + " " + str(positions[i])+"\n")
        #go through the js
        f.write("Interactions:"+"\n")
        for line in js:
            f.write(str(line)+"\n")
        return
