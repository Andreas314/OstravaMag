#!/bin/env python3
import numpy as np
import sys

def write_J_KKR(location, structure, output, Jr, max_R, lattice, mag_elem, moments):
    positions = []
    labels = []
    cartesian_switch = False
    scale = np.array([1, 1, 1])
    lattice_constant = 1
    a = 1
    with open(location + "/" + structure + ".pot", "r") as f:
        for line in f:
            if line.lower().startswith("alat"):
                lattice_constant = float(line.strip().split()[1]) * 0.529178 
            if line.lower().startswith("sites"):
                while True:
                    line = next(f)
                    if not line.split() or line.startswith("*"):
                        break
                    else:
                        words = line.strip().split()
                        if words[0].lower() == "cartesian":
                            if words[1] == "T":
                                cartesian_switch = True
                            continue
                        if words[0].lower() == "basscale":
                            for ii in range(3):
                                scale[ii] = float(words[ii + 1])
                            continue
                        if words[0].lower() == "iq":
                            continue
                        if cartesian_switch:
                            a = lattice_constant
                        if int(words[0])-1 in mag_elem:
                            coords = [0] * 3
                            for ii in range(3):
                                coords[ii] = float(words[1 + ii]) * scale[ii] * a
                            positions.append(np.array(coords))
                            labels.append(words[0])
    N = len(positions)
    positions = np.array(positions)
    if not cartesian_switch:
        for i in range(len(positions)):
            positions[i] = np.matmul((lattice), positions[i])
    js = []
    for ii in range(N):
        for jj in range(N):
            for Rx in range(-max_R[0], max_R[0] + 1):
                for Ry in range(-max_R[1], max_R[1] + 1):
                    for Rz in range(-max_R[2], max_R[2] + 1):
                        new_line = [0]*5
                        J_iso = Jr[Rx + max_R[0], Ry + max_R[1], Rz + max_R[2], ii, jj] * 1000
                        R = (Rx, Ry, Rz)
                        distance_cell = np.matmul(np.transpose(lattice), R)
                        new_line[0] = labels[ii]
                        new_line[1] = labels[jj]
                        new_line[2] = J_iso * moments[ii] * moments[jj]
                        new_line[4] = R
                        distance_cell = np.matmul(np.transpose(lattice), R)
                        distance_atom = positions[jj] - positions[ii]
                        distance = (distance_cell + distance_atom)
                        new_line[3] = np.linalg.norm(distance)
                        if new_line[3] != 0:
                            js.append(new_line)
    js.sort(key=lambda x: x[3])
    final_write(output, lattice, N, labels, positions, js)


def write_J_GPAW(location, structure, output, Jr, max_R, lattice, mag_elem):
    positions = []
    labels = []
    with open(location + "/" + structure + ".txt", "r") as f:
        for line in f:
            if line.startswith("Positions:"):
                while True:
                    line = next(f)
                    if not line.split():
                        break
                    else:
                        words = line.strip().split()
                        if int(words[0]) in mag_elem:
                            coords = [0] * 3
                            for ii in range(3):
                                coords[ii] = float(words[2 + ii])
                            positions.append(np.array(coords))
                            labels.append(words[1] + words[0])
    N = len(positions)
    js = []
    for ii in range(N):
        for jj in range(N):
            for Rx in range(-max_R[0], max_R[0] + 1):
                for Ry in range(-max_R[1], max_R[1] + 1):
                    for Rz in range(-max_R[2], max_R[2] + 1):
                        new_line = [0]*5
                        J_iso = Jr[Rx + max_R[0], Ry + max_R[1], Rz + max_R[2], ii, jj] * 1000
                        R = (Rx, Ry, Rz)
                        distance_cell = np.matmul(np.transpose(lattice), R)
                        new_line[0] = labels[ii]
                        new_line[1] = labels[jj]
                        new_line[2] = J_iso / 2
                        new_line[4] = R
                        distance_cell = np.matmul(np.transpose(lattice), R)
                        distance_atom = positions[jj] - positions[ii]
                        distance = (distance_cell + distance_atom)
                        new_line[3] = np.linalg.norm(distance)
                        if new_line[3] != 0:
                            js.append(new_line)
    js.sort(key=lambda x: x[3])
    final_write(output, lattice, N, labels, positions, js)


def write_J_OpenMX(location, structure, output, mag_elem, moments):
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
    units_switch = True
    with open(location + "/" + structure + ".out", "r") as f:
        for line in f:
            word = line.strip().split()
            if word == []:
                continue
            if word[0] == "Atoms.SpeciesAndCoordinates.Unit":
                if word[1].lower() == "ang":
                    units_switch = False
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
            if units_switch:
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
                if i not in mag_elem or j not in mag_elem:
                    continue
                new_line[0] = labels[i]
                new_line[1] = labels[j]
                new_line[2] = float(word[5]) * moments[i] * moments[j]
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
    final_write(output, lattice, bumber, labels, positions, js)

def final_write(output, lattice, bumber, labels, positions, js):
    #go from relative to absolute positions
    with open(output, "w") as f:
        f.write("Lattice:"+"\n")
        f.write(str(lattice)+"\n")
        f.write("Positions:"+"\n")
        for i in range(bumber):
            f.write(str(labels[i]) + " " + str(positions[i])+"\n")
        #go through the js
        f.write("Interactions:"+"\n")
        for line in js:
            f.write(str(line)+"\n")
        return
