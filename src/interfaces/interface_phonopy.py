import numpy as np
import os

def find_fe_v(directory):
    return "fe-v.dat" in os.listdir(directory);

def load_e_v(directory):
    volumes = []
    energies = []
    with open(directory + "/e-v.dat", "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            parts = line.split()
            volumes.append(float(parts[0]))
            energies.append(float(parts[1]))
    volumes = np.array(volumes)
    energies = np.array(energies)
    return volumes, energies

def load_fe_v(directory):
    temperatures = []
    free_energies = []
    with open(directory + "/fe-v.dat", "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
 
            # Data line: first column is temperature, rest are free energies
            parts = line.split()
            temperatures.append(float(parts[0]))
            free_energies.append([float(x) for x in parts[1:]])
 
    temperatures = np.array(temperatures)
    free_energies = np.array(free_energies)
    return temperatures, free_energies

def write_fe_v(filename, volumes, temperatures, free_energies, fmt="%.8f"):
    volumes = np.asarray(volumes)
    temperatures = np.asarray(temperatures)
    free_energies = np.asarray(free_energies)
    with open(filename, "w") as f:
        f.write("# volume: " + " ".join(fmt % v for v in volumes) + " \n")
        f.write("#    T(K)     Free energies\n")
        for T, row in zip(temperatures, free_energies):
            values = " ".join(fmt % val for val in row)
            f.write(f"{fmt % T}  {values} \n")

