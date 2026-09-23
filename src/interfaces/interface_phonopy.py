import numpy as np

def load_fe_v(filename):
    volumes = None
    temperatures = []
    free_energies = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Read the volume header
            if line.startswith("#"):
                if "volume" in line.lower():
                    after_colon = line.split(":", 1)[1] if ":" in line else line.lstrip("#")
                    volumes = np.array([float(x) for x in after_colon.split()])
                continue
 
            # Data line: first column is temperature, rest are free energies
            parts = line.split()
            temperatures.append(float(parts[0]))
            free_energies.append([float(x) for x in parts[1:]])
 
    if volumes is None:
        raise ValueError(f"No '# volume:' header line found in {filename}")
 
    temperatures = np.array(temperatures)
    free_energies = np.array(free_energies)
 
    if free_energies.shape[1] != volumes.shape[0]:
        raise ValueError(
            f"Number of volume columns in data ({free_energies.shape[1]}) "
            f"does not match number of volumes in header ({volumes.shape[0]})"
        )
    return volumes, temperatures, free_energies

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

