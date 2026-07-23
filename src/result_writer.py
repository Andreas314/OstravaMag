from src.mag_compute import mag_compute
import matplotlib.pyplot as plt
import numpy as np

class result_writer:
    def __init__(self, comp):
        self.comp = comp
        self.bands_plot = self.comp.dir + "/magnon_bands.png" 
        self.bands_data = self.comp.dir + "/magnon_bands.dat"
        self.tdos_plot = self.comp.dir + "/magnon_tdos.png"
        self.tdos_data = self.comp.dir + "/magnon_tdos.dat"
        self.pdos_plot = self.comp.dir + "/magnon_pdos.png"
        self.pdos_data = self.comp.dir + "/magnon_pdos.dat"
        self.projected_dos_plot = self.comp.dir + "/magnon_projected_dos.png"
        self.projected_dos_data = self.comp.dir + "/magnon_projected_dos.dat"
        self.energy_on_grid = self.comp.dir+"/energy_on_grid.dat"
        self.J_data = self.comp.dir + "/Jr_GPAW.dat"
        self.thermo_plot = self.comp.dir + "/thermo.png"
        self.thermo_data = self.comp.dir + "/thermo.dat"

 
    def write_bands(self):
        plt.xticks(self.comp.letters_position, self.comp.letters)
        plt.ylabel("E (meV)")
        for pos in self.comp.letters_position:
            plt.axvline(x = pos, color = "grey", lw = 0.2)
        for band in range(self.comp.magnon_energies.shape[1]):
            plt.plot(self.comp.kpath_x, np.real(self.comp.magnon_energies[:, band]))
        plt.gca().set_ylim(bottom=0)
        plt.savefig(self.bands_plot, dpi=300)
        plt.clf()
        with open(self.bands_data, "w") as f:
            for band in range(self.comp.magnon_energies.shape[1]):
                for i, k in enumerate(self.comp.kpath_x):
                    if (self.comp.magnon_energies[i, band] > 0 and np.imag(self.comp.magnon_energies[i, band]) < 1e-3):
                        f.write(f"{k:8.5f}" + "\t" + f"{np.real(self.comp.magnon_energies[i, band]):8.4f}"+"\n")
                f.write("\n")


    def write_tdos(self):
        plt.xlabel("E (meV)")
        plt.ylabel("DOS (states/meV)")
        plt.plot(self.comp.dos_E, self.comp.magnon_dos)
        plt.savefig(self.tdos_plot, dpi=300)
        plt.clf()
        with open(self.tdos_data, "w") as f:
            for i, E0 in enumerate(self.comp.dos_E):
                f.write(f"{E0:8.4f}" + "\t" + f"{self.comp.magnon_dos[i]:8.7f}"+"\n")


    def write_energy_grid(self):
       with open(self.energy_on_grid, "w") as f:
                for i, k in enumerate(self.comp.k_grid):
                    number = ""
                    for band in range(self.comp.magnon_energies.shape[1]):
                        if (self.comp.magnon_energies[i, band] > 0):
                            number += f"\t{np.real(self.comp.magnon_energies[i,band]):8.4f}"
                    f.write(f"{k[0]:8.5f}\t{k[1]:8.5f}\t{k[2]:8.5f}" + number +"\n")
                f.write("\n")


    def write_pdos(self):
        plt.xlabel("E (meV)")
        plt.ylabel("DOS (states/meV)")
        for i in range(len(self.comp.partial_magnon_dos)):
            plt.plot(self.comp.dos_E, self.comp.partial_magnon_dos[i])
        plt.savefig(self.pdos_plot, dpi=300)
        plt.clf()
        
        with open(self.pdos_data, "w") as f:
            for i, E0 in enumerate(self.comp.dos_E):
                text = f"{E0:8.4f}" + "\t"
                for branch in range(len(self.comp.partial_magnon_dos)):
                    text += f"{self.comp.partial_magnon_dos[branch][i]:8.7f}\t"
                text += "\n"
                f.write(text)


    def write_projected_dos(self):
      plt.xlabel("E (meV)")
      plt.ylabel("DOS (states/meV)")
      for i in range(len(self.comp.projected_magnon_dos)):
          plt.plot(self.comp.dos_E, self.comp.projected_magnon_dos[i], label = self.comp.mag_elem_map[i])
      plt.legend()
      plt.savefig(self.projected_dos_plot, dpi=300)
      plt.clf()
      
      with open(self.projected_dos_data, "w") as f:
          for i, E0 in enumerate(self.comp.dos_E):
              text = f"{E0:8.4f}" + "\t"
              for atom in range(len(self.comp.projected_magnon_dos)):
                  text += f"{self.comp.projected_magnon_dos[atom][i]:8.7f}\t"
              text += "\n"
              f.write(text)

         

    def write_J(self):
        with open(self.J_data, "w") as f:
             f.write(f"i j Rx Ry Rz J_ij (meV)"+"\n")
             for Rx in range(-self.comp.max_R[0], self.comp.max_R[0] + 1):
                 for Ry in range(-self.comp.max_R[1], self.comp.max_R[1] + 1):
                     for Rz in range(-self.comp.max_R[2], self.comp.max_R[2] + 1):
                         J = self.comp.Jr[Rx + self.comp.max_R[0], Ry + self.comp.max_R[1], Rz + self.comp.max_R[2]]
                         for i in range(np.shape(J)[0]):
                             for j in range(np.shape(J)[0]):
                                 if not((i == j) and (Rx == 0 and Ry == 0 and Rz == 0)):
                                    f.write(f"{i:d} {j:d} {Rx:d} {Ry:d} {Rz:d} {J[i, j] * 1000:6.5f}"+"\n")


    def write_thermo(self):
       with open(self.thermo_data, "w") as f:
               f.write(f"T (K)\tU (meV)\tF (meV)\tS (J/K/mol)\tC (J/K/mol)\n")
               for i in range(len(self.comp.interval_T)):
                   output = f"{self.comp.interval_T[i]:9.6f}\t"
                   for j in range(4):
                       output += f"{self.comp.thermo_quant[i, j]:9.6f}\t"
                   output += "\n"
                   f.write(output)
       fig, axes = plt.subplots(2, 2, figsize=(10, 8))
       fig.suptitle("Magnon Thermodynamic Quantities", fontsize=14)
       
       plots = [
           (0, 0, r"U (meV)",                        "Internal Energy",  0),
           (0, 1, r"F (meV)",                        "Free Energy",      1),
           (1, 0, r"S (J$\cdot$K$^{-1}$$\cdot$mol$^{-1}$)",       "Entropy",          2),
           (1, 1, r"$C_V$ (J$\cdot$K$^{-1}$$\cdot$mol$^{-1}$)",   "Heat Capacity",    3),
       ]
       
       for row, col, ylabel, title, idx in plots:
           ax = axes[row, col]
           ax.plot(self.comp.interval_T, self.comp.thermo_quant[:, idx])
           ax.set_xlabel("T (K)")
           ax.set_ylabel(ylabel)
           ax.set_title(title)
       
       plt.tight_layout()
       plt.savefig(self.thermo_plot, dpi=300)
       plt.clf()
