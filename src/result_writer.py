from src.mag_compute import mag_compute
from src.interfaces.J_writers import write_J_OpenMX, write_J_GPAW, write_J_KKR
import matplotlib.pyplot as plt
import numpy as np

class result_writer:
    def __init__(self, comp):
        self.comp = comp
        self.bands_plot = self.comp.dir + "/magnon_bands.png" 
        self.bands_data = self.comp.dir + "/magnon_bands.dat"
        self.bands_data_tics = self.comp.dir + "/magnon_bands_tics.dat"
        self.tdos_plot = self.comp.dir + "/magnon_tdos.png"
        self.tdos_data = self.comp.dir + "/magnon_tdos.dat"
        self.pdos_plot = self.comp.dir + "/magnon_pdos.png"
        self.pdos_data = self.comp.dir + "/magnon_pdos.dat"
        self.projected_dos_plot = self.comp.dir + "/magnon_projected_dos.png"
        self.projected_dos_data = self.comp.dir + "/magnon_projected_dos.dat"
        self.energy_on_grid = self.comp.dir+"/energy_on_grid.dat"
        self.J_data = self.comp.dir + "/Jr.dat"
        self.thermo_plot = self.comp.dir + "/thermo.png"
        self.thermo_data = self.comp.dir + "/thermo.dat"

 
    def write_bands(self):
        plt.xticks(self.comp.letters_position, self.comp.letters)
        plt.ylabel("E (meV)")
        for pos in self.comp.letters_position:
            plt.axvline(x = pos, color = "grey", lw = 0.2)
        for band in range(self.comp.magnon_energies.shape[1]):
            plt.plot(self.comp.kpath_x, np.real(self.comp.magnon_energies[:, band]), "r")
        plt.gca().set_ylim(bottom=0)
        plt.savefig(self.bands_plot, dpi=300)
        plt.clf()
        with open(self.bands_data, "w") as f:
            for band in range(self.comp.magnon_energies.shape[1]):
                for i, k in enumerate(self.comp.kpath_x):
                    if (self.comp.magnon_energies[i, band] > 0 and np.imag(self.comp.magnon_energies[i, band]) < 1e-3):
                        f.write(f"{k:8.5f}" + "\t" + f"{np.real(self.comp.magnon_energies[i, band]):8.4f}"+"\n")
                f.write("\n")
        with open(self.bands_data_tics, "w") as f:
            for i in range(len(self.comp.letters)):
                f.write(self.comp.letters[i] + " " + str(self.comp.letters_position[i]) + "\n")


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
        if self.comp.data_type == "OpenMX":
            write_J_OpenMX(self.comp.data_name, self.comp.structure_name, self.J_data, self.comp.mag_elem_map)
        elif self.comp.data_type == "GPAW":
            write_J_GPAW(self.comp.data_name, self.comp.structure_name, self.J_data, self.comp.Jr, self.comp.max_R, self.comp.cell, self.comp.mag_elem_map)
        elif self.comp.data_type == "KKR":
            write_J_KKR(self.comp.data_name, self.comp.structure_name, self.J_data, self.comp.Jr, self.comp.max_R, self.comp.cell, self.comp.mag_elem_map)
        else:
            pass


    def write_thermo(self):
       N_avo = 6.02214076e23
       mEv_to_J =  1.60218e-22
       with open(self.thermo_data, "w") as f:
               f.write(f"#T (K)\tU (meV)\tF (meV)\tS (J/K/mol)\tC (J/K/mol)\n")
               for i in range(len(self.comp.interval_T)):
                   output = f"{self.comp.interval_T[i]:9.6f}\t"
                   output += f"{self.comp.thermo_quant[i, 1] * N_avo * mEv_to_J / 1000:9.14f}\t"
                   for j in range(2):
                       output += f"{self.comp.thermo_quant[i, j+2]:9.14f}\t"
                   output += f"{self.comp.thermo_quant[i, 0] * N_avo * mEv_to_J / 1000:9.14f}\t"
                   output += "\n"
                   f.write(output)
       fig, axes = plt.subplots(2, 2, figsize=(10, 8))
       fig.suptitle("Magnon Thermodynamic Quantities", fontsize=14)
       
       plots = [
           (0, 0, r"F (kJ$\cdot$mol$^{-1}$)",                        "Free Energy",      1),
           (1, 0, r"S (J$\cdot$K$^{-1}$$\cdot$mol$^{-1}$)",       "Entropy",          2),
           (1, 1, r"$C_V$ (J$\cdot$K$^{-1}$$\cdot$mol$^{-1}$)",   "Heat Capacity",    3),
           (0, 1, r"E (kJ$\cdot$mol$^{-1}$)",                        "Internal Energy",  0),
       ]
       
       for row, col, ylabel, title, idx in plots:
           ax = axes[row, col]
           if idx > 1:
                coeff = N_avo * mEv_to_J
           else:
                coeff = N_avo * mEv_to_J / 1000
           ax.plot(self.comp.interval_T, self.comp.thermo_quant[:, idx] * coeff)
           ax.set_xlabel("T (K)")
           ax.set_ylabel(ylabel)
           ax.set_title(title)
       
       plt.tight_layout()
       plt.savefig(self.thermo_plot, dpi=300)
       plt.clf()
