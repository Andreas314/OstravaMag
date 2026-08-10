import numpy as np
import matplotlib.pyplot as plt
class qha_compute:
    
    def __init__(self, num_V, T_inf):
        T_max = T_inf[0]
        num_T = T_inf[1]
        self.T = np.linspace(0, T_max, num_T)
        self.V = np.zeros(shape = (num_V))
        self.F = np.zeros(shape = (num_T, num_V))
        self.counter_F = 0
        self.counter_V = 0
        self.U = np.zeros(shape = (num_V))
        self.F_tot = np.zeros(shape = (num_T, num_V))
        self.V_0 = np.zeros(shape = (num_T)) 
        self.G = np.zeros(shape = (num_T)) 

    def add_energy(self, energy):
        self.F[:, self.counter_F] = energy
        self.counter_F += 1

    def add_volume(self, volume):
        self.V[self.counter_V] = volume
        self.counter_V += 1

    def add_tot_energies(self, energies):
        if len(energies) != len(self.U):
            raise RuntimeError("Wrong number of total energies in the OstravaMag_qha.txt file!")
        self.U = energies

    def run_calculation(self):
        self._sort_arrays()
        self._sum_contributions()
        self._find_minima()
        self._compute_cp()

    def _sort_arrays(self):
        idx = np.argsort(self.V)
        self.V = self.V[idx]
        self.F = self.F[:, idx]

    def _sum_contributions(self):
        self.F_tot = self.F + self.U

    def _find_minima(self):
        V_array = np.linspace(min(self.V), max(self.V), 1000)
        quart = int(len(self.T) / 4)
        ind = np.array([0, 1, 2, 3, 4]) * quart
        ind[-1] -= 1
        self.to_plot = []
        for i, entry in enumerate(self.F_tot):
            p = np.polyfit(self.V, entry, 2)
            minimum = -p[1] / p[0] / 2
            self.V_0[i] = minimum
            polynomial = np.poly1d(p)
            self.G[i] = polynomial(minimum)
            if p[0] < 0:
                print("WARNING: No minimum at ", self.T[i])
            if i in ind:
                pp = np.poly1d(p)
                self.to_plot.append(pp(V_array))
    
    def _compute_cp(self):
        N_avo = 6.02214076e23
        mEv_to_J =  1.60218e-22
        self.Cp = -self.T * np.gradient(np.gradient(self.G, max(self.T / len(self.T))), max(self.T / len(self.T)))
        self.Cp *= N_avo * mEv_to_J

    def write_results(self, directory):
        V_array = np.linspace(min(self.V), max(self.V), 1000)
        quart = int(len(self.T) / 4)
        ind = np.array([0, 1, 2, 3, 4]) * quart
        ind[-1] -= 1
        for j, i in enumerate(ind):
            plt.plot(self.V, self.F_tot[i] / 1000, "o", label = str(int(self.T[i])) + "K")
            plt.plot(V_array, self.to_plot[j] / 1000, "--")

        plt.legend()
        plt.xlabel(r"V ($\AA^3$)")
        plt.ylabel(r"Free energy (eV)")
        plt.savefig(directory + "/qha_fits.png")
        plt.clf()
        
        plt.plot(self.T, self.V_0)
        plt.ylabel(r"V_0 ($\AA^3$)")
        plt.xlabel(r"T (K)")
        plt.savefig(directory + "/qha_V0.png")
        plt.clf()
        
        plt.plot(self.T[:-2], self.Cp[:-2])
        plt.xlabel(r"T (K)")
        plt.ylabel(r"$C_p$ (J$\cdot$K$^{-1}$$\cdot$mol$^{-1}$)")
        plt.savefig(directory + "/qha_Cp.png")

        with open(directory + "/qha_results.dat", "w") as f:
                f.write(f"T (K)\tV_0 (Ang)\tCp (J/K/mol)\n")
                for i in range(len(self.T)):
                    output = f"{self.T[i]:9.6f}\t"
                    output += f"{self.V_0[i]:9.6f}\t"
                    output += f"{self.Cp[i]:9.6f}\t"
                    output += "\n"
                    f.write(output)
