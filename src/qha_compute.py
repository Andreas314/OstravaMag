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

    def sort_arrays(self):
        idx = np.argsort(self.V)
        self.V = self.V[idx]
        self.F = self.F[:, idx]

    def sum_contributions(self):
        self.F_tot = self.F + self.U

    def find_minima(self):
        V_array = np.linspace(min(self.V), max(self.V), 1000)
        quart = int(len(self.T) / 4)
        ind = np.array([0, 1, 2, 3, 4]) * quart
        ind[-1] -= 1
        to_plot = []
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
                to_plot.append(pp(V_array))
        for j, i in enumerate(ind):
            plt.plot(self.V, (self.F_tot[i] - pp(minimum)) / 1000, "o", label = str(int(self.T[i])) + "K")
            plt.plot(V_array, (to_plot[j] - pp(minimum)) / 1000, "--")

        plt.legend()
        plt.xlabel(r"V ($\AA^3$)")
        plt.ylabel(r"Free energy (eV)")
        plt.savefig("Pokus_1.png")
        plt.clf()
        plt.plot(self.T, self.V_0)
        plt.ylabel(r"V_0 ($\AA^3$)")
        plt.xlabel(r"T (K)")
        plt.savefig("Pokus.png")
        plt.clf()
        derivative = self.T * np.gradient(np.gradient(self.G, max(self.T / len(self.T))), max(self.T / len(self.T)))
        plt.plot(self.T, derivative)
        plt.savefig("HokusPokus.png")

