import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
class qha_compute:
    
    def __init__(self, num_V, V_pho, T_pho, F_pho):
        self.V_pho = V_pho
        self.F_pho = F_pho
        self.T = T_pho
        num_T = len(T_pho)
        self.V = np.zeros(shape = (num_V))
        self.F = np.zeros(shape = (num_T, num_V))
        self.counter_F = 0
        self.counter_V = 0

    def add_energy(self, energy):
        energy = energy - energy[0]
        self.F[:, self.counter_F] = energy / 1000
        self.counter_F += 1

    def add_volume(self, volume):
        self.V[self.counter_V] = volume
        self.counter_V += 1

    def run_calculation(self):
        num_T = len(self.T)
        to_plot = [0, int(num_T / 4), int(num_T / 2), int(3 * num_T / 4), num_T -1]
        self._sort_arrays()
        tol = np.mean(np.diff(self.V)) / 50
        mask = np.logical_and(self.V_pho > (min(self.V) - tol), self.V_pho < (max(self.V) + tol))
        if not all(mask):
            print("WARNING: Some values of volume are being cut off!")
            self.V_pho = self.V_pho[mask]
            self.F_pho = self.F_pho[:, mask]
        for i, T in enumerate(self.T):
            spline = CubicSpline(self.V, self.F[i, :])
            values = spline(self.V_pho)
            self.F_pho[i, :] += values
            if i in to_plot:
                V = np.linspace(min(self.V), max(self.V), 1000)
                plt.plot(self.V, self.F[i, :], "o", label = "_nolegend_")
                plt.plot(V, spline(V),label = f"{self.T[i]:4.2f} K")
        plt.legend()
        plt.xlabel(r"V ($\AA$)")
        plt.ylabel("F (eV)")
        plt.savefig("pokus.png")



    def _sort_arrays(self):
        
        idx = np.argsort(self.V)
        self.V = self.V[idx]
        self.F = self.F[:, idx]

        idx = np.argsort(self.V_pho)
        self.V_pho = self.V_pho[idx]
        self.F_pho = self.F_pho[:, idx]
