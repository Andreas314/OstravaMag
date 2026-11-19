import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import CubicSpline
class qha_compute:
    
    def __init__(self, num_V, V_pho, T_pho, F_pho, n_cells):
        self.V_pho = V_pho
        self.F_pho = F_pho
        self.T = T_pho
        num_T = len(T_pho)
        self.V = np.zeros(shape = (num_V))
        self.U = np.zeros(shape = (num_V))
        self.F = np.zeros(shape = (num_T, num_V))
        self.counter_F = 0
        self.counter_V = 0
        self.V_0 = np.zeros(shape = (num_T)) 
        self.G = np.zeros(shape = (num_T))
        self.N = -1
        self.n_cells = n_cells

    def add_energy(self, energy):
        energy = energy - energy[0]
        self.F[:, self.counter_F] = energy / 1000
        self.counter_F += 1

    def add_volume(self, volume):
        self.V[self.counter_V] = volume / self.n_cells
        self.counter_V += 1

    def add_tot_energies(self, energies):
        self.U = np.array(energies)

    def run_interface(self):
        self._sort_arrays()
        self._sort_phonopy_arrays()
        num_T = len(self.T)
        to_plot = [0, int(num_T / 4), int(num_T / 2), int(3 * num_T / 4), num_T -1]
        tol = np.mean(np.diff(self.V)) / 50
        mask = np.logical_and(self.V_pho > (min(self.V) - tol), self.V_pho < (max(self.V) + tol))
        if not all(mask):
            print("WARNING: Some values of volume are being cut off!")
            self.V_pho = self.V_pho[mask]
            self.F_pho = self.F_pho[:, mask]
        for i, T in enumerate(self.T):
            spline = CubicSpline(self.V, self.F[i, :])
            values = spline(self.V_pho) * self.N / self.n_cells
            self.F_pho[i, :] += values
    
    def run_calculation(self, phonopy):
        self._sort_arrays()
        self._sum_contributions(phonopy)
        self._find_minima()
        self._compute_cp()

    def _sum_contributions(self, phonopy):
        if phonopy:
            self.F_tot = self.F_pho
        else:
            tol = np.mean(np.diff(self.V)) / 50
            mask = np.logical_and(self.V_pho > (min(self.V) - tol), self.V_pho < (max(self.V) + tol))
            if not all(mask):
                print("WARNING: Some values of volume are being cut off!")
                self.V_pho = self.V_pho[mask]
                self.U = self.U[mask]
            self.F_tot = np.zeros(shape = (len(self.T), len(self.V_pho)))
            for i, T in enumerate(self.T):
                spline = CubicSpline(self.V, self.F[i, :])
                values = spline(self.V_pho) * self.N / self.n_cells
                self.F_tot[i, :] = self.U + values
        self.V = self.V_pho


    def _find_minima(self):
        V_array = np.linspace(min(self.V), max(self.V), 1000)
        V_lo, V_hi = min(self.V), max(self.V)
        quart = int(len(self.T) / 4)
        ind = np.array([0, 1, 2, 3, 4]) * quart
        ind[-1] -= 1
        self.to_plot = []
    
        # extra storage for the maxima (optional, but useful for diagnostics)
        self.V_max = np.full(len(self.F_tot), np.nan)
        self.G_max = np.full(len(self.F_tot), np.nan)
    
        for i, entry in enumerate(self.F_tot):
            p = np.polyfit(self.V, entry, 3)
            poly = np.poly1d(p)
            dp = np.polyder(p)          # quadratic
            d2p = np.polyder(dp)        # linear
    
            # all roots of the derivative, keep only the real ones
            roots = np.roots(dp)
            real_roots = roots[np.abs(roots.imag) < 1e-10].real
    
            minima, maxima = [], []
            for r in real_roots:
                curvature = np.polyval(d2p, r)
                if curvature > 0:
                    minima.append(r)
                elif curvature < 0:
                    maxima.append(r)
                # curvature == 0 -> inflection point, neither
    
            # prefer stationary points inside the sampled volume range
            minima_in = [r for r in minima if V_lo <= r <= V_hi]
            maxima_in = [r for r in maxima if V_lo <= r <= V_hi]
    
            if minima_in:
                # a cubic has at most one local min, but stay general:
                # pick the one with the lowest free energy
                best = min(minima_in, key=poly)
                self.V_0[i] = best
                self.G[i] = poly(best)
            elif minima:
                print(f"WARNING: minimum at T = {self.T[i]} lies outside the "
                      f"volume range ({minima[0]:.4f} not in [{V_lo:.4f}, {V_hi:.4f}])")
                self.V_0[i] = minima[0]
                self.G[i] = poly(minima[0])
            else:
                print(f"WARNING: No minimum at T = {self.T[i]}")
                self.V_0[i] = np.nan
                self.G[i] = np.nan
    
            if i in ind:
                self.to_plot.append(poly(V_array))

    def _compute_cp(self):
        N_avo = 6.02214076e23
        Ev_to_J =  1.60218e-19
        self.Cp = -self.T * np.gradient(np.gradient(self.G, max(self.T / len(self.T))), max(self.T / len(self.T)))
        self.Cp *= N_avo * Ev_to_J


    def _sort_arrays(self):
        idx = np.argsort(self.V)
        self.V = self.V[idx]
        self.F = self.F[:, idx]


    def _sort_phonopy_arrays(self):
        idx = np.argsort(self.V_pho)
        self.V_pho = self.V_pho[idx]
        self.F_pho = self.F_pho[:, idx]

    def write_results(self, directory):
        if not os.path.exists(directory):
            os.makedirs(directory)
        V_array = np.linspace(min(self.V), max(self.V), 1000)
        quart = int(len(self.T) / 4)
        ind = np.array([0, 1, 2, 3, 4]) * quart
        ind[-1] -= 1
        for j, i in enumerate(ind):
            plt.plot(self.V, self.F_tot[i], "o", label = str(int(self.T[i])) + "K")
            plt.plot(V_array, self.to_plot[j], "--")
    
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
        plt.clf()
        
        plt.plot(self.T, self.G)
        plt.xlabel(r"T (K)")
        plt.ylabel(r"G (eV)")
        plt.savefig(directory + "/qha_G.png")
    
        with open(directory + "/qha_results.dat", "w") as f:
                f.write(f"#T (K)\tV_0 (Ang^3)\tCp (J/K/mol) G(eV)\n")
                for i in range(len(self.T)):
                    output = f"{self.T[i]:9.6f}\t"
                    output += f"{self.V_0[i]:9.6f}\t"
                    output += f"{self.Cp[i]:9.6f}\t"
                    output += f"{self.G[i]:9.6f}\t"
                    output += "\n"
                    f.write(output)
        idx = np.argsort(self.V_pho)
