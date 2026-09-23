import numpy as np
import matplotlib.pyplot as plt
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
        self.F[:, self.counter_F] = energy
        self.counter_F += 1

    def add_volume(self, volume):
        self.V[self.counter_V] = volume
        self.counter_V += 1

    def run_calculation(self):
        self._sort_arrays()

    def _sort_arrays(self):
        
        idx = np.argsort(self.V)
        self.V = self.V[idx]
        self.F = self.F[:, idx]

        idx = np.argsort(self.V_pho)
        self.V_pho = self.V_pho[idx]
        self.F_pho = self.F_pho[:, idx]
