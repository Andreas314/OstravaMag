import numpy as np

class qha_compute:
    
    def __init__(self, num_V, T_inf):
        T_max = T_inf[0]
        num_T = T_inf[1]
        self.T = np.linspace(0, T_max, num_T)
        self.V = np.zeros(shape = (num_V))
        self.F = np.zeros(shape = (num_T, num_V))
        self.counter_F = 0
        self.counter_V = 0
    
    def add_energy(self, energy):
        self.F[:, self.counter_F] = energy
        self.counter_F += 1

    def add_volume(self, volume):
        self.V[self.counter_V] = volume
        self.counter_V += 1

