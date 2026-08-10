import numpy as np
from gpaw.response.heisenberg import get_q0_index

kB = 8.617333252e-5 #meV/K

class tp_compute:
    
    def __init__(self, Jq, magmoms, k_grid):
        self.Jq = Jq / 2
        self.magmoms = abs(magmoms[:, 2])
        self.k_grid = k_grid
        q0 = get_q0_index(self.k_grid)
        self.Nq = np.shape(self.k_grid)[0]
        self.J0 = self.Jq[q0]
        self.Na = len(self.J0)

    def compute_MFA(self):
        # D_A = 2*S_A*(S_A+1)/3  -- the prefactor relating <sigma_A> to J
        D = 2.0 * self.magmoms * (self.magmoms + 1.0) / 3.0 * np.eye(self.Na)
        E = np.outer(1 / self.magmoms, 1 / self.magmoms)
    
        # Build the Hermitian matrix H_AB = sqrt(D_A) * J0_AB * sqrt(D_B)
        # so that its largest eigenvalue equals kB * T_C^MFA
        H = D @ self.J0 * E

        eigvals, eigvecs = np.linalg.eigh(H)
    
        idx = np.argmax(eigvals)
        Tc_mfa = eigvals[idx] / kB   # kB in meV/K, J0 assumed in meV -> T in K
    
    
        self.Tc_mfa = Tc_mfa
    
        print(Tc_mfa)
        
         
        
