import numpy as np
from gpaw.response.heisenberg import get_q0_index

kB = 8.617333252e-5 #meV/K

class tp_compute:
    
    def __init__(self, epsilon, max_it, error_min):
        self.epsilon = epsilon * kB
        self.max_it = max_it
        self.error_min = error_min

    def set_exchange(self, Jq, magmoms, k_grid):
        #Equation A3 shows that we do not need to transform the elements of Jq
        self.Jq = Jq / 2

        self.k_grid = k_grid
        self.q0 = get_q0_index(self.k_grid)
        self.Nq = np.shape(self.k_grid)[0]
        self.J0 = self.Jq[self.q0]
        self.Na = len(self.J0)
        
    def _compute_MFA(self):
    
        H = self.J0

        eigvals, eigvecs = np.linalg.eigh(H)
    
        idx = np.argmax(eigvals)
        Tc_mfa = eigvals[idx] / kB * 2 / 3
        
        v = eigvecs[:, idx]
        v_mfa = v / np.linalg.norm(v)
    
        self.v_mfa = v_mfa
        self.Tc_mfa = Tc_mfa
    
    def _integral_at_shift(self, s, z):
        diag_part = (self.J0 @ s) * np.eye(self.Na)
        if self.Na != 1:
            rest_part = s[None, :, None] * self.Jq
            N = diag_part[None, :, :] - rest_part + z * np.eye(self.Na)[None, :, :]
            integrant = np.diagonal(np.linalg.inv(N), axis1=1, axis2=2)
        else:
            rest_part = s * self.Jq
            N = diag_part - rest_part + z * np.eye(self.Na)
            integrant = 1 / N
        return integrant.sum(axis=0) / self.Nq
    
    def _run_iteration(self, s):
        eps = self.epsilon
        f_ieps      = self._integral_at_shift(s, 1j*eps)
        f_2ieps     = self._integral_at_shift(s, 2j*eps)
        f_eps_ieps  = self._integral_at_shift(s, eps + 1j*eps)
        f_meps_ieps = self._integral_at_shift(s, -eps + 1j*eps)
    
        integral = 4*f_ieps - f_2ieps - f_eps_ieps - f_meps_ieps
        integral = integral.real   # discard small numerical imaginary residue
    
        end_iteration =  2 / 3 / integral
        return end_iteration

    def calculate(self):
        self._compute_MFA()
        s_old = self.v_mfa
        it = 0
        error = 1
        norm_last = 0
        while error > self.error_min and it < self.max_it:
            it += 1
            s_new = self._run_iteration(s_old)
            norm_new = np.linalg.norm(s_new)
            s_new /= norm_last
            error = abs(norm_new - norm_last)
            norm_last = norm_new
        self.error = error
        self.iterations = it
        self.Tc_rpa = norm_last / kB

    def print_results(self):
        print(r"MFA: Tc=", self.Tc_mfa, "K")
        print(r"RPA: Tc=", self.Tc_rpa, "K")
        print("RPA finished with", self.iterations, "iterations and error of",self.error)


        
        

