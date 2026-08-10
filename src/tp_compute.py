import numpy as np
from gpaw.response.heisenberg import get_q0_index

kB = 8.617333252e-5 #meV/K

class tp_compute:
    
    def __init__(self, Jq, magmoms, k_grid):
        self.Jq = Jq / 2
        self.magmoms = abs(magmoms[:, 2])
        E = np.outer(1 / self.magmoms, 1 / self.magmoms)
        self.Jq *= E

        self.k_grid = k_grid
        self.q0 = get_q0_index(self.k_grid)
        self.Nq = np.shape(self.k_grid)[0]
        self.J0 = self.Jq[self.q0]
        self.Na = len(self.J0)
        
        self.epsilon = 0.3 * kB
        self.it_max = 1000
        self.error_min = 1e-8

    def _compute_MFA(self):
        D = 2.0 * self.magmoms * (self.magmoms + 1.0) / 3.0 * np.eye(self.Na)
    
        H = D @ self.J0

        eigvals, eigvecs = np.linalg.eigh(H)
    
        idx = np.argmax(eigvals)
        Tc_mfa = eigvals[idx] / kB   # kB in meV/K, J0 assumed in meV -> T in K
        
        v = eigvecs[:, idx] @ D
        v_mfa = v / np.linalg.norm(v)
    
        self.v_mfa = v_mfa
        self.Tc_mfa = Tc_mfa
    
    def _integral_at_shift(self, s, z):
        """Evaluate BZ-averaged diagonal of [N(q) + z*I]^{-1}, for complex shift z."""
        diag_part = (self.J0 @ s) * np.eye(self.Na)
        if self.Na != 1:
            rest_part = s[None, :, None] * self.Jq
            N = diag_part[None, :, :] - rest_part + z * np.eye(self.Na)[None, :, :]
            integrant = np.diagonal(np.linalg.inv(N), axis1=1, axis2=2)
        else:
            rest_part = s * self.Jq
            N = diag_part - rest_part + z * np.eye(self.Na)
            integrant = 1 / N
        return integrant.sum(axis=0) / self.Nq   # shape (Na,)
    
    def _run_iteration(self, s):
        eps = self.epsilon
        f_ieps      = self._integral_at_shift(s, 1j*eps)
        f_2ieps     = self._integral_at_shift(s, 2j*eps)
        f_eps_ieps  = self._integral_at_shift(s, eps + 1j*eps)
        f_meps_ieps = self._integral_at_shift(s, -eps + 1j*eps)
    
        integral = 4*f_ieps - f_2ieps - f_eps_ieps - f_meps_ieps
        integral = integral.real   # discard small numerical imaginary residue
    
        end_iteration = 2 * (self.magmoms + 1) * self.magmoms / 3 / integral
        return end_iteration

    def calculate(self):
        self._compute_MFA()
        s_old = self.v_mfa
        it = 0
        error = 1
        while error > self.error_min and it < self.it_max:
            it += 1
            s_new = self._run_iteration(s_old)
            norm_last = np.linalg.norm(s_new)
            s_new /= norm_last
            error = np.linalg.norm(s_new - s_old)
            s_old = s_new
        self.error = error
        self.iterations = it
        self.Tc_rpa = norm_last / kB

    def print_results(self):
        print(r"MFA: Tc=", self.Tc_mfa, "K")
        print(r"RPA: Tc=", self.Tc_rpa, "K")
        print("RPA finished with", self.iterations, "iterations and error of",self.error)


        
        

