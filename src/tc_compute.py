import numpy as np
from gpaw.response.heisenberg import get_q0_index

kB = 8.617333252e-5 #meV/K

class tc_compute:
    
    def __init__(self, epsilon, max_it, error_min, mixing):
        self.epsilon = epsilon * kB
        self.max_it = max_it
        self.error_min = error_min
        self.mixing = mixing

    def set_exchange(self, Jq, magmoms, k_grid):
        #Equation A3 shows that we do not need to transform the elements of Jq
        self.magmoms=magmoms
        S = np.linalg.norm(magmoms, axis = 1)
        unit_vectors = np.dot(np.diag(1 / S), magmoms)
        self.Jq = Jq / 2 *  np.inner(unit_vectors,  unit_vectors)
        self.k_grid = k_grid
        self.q0 = get_q0_index(self.k_grid)
        self.Nq = np.shape(self.k_grid)[0]
        self.J0 = self.Jq[self.q0]
        self.Na = len(self.J0)
        self.Tc_mfa = 0
        self.Tc_rpa = 0
    def _compute_MFA(self, apply_correction):
        S = np.linalg.norm(self.magmoms, axis=1) / 2
        correction = (S + 1) / S
        D = np.sqrt(correction)                     # 1D, shape (Na,)
        if apply_correction:
            H = D[:, None] * self.J0 * D[None, :]        # (Na, Na)
        else:
            H = self.J0
        eigvals, eigvecs = np.linalg.eigh(H)
        
        idx = np.argmax(eigvals)
        self.Tc_mfa = eigvals[idx] / kB * 2 / 3
        if apply_correction:
            v = correction * eigvecs[:, idx]                      # back-transform to the eigenvector of diag(c) @ J0
            v = (v)                                # signs are already inside J
        else:
            v = eigvecs[:, idx]                      # back-transform to the eigenvector of diag(c) @ J0
        self.v_mfa = v / np.linalg.norm(v)
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
    
    def _run_iteration(self, s,  apply_correction):
        eps = self.epsilon
        f_ieps      = self._integral_at_shift(s, 1j*eps)
        f_2ieps     = self._integral_at_shift(s, 2j*eps)
        f_eps_ieps  = self._integral_at_shift(s, eps + 1j*eps)
        f_meps_ieps = self._integral_at_shift(s, -eps + 1j*eps)
    
        integral = 4*f_ieps - f_2ieps - f_eps_ieps - f_meps_ieps
        integral = integral.real   # discard small numerical imaginary residue
        S = np.linalg.norm(self.magmoms, axis = 1) / 2
        correction = (S + 1) / S
        end_iteration =  2 / 3 / integral
        if apply_correction:
            end_iteration = np.dot(end_iteration, np.diag(correction))
        return end_iteration

    def calculate(self):
        for ii in [False, True]:
            self._compute_MFA(ii)
            s_old = self.v_mfa
            it = 0
            error = 1
            norm_last = 0
            mix = self.mixing
            while error > self.error_min and it < self.max_it:
                it += 1
                s_new = self._run_iteration(s_old, ii)
                norm_new = np.linalg.norm(s_new)
                s_trial = s_new / norm_new
                s_old = (1 - mix) * s_old + mix * s_trial
                s_old /= np.linalg.norm(s_old)
                error = abs(norm_new - norm_last) / kB
                norm_last = norm_new
            self.error = error
            self.iterations = it
            self.Tc_rpa = norm_last / kB
            if ii:
                print("Results with the (1+S)/S correction:")
            else:
                print("Results without the (1+S)/S correction:")
            print(f"MFA: Tc= {self.Tc_mfa:4.2f} K")
            print(f"RPA: Tc= {self.Tc_rpa:4.2f} K")
            print("RPA finished with", self.iterations, "iterations and error of",self.error)


        
        

