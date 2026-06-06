from read_J import read_data, read_kpath
import numpy as np
import matplotlib.pyplot as plt
import os
from gpaw.response.heisenberg import get_q0_index
class mag_compute:
    def __init__(self, arguments):
        
        ### variables from the command line
        ##switches
        self.compute_magnon_bands = arguments.compute_magnon_bands
        self.compute_magnon_dos = arguments.compute_magnon_dos
        self.compute_magnon_thermo = arguments.compute_magnon_thermo
        ##
        ##others
        self.max_R = arguments.max_R
        self.data_name = arguments.directory_name
        self.structure_name = arguments.root_name
        self.to_exclude = np.array(arguments.exclude_elements, dtype = int)
        self.data_type = arguments.type
        self.kpath_name = arguments.kpath_name
        self.max_T = arguments.max_T
        self.num_T = arguments.num_T
        self.n_points = arguments.n_points
        self.smearing = arguments.smearing
        self.k_grid_division = arguments.k_grid
        ##
        ###

        ### variables related to the kpath
        self.kpath = np.array([])
        self.letters = []
        self.letters_position = [0]
        self.kpath_x = np.array([])
        self.gamma = 0
        self.path_len = 0
        ###
        
        ### variables loaded from the input file
        self.Jr = np.array([])
        self.cell = np.array([]) #cell vectors
        self.magmoms = np.array([]) #magnetic moments inside the cell (only > 0)
        self.mag_elem_map = np.array([]) #indices of the original magnetic elements from GPAW input
        ###
        
        ## reciprocal space variables
        self.Jq = np.array([]) #Js in reciprocal space for kpath
        self.reciprocal = np.array([]) #basis of the reciprocal sapce
        self.k_grid = np.array([]) #actual full gdrid
        self.weights = np.array([]) #weights of the kpoints in the GPAW calculation
        ##

        self.magnon_energies = np.array([]) #magnon bands
        self.interval_T = np.linspace(0, self.max_T, self.num_T)
        self.thermo_quant = np.zeros(shape = (self.num_T, 4)) #Internal energy, Free energy, Entropy, Heat capacity
        self.K_B = 0.08617 #meV/K
    def _compute_reciprocal(self):
        #function to compute the reciprocal basis
        a1, a2, a3 = self.cell
        volume = np.dot(a1, np.cross(a2, a3))
    
        if np.isclose(volume, 0):
            raise ValueError("Lattice vectors are linearly dependent")
    
        b1 = 2 * np.pi * np.cross(a2, a3) / volume
        b2 = 2 * np.pi * np.cross(a3, a1) / volume
        b3 = 2 * np.pi * np.cross(a1, a2) / volume
        self.reciprocal = np.vstack((b1, b2, b3))
    
    def _make_grid(self):
        #generate gamma centred kmesh
        self.k_grid = np.zeros( shape = (np.prod(np.array(self.k_grid_division)), 3), dtype = float)
        ii = 0
        for kx in range(self.k_grid_division[0]):
            vx = kx / self.k_grid_division[0] * np.array([1,0,0])
            for ky in range(self.k_grid_division[1]):
                vy = ky / self.k_grid_division[1] * np.array([0,1,0])
                for kz in range(self.k_grid_division[2]):
                    vz = kz / self.k_grid_division[2] * np.array([0,0,1])
                    self.k_grid[ii] = vx + vy + vz
                    ii += 1

    def _make_path(self, num_points, kpoints):
        #make the kpath from the input
        #also handles the thing about having the gamma point inside the path
        steps = len(kpoints)
        self.path_len = round((num_points + 1) * steps / 2 )
        self.kpath = np.zeros(shape = (self.path_len, 3), dtype = float)
        last = 0
        for ii in range(0, steps, 2):
            begin = kpoints[ii]
            end = kpoints[ii + 1];
            diff = (end - begin) / num_points
            for jj in range(num_points + 1):
                #handle the lenght scale of the x axis
                self.kpath_x = np.append(self.kpath_x, last)
                last += np.dot(diff, diff)
                #handle the actual path
                index = int(ii / 2 * (num_points+1) + jj)
                self.kpath[index] = begin + diff * jj
            #positions of the ticks
            self.letters_position.append(self.kpath_x[-1])
        #test of gamma point is inside the kpath
        if np.any(np.all(np.isin(self.kpath, np.array([0, 0, 0])), axis = 1)):
            mask = np.all(self.kpath == [0, 0, 0], axis=1)
            self.gamma = np.where(mask)[0][0]
        #if not, add gamma point to the end of the kpath
        else:
            self.gamma = round((num_points + 1) * steps / 2)
            self.kpath = np.vstack((self.kpath, np.array([0,0,0])))

    def read(self):
        #reads from all the files
        if (self.data_type != "GPAW"):
            self.Jr, self.cell, self.magmoms = read_data(self.data_type, self.data_name, self.max_R, self.structure_name, self.to_exclude)
        else:
            self.Jq, self.k_grid, self.mag_elem_map, self.cell, self.magmoms = read_data(self.data_type, self.data_name, self.max_R, self.structure_name, self.to_exclude)
        if self.compute_magnon_bands:
            num_points, kpoints, self.letters = read_kpath(self.kpath_name)
            self._make_path(num_points, kpoints)

    def _fourier_transform(self, _type):
        if _type == "path":
            self.Jq = np.zeros(shape = (len(self.kpath), np.shape(self.Jr)[-1], np.shape(self.Jr)[-1]), dtype = complex)
            k_size = len(self.kpath)
        elif _type == "grid":
            self.Jq = np.zeros(shape = (len(self.k_grid), np.shape(self.Jr)[-1], np.shape(self.Jr)[-1]), dtype = complex)
            k_size = len(self.k_grid)
        else:
            return 0
        for q_ind in range(k_size):
            q = 0
            if _type == "path":
                q = self.kpath[q_ind]
            else:
                q = self.k_grid[q_ind]
            #computet he reciprocal vector for the phase calculation
            q_ = np.dot(q, self.reciprocal)
            for Rx in range(-self.max_R[0], self.max_R[0] + 1):
                for Ry in range(-self.max_R[1], self.max_R[1] + 1):
                    for Rz in range(-self.max_R[2], self.max_R[2] + 1):
                        #computet he real space difference vector for the phase calculation
                        R_ = np.dot(np.array([Rx, Ry, Rz]), self.cell)  
                        #shift the indices so that they start in zero
                        self.Jq[q_ind] += self.Jr[Rx + self.max_R[0], Ry + self.max_R[1], Rz +self.max_R[2]] * np.exp(1j * np.dot(q_, R_))

    def _inverse_fourier_transform(self):
            x_size = 2 * self.max_R[0] + 1
            y_size = 2 * self.max_R[1] + 1
            z_size = 2 * self.max_R[2] + 1
            self.Jr = np.zeros(shape = (x_size, y_size, z_size, np.shape(self.Jq)[-1], np.shape(self.Jq)[-1]), dtype = float)
            for Rx in range(-self.max_R[0], self.max_R[0] + 1):
                for Ry in range(-self.max_R[1], self.max_R[1] + 1):
                    for Rz in range(-self.max_R[2], self.max_R[2] + 1):
                        #computet he real space difference vector for the phase calculation
                        R_ = np.dot(np.array([Rx, Ry, Rz]), self.cell)  
                        #shift the indices so that they start in zero
                        tmp = np.zeros(shape = (np.shape(self.Jq)[-1], np.shape(self.Jq)[-1]), dtype = complex)
                        for q_ind, q in enumerate(self.k_grid):
                             q_ = np.dot(q, self.reciprocal)
                             tmp  +=  (self.Jq[q_ind]* np.exp(-1j * np.dot(q_, R_))) / len(self.k_grid)
                        self.Jr[Rx + self.max_R[0], Ry + self.max_R[1], Rz +self.max_R[2]] = np.real(tmp)
    
    def _magnon_energies(self, _type):
        if _type == "path":
            q0 = self.gamma
            Nq = self.path_len
        if _type == "grid":
            q0 = get_q0_index(self.k_grid)
            Nq = np.shape(self.k_grid)[0]
        J0 = self.Jq[q0]
        Na = len(J0)
        n = np.sign(self.magmoms[:, 2])
        S = np.abs(self.magmoms[:, 2])
    
        def KroneckerDelta(i,j):
            if i == j:
                return 1
            else:
                return 0
    
        Sz = np.zeros((2 * Na, 2 * Na), int)
        for a in range(2 * Na):
            for b in range(2 * Na):
                if a < Na and b < Na:
                    Sz[a,b] = KroneckerDelta(a,b)
                elif a >= Na and b >= Na:
                    Sz[a,b] = - KroneckerDelta(a,b)
        sJ0 = np.zeros((Na), complex)
        for a in range(Na):
            for b in range(Na):
                sJ0[a] += J0[a,b] 
    
        Aq = np.zeros((Nq, Na, Na), complex)
        for iq in range(Nq):
            for a in range(Na):
                for b in range(Na):
                        Aq[iq, a, b] += 2 * ((KroneckerDelta(a,b) * sJ0[a] / S[a]) 
                              - (KroneckerDelta(n[a],n[b]) * self.Jq[iq,a,b] / np.sqrt(S[a] * S[b])))
        Bq = np.zeros((Nq, Na, Na), complex)
        for iq in range(Nq):
            for a in range(Na):
                for b in range(Na):
                        Bq[iq, a, b] -= 2* (KroneckerDelta(n[a],-n[b]) * self.Jq[iq,a,b] / np.sqrt(S[a] * S[b]))
        Hq = np.zeros((Nq, 2 * Na, 2 * Na), complex)
        for iq in range(Nq):
            for a in range(2 * Na):
                for b in range(2 * Na):
                    if a < Na and b < Na:
                        Hq[iq, a, b] += Aq[iq,a,b]
                    elif a < Na and b >= Na:
                        Hq[iq,a,b] += Bq[iq,a,b - Na]
                    elif a >= Na and b < Na:
                        Hq[iq,a,b] += np.conj(Bq[iq,b,a - Na])
                    else:
                        Hq[iq,a,b] += (Aq[iq,b - Na,a - Na])
            Hq[iq,:,:] = np.dot(Sz,Hq[iq,:,:])
        E_qn, v_qan = np.linalg.eig(Hq)
        self.magnon_energies = np.zeros((Nq, 2 * Na),complex)
        for i in range(len(E_qn)):
            #self.magnon_energies[i] = E_qn[i][np.real(E_qn[i]) >= -1e-2]
            self.magnon_energies[i] = E_qn[i][:]
        self.magnon_energies = np.sort(self.magnon_energies, axis=1) * 1000
    
    def _compute_dos(self):
        self.dos_E = np.linspace(-3 * self.smearing, max(np.real(self.magnon_energies[:, -1])) + 3*self.smearing, self.n_points)
        self.magnon_dos = np.zeros(shape = (len(self.dos_E)), dtype = float) #magnon bands
        ###
        J0 = self.Jq[0]
        N = len(J0)
        ###
        #TODO declare this
        self.partial_magnon_dos_total = np.zeros(shape = (N, len(self.dos_E)), dtype = float) #magnon bands
        num_kpoints = np.shape(self.magnon_energies)[0]
        E_vals = np.real(self.magnon_energies).flatten()
        imag = (np.imag(self.magnon_energies).flatten() > 1e-3)
        if any(imag):
            print("WARNING: " + str(sum(imag)) + " entries with non-zero imaginary component out of " + str(len(E_vals)))
        E_vals = E_vals[((E_vals > 1e-2) * (~imag))]
        for i, E0 in enumerate(self.dos_E):
            self.magnon_dos[i] = np.sum(
                np.exp(-(E_vals  - E0)**2 / (2 * self.smearing**2))
            )
        #TODO declare this
        self.magnon_dos /= (self.smearing * np.sqrt(2 * np.pi)) * num_kpoints
        self.dos_normalization = np.trapz(self.magnon_dos* (self.dos_E[1] - self.dos_E[0]))
        self.magnon_dos /= self.dos_normalization
        print(f"Integral of total DOS: {self.dos_normalization}")


    def _compute_thermo(self):
        mask = self.dos_E >= 0
        self.dos_E = self.dos_E[mask]
        self.magnon_dos = self.magnon_dos[mask]
        
        dE = self.dos_E[1] - self.dos_E[0]
        g = self.magnon_dos * self.dos_normalization #do not output per magnetic atom
        E = self.dos_E
        kB = self.K_B
    
        # zero point energy — same for all T
        U_zp = np.trapz(g * E * 0.5, dx=dE)
    
        for i, T in enumerate(self.interval_T):
            if T == 0:
                self.thermo_quant[i, 0] = U_zp
                self.thermo_quant[i, 1] = U_zp
                self.thermo_quant[i, 2] = 0.0
                self.thermo_quant[i, 3] = 0.0
                continue
        
            kBT = kB * T
            x = E / kBT
        
            # clip x to a safe range BEFORE any exp/expm1 calls
            # x > ~700 means n ≈ 0 and all thermal terms vanish anyway
            x_safe = np.clip(x, 1e-10, 500.0)
        
            n = 1.0 / np.expm1(x_safe)  # no overflow since x_safe <= 500
        
            # internal energy
            self.thermo_quant[i, 0] = np.trapz(g * E * (0.5 + n), dx=dE)
        
            # free energy: log(1 - exp(-x)) = log(expm1(x)) - x, stable for all x_safe
            log_term = np.log(np.expm1(x_safe)) - x_safe
            self.thermo_quant[i, 1] = np.trapz(g * (0.5 * E + kBT * log_term), dx=dE)
        
            # entropy
            entropy_term = x_safe / np.expm1(x_safe) - np.log(np.expm1(x_safe)) + x_safe
            self.thermo_quant[i, 2] = kB * np.trapz(g * entropy_term, dx=dE)
        
            # heat capacity: x^2 * e^x / (e^x - 1)^2 = x^2 * n * (n+1)
            # avoids exp entirely — n already computed, and (n+1) is trivial
            cv_term = x_safe**2 * n * (n + 1.0)
            self.thermo_quant[i, 3] = kB * np.trapz(g * cv_term, dx=dE)
        N_avo = 6.02214076e23
        mEv_to_J =  1.60218e-22
        self.thermo_quant *= mEv_to_J * N_avo #change units to J/mol


    def run_computations(self):
        self.dir = "./OstravaMag_results_"+self.data_type
        if not os.path.exists(self.dir):
            os.makedirs(self.dir)
        self._compute_reciprocal()
        if self.data_type != "GPAW":
             if self.compute_magnon_bands:
                 self._fourier_transform("path")
                 self._magnon_energies("path")
                 self._write("path")
             if self.compute_magnon_dos:
                 self._make_grid()
                 self._fourier_transform("grid")
                 self._magnon_energies("grid")
                 self._compute_dos()
                 self._write("grid")
        else:
            self._inverse_fourier_transform()
            self._write("Jr")
            if self.compute_magnon_dos:
                self._magnon_energies("grid")
                self._compute_dos()
                self._write("grid")
            if self.compute_magnon_bands:
                self._fourier_transform("path")
                self._magnon_energies("path")
                self._write("path")
        if self.compute_magnon_thermo:
            self._compute_thermo()
            self._write("thermo")


    def _write(self, type_write):
        if type_write == "path":
            plt.xticks(self.letters_position, self.letters)
            plt.ylabel("E (meV)")
            for pos in self.letters_position:
                plt.axvline(x = pos, color = "grey", lw = 0.2)
            for band in range(self.magnon_energies.shape[1]):
                plt.plot(self.kpath_x, np.real(self.magnon_energies[:, band]))
            plt.gca().set_ylim(bottom=0)
            plt.savefig(self.dir+"/magnon_bands.png", dpi=300)
            plt.clf()
            with open(self.dir+"/magnon_bands.dat", "w") as f:
                for band in range(self.magnon_energies.shape[1]):
                    for i, k in enumerate(self.kpath_x):
                        if (self.magnon_energies[i, band] > 0 and np.imag(self.magnon_energies[i, band]) < 1e-3):
                            f.write(f"{k:8.5f}" + "\t" + f"{np.real(self.magnon_energies[i, band]):8.4f}"+"\n")
                    f.write("\n")


        elif type_write == "grid":
            plt.xlabel("E (meV)")
            plt.ylabel("DOS (states/meV)")
            plt.plot(self.dos_E, self.magnon_dos)
            plt.savefig(self.dir+"/magnon_dos.png", dpi=300)
            plt.clf()
            with open(self.dir+"/magnon_dos.dat", "w") as f:
                for i, E0 in enumerate(self.dos_E):
                    f.write(f"{E0:8.4f}" + "\t" + f"{self.magnon_dos[i]:8.7f}"+"\n")
            with open(self.dir+"/magnons_on_grid.dat", "w") as f:
                    for i, k in enumerate(self.k_grid):
                        number = ""
                        for band in range(self.magnon_energies.shape[1]):
                            if (self.magnon_energies[i, band] > 0):
                                number += f"\t{np.real(self.magnon_energies[i,band]):8.4f}"
                        f.write(f"{k[0]:8.5f}\t{k[1]:8.5f}\t{k[2]:8.5f}" + number +"\n")
                    f.write("\n")
        elif type_write == "Jr":
            with open(self.dir + "/Jr_J2mag_GPAW.out", "w") as f:
                 f.write(f"i j Rx Ry Rz J_ij (meV)"+"\n")
                 for Rx in range(-self.max_R[0], self.max_R[0] + 1):
                     for Ry in range(-self.max_R[1], self.max_R[1] + 1):
                         for Rz in range(-self.max_R[2], self.max_R[2] + 1):
                             J = self.Jr[Rx + self.max_R[0], Ry + self.max_R[1], Rz +self.max_R[2]]
                             for i in range(np.shape(J)[0]):
                                 for j in range(np.shape(J)[0]):
                                     if not((i == j) and (Rx == 0 and Ry == 0 and Rz == 0)):
                                        f.write(f"{i:d} {j:d} {Rx:d} {Ry:d} {Rz:d} {J[i, j] * 1000:6.5f}"+"\n")


        elif type_write == "thermo":
            with open(self.dir+"/thermo.dat", "w") as f:
                    f.write(f"T (K)\tU (J)\tF (J)\tS (J/K/mol)\tC (J/K/mol)\n")
                    for i in range(len(self.interval_T)):
                        output = f"{self.interval_T[i]:9.6f}\t"
                        for j in range(4):
                            output += f"{self.thermo_quant[i, j]:9.6f}\t"
                        output += "\n"
                        f.write(output)
        fig, axes = plt.subplots(2, 2, figsize=(10, 8))
        fig.suptitle("Magnon Thermodynamic Quantities", fontsize=14)
        
        plots = [
            (0, 0, r"U (J$\cdot$mol$^{-1}$)",                        "Internal Energy",  0),
            (0, 1, r"F (J$\cdot$mol$^{-1}$)",                        "Free Energy",      1),
            (1, 0, r"S (J$\cdot$K$^{-1}$$\cdot$mol$^{-1}$)",       "Entropy",          2),
            (1, 1, r"$C_V$ (J$\cdot$K$^{-1}$$\cdot$mol$^{-1}$)",   "Heat Capacity",    3),
        ]
        
        for row, col, ylabel, title, idx in plots:
            ax = axes[row, col]
            ax.plot(self.interval_T, self.thermo_quant[:, idx])
            ax.set_xlabel("T (K)")
            ax.set_ylabel(ylabel)
            ax.set_title(title)
        
        plt.tight_layout()
        plt.savefig(self.dir + "/thermo.png", dpi=300)
        plt.clf()

    
        
