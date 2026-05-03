from read_J import read_data, read_kpath
import numpy as np
import matplotlib.pyplot as plt
import os

class mag_compute:
    def __init__(self, arguments):
        
        ### variables from the command line
        ##switches
        self.compute_magnon_bands = arguments.compute_magnon_bands
        self.compute_magnon_dos = arguments.compute_magnon_dos
        ##
        ##others
        self.max_R = arguments.max_R
        self.data_name = arguments.name
        self.structure_name = arguments.name_structure
        self.to_exclude = np.array(arguments.exclude_elements, dtype = int)
        self.data_type = arguments.type
        self.kpath_name = arguments.kpath_name
        self.max_E = arguments.max_E
        self.min_E = arguments.min_E
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
        ###
        
        ## reciprocal space variables
        self.Jq = np.array([]) #Js in reciprocal space for kpath
        self.Jq_grid = np.array([]) #Js in reciprocal space for kgrid
        self.reciprocal = np.array([]) #basis of the reciprocal sapce
        self.k_grid = np.array([]) #actual full gdrid
        ##

        self.magnon_energies = np.array([]) #magnon bands
        self.dos_E = np.linspace(self.min_E, self.max_E, self.n_points)
        self.magnon_dos = np.zeros(shape = (len(self.dos_E)), dtype = float) #magnon bands
    
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
            vx = kx / self.k_grid_division[0] * self.reciprocal[0]
            for ky in range(self.k_grid_division[1]):
                vy = ky / self.k_grid_division[1] * self.reciprocal[1]
                for kz in range(self.k_grid_division[2]):
                    vz = kz / self.k_grid_division[2] * self.reciprocal[2]
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
        self.Jr, self.cell, self.magmoms = read_data(self.data_type, self.data_name, self.max_R, self.structure_name, self.to_exclude)
        if self.compute_magnon_bands:
            num_points, kpoints, self.letters = read_kpath(self.kpath_name)
            self._make_path(num_points, kpoints)

    def _fourier_transform(self, _type):
        if _type == "path":
            self.Jq = np.zeros(shape = (len(self.kpath), np.shape(self.Jr)[-1], np.shape(self.Jr)[-1]), dtype = complex)
            k_size = len(self.kpath)
        elif _type == "grid":
            self.Jq_grid = np.zeros(shape = (len(self.k_grid), np.shape(self.Jr)[-1], np.shape(self.Jr)[-1]), dtype = complex)
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
                        if _type == "path":
                            self.Jq[q_ind] += self.Jr[Rx + self.max_R[0], Ry + self.max_R[1], Rz +self.max_R[2]] * np.exp(1j * np.dot(q_, R_))
                        else:
                            self.Jq_grid[q_ind] += self.Jr[Rx + self.max_R[0], Ry + self.max_R[1], Rz +self.max_R[2]] * np.exp(1j * np.dot(q_, R_))
    
    def _magnon_energies(self, _type):
        if _type == "path":
            q0 = self.gamma
            Nq = self.path_len
        if _type == "grid":
            q0 = 0
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
                    if _type == "path":
                        Aq[iq, a, b] += 2 * ((KroneckerDelta(a,b) * sJ0[a] / S[a]) 
                              - (KroneckerDelta(n[a],n[b]) * self.Jq[iq,a,b] / np.sqrt(S[a] * S[b])))
                    else:
                        Aq[iq, a, b] += 2 * ((KroneckerDelta(a,b) * sJ0[a] / S[a]) 
                              - (KroneckerDelta(n[a],n[b]) * self.Jq_grid[iq,a,b] / np.sqrt(S[a] * S[b])))
        Bq = np.zeros((Nq, Na, Na), complex)
        for iq in range(Nq):
            for a in range(Na):
                for b in range(Na):
                    if _type == "path":
                        Bq[iq, a, b] -= 2* (KroneckerDelta(n[a],-n[b]) * self.Jq[iq,a,b] / np.sqrt(S[a] * S[b]))
                    else:
                        Bq[iq, a, b] -= 2* (KroneckerDelta(n[a],-n[b]) * self.Jq_grid[iq,a,b] / np.sqrt(S[a] * S[b]))
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
                        Hq[iq,a,b] += Aq[iq,b - Na,a - Na]
            Hq[iq,:,:] = np.dot(Sz,Hq[iq,:,:])
        E_qn, v_qan = np.linalg.eig(Hq)
        self.magnon_energies = np.zeros((Nq, 2 * Na),complex)
        for i in range(len(E_qn)):
            #self.magnon_energies[i] = E_qn[i][np.real(E_qn[i]) >= -1e-2]
            self.magnon_energies[i] = E_qn[i][:]
        self.magnon_energies = np.sort(self.magnon_energies, axis=1) * 1000
    
    def _compute_dos(self):
        num_kpoints = np.shape(self.magnon_energies)[0]
        E_vals = np.real(self.magnon_energies).flatten()
        E_vals = E_vals[E_vals >= 0]
        for i, E0 in enumerate(self.dos_E):
            self.magnon_dos[i] = np.sum(
                np.exp(-(E_vals  - E0)**2 / (2 * self.smearing**2))
            )
        
        self.magnon_dos /= (self.smearing * np.sqrt(2 * np.pi)) * num_kpoints

    def run_computations(self):
        self.dir = "./J2mag_results_"+self.data_type
        if not os.path.exists(self.dir):
            os.makedirs(self.dir)
        self._compute_reciprocal()
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

    def _write(self, type_write):
        if type_write == "path":
            plt.xticks(self.letters_position, self.letters)
            plt.ylabel("E (meV)")
            for band in range(self.magnon_energies.shape[1]):
                plt.plot(self.kpath_x, np.real(self.magnon_energies[:, band]))
            plt.gca().set_ylim(bottom=0)
            plt.savefig(self.dir+"/bandstructure.png", dpi=300)
            plt.clf()
            with open(self.dir+"/bandstructure.dat", "w") as f:
                for band in range(self.magnon_energies.shape[1]):
                    for i, k in enumerate(self.kpath_x):
                        if (self.magnon_energies[i, band] > 0):
                            f.write(f"{k:8.5f}" + "\t" + f"{np.real(self.magnon_energies[i, band]):8.4f}"+"\n")
                    f.write("\n")


        elif type_write == "grid":
            plt.xlabel("E (meV)")
            plt.ylabel("DOS (states/meV)")
            plt.plot(self.dos_E, self.magnon_dos)
            plt.savefig(self.dir+"/dos.png", dpi=300)
            plt.clf()
            with open(self.dir+"/dos.dat", "w") as f:
                for i, E0 in enumerate(self.dos_E):
                    f.write(f"{E0:8.4f}" + "\t" + f"{self.magnon_dos[i]:8.7f}"+"\n")
        
