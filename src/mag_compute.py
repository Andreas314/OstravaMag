from read_J import read_data, read_kpath
import numpy as np
import matplotlib.pyplot as plt
class mag_compute:
    def __init__(self, arguments):
        
        ### variables from the command line
        self.max_R = arguments.max_R
        self.data_name = arguments.name
        self.data_type = arguments.type
        self.compute_magnon_bands = arguments.compute_magnon_bands
        self.kpath_name = arguments.kpath_name
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
        self.Jq = np.array([]) #Js in reciprocal space
        self.reciprocal = np.array([]) #basis of the reciprocal sapce
        ##

        self.magnon_energies = np.array([]) #magnon bands

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
        self.Jr, self.cell, self.magmoms = read_data(self.data_type, self.data_name, self.max_R)
        if self.compute_magnon_bands:
            num_points, kpoints, self.letters = read_kpath(self.kpath_name)
            self._make_path(num_points, kpoints)

    def _fourier_transform_kpath(self):
        self._compute_reciprocal()
        self.Jq = np.zeros(shape = (len(self.kpath), np.shape(self.Jr)[-1], np.shape(self.Jr)[-1]), dtype = complex)
        for q_ind in range(len(self.kpath)):
            q = self.kpath[q_ind]
            #computet he reciprocal vector for the phase calculation
            q_ = np.dot(q, self.reciprocal)
            for Rx in range(-self.max_R[0], self.max_R[0] + 1):
                for Ry in range(-self.max_R[1], self.max_R[1] + 1):
                    for Rz in range(-self.max_R[2], self.max_R[2] + 1):
                        #computet he real space difference vector for the phase calculation
                        R_ = np.dot(np.array([Rx, Ry, Rz]), self.cell)  
                        #shift the indices so that they start in zero
                        self.Jq[q_ind] += self.Jr[Rx + self.max_R[0], Ry + self.max_R[1], Rz +self.max_R[2]] * np.exp(1j * np.dot(q_, R_))
    
    def _magnon_energies(self):
        q0 = self.gamma
        J0 = self.Jq[q0]
        Na = len(J0)
        Nq = self.path_len
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
                    Bq[iq, a, b] -= 2 * (KroneckerDelta(n[a],-n[b]) * self.Jq[iq,a,b] / np.sqrt(S[a] * S[b]))            
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

    def run_computations(self):
        if self.compute_magnon_bands:
            self._fourier_transform_kpath()
            self._magnon_energies()
            plt.xticks(self.letters_position, self.letters)
            for band in range(self.magnon_energies.shape[1]):
                plt.plot(self.kpath_x, np.abs(self.magnon_energies[:, band]))
            plt.gca().set_ylim(bottom=0)
            plt.savefig("bandstructure.png", dpi=300)
        
