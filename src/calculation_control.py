
from src.mag_compute import mag_compute
from src.result_writer import result_writer
from src.qha_compute import qha_compute
from src.tc_compute import tc_compute
from src.interfaces.interface_phonopy import load_fe_v, write_fe_v, find_fe_v, load_e_v
import os
import numpy as np

specials = ["GPAW"]

class calculation_control:
    def __init__(self, arguments):
        self.compute_magnon_bands = arguments.compute_magnon_bands
        self.compute_magnon_dos = arguments.compute_magnon_dos
        self.compute_magnon_thermo = arguments.compute_magnon_thermo
        self.compute_magnon_projected_dos = arguments.compute_magnon_projected_dos
        self.compute_tc = arguments.compute_tc
        self.comp = mag_compute(arguments)
        self.comp.set_up()
        self.writer = result_writer(self.comp)
        self.tc = tc_compute(arguments.epsilon, arguments.max_it, arguments.tol, arguments.mixing)
        data_type = arguments.type
        self.special = data_type in specials
    
    def read_data(self):
        self.comp.read()
        if self.compute_magnon_bands:
            self.comp.read_path()
    

    def _run_tc_calc(self):
        if self.compute_tc:
            if not self.compute_magnon_dos and not self.special:
                self.comp.make_grid()
            self.comp.make_fourier_transform("grid")
            self.tc.set_exchange(self.comp.Jq, self.comp.magmoms, self.comp.k_grid)
            self.tc.calculate()

    def _run_dos_calc(self):
        if self.compute_magnon_dos:
            if not self.special:
                self.comp.make_grid()
            self.comp.make_fourier_transform("grid")
            self.comp.compute_magnon_energies("grid")
            self.comp.compute_dos()
            self.writer.write_tdos()
            self.writer.write_pdos()
            if self.compute_magnon_projected_dos:
                self.comp.compute_projected_dos()
                self.writer.write_projected_dos()
            if self.compute_magnon_thermo:
                self.comp.compute_thermo()
                self.writer.write_thermo()

    def _run_ift(self):
        if self.comp.data_type == "GPAW":
            self.comp.make_inverse_fourier_transform()
        self.writer.write_J()

    def _run_band_calc(self):
        if self.compute_magnon_bands:
            self.comp.make_fourier_transform("path")
            self.comp.compute_magnon_energies("path")
            self.writer.write_bands()

    def run_calculations(self):
            self._run_ift()
            self._run_dos_calc()
            self._run_tc_calc()
            self._run_band_calc()
 
class calculation_control_qha(calculation_control):
    def __init__(self, arguments):
        self.compute_magnon_bands = arguments.compute_magnon_bands
        self.compute_magnon_dos = True
        self.compute_magnon_thermo = True
        self.compute_magnon_projected_dos = arguments.compute_magnon_projected_dos
        self.phonopy = False
        self.compute_tc = arguments.compute_tc
        self.tc = tc_compute(arguments.epsilon, arguments.max_it, arguments.tol, arguments.mixing)
        self.directory = arguments.directory_name_qha
        self.comp = mag_compute(arguments)
        self.n_cells = arguments.n_cells
        data_type = arguments.type
        self.special = data_type in specials
        self.res_name_phonopy = self.directory + "/fe-v-OstravaMag.dat"
        self.res_name = "OstravaMag_results_qha"
        self.T = np.linspace(0, arguments.max_T, arguments.num_T) 

    def make_directories_list(self):
        self.subdirs = []
        for item in os.listdir(self.directory):
            if item == self.res_name:
                continue
            if os.path.isdir(os.path.join(self.directory, item)):
                self.subdirs.append(os.path.join(self.directory, item))
        return len(self.subdirs)
    
    def read_files(self):
        self._read_phonopy_data()
    
    def _read_phonopy_data(self):
        V, U = load_e_v(self.directory)
        self.phonopy = find_fe_v(self.directory)
        if self.phonopy:
            print("fe-v.dat found, reading temperature field and energies from it")
            T, F = load_fe_v(self.directory)
            if len(V) != len(F[0, :]):
                raise RuntimeError("Different number of entries in fe-v.dat and e-v.dat!")
        else:
            print("fe-v.dat not found, reading temperature field from arguments and energies from e-v.dat")
            T = self.T
            F = []
        self.qha = qha_compute(len(self.subdirs), V, T, F, self.n_cells)
        self.qha.add_tot_energies(U)
        self.comp.max_T = max(T)
        self.comp.num_T = len(T)
        self.comp.interval_T = T
        self.comp.thermo_quant = np.zeros(shape = (len(T), 4))
        
    def advance(self):
       current_dir = self.subdirs.pop()
       self.comp.data_name = current_dir
       self.comp.set_up()
       self.writer = result_writer(self.comp)
       return current_dir
    
    def add_entries(self):
        if self.qha.N == -1:
            self.qha.N = len(self.comp.magmoms)
        elif self.qha.N != len(self.comp.magmoms):
            raise RuntimeError("Different number of magnetic elements in the files!")
        self.qha.add_energy(self.comp.thermo_quant[:, 1])
        self.qha.add_volume(self.comp.compute_volume())
    
    def run_qha(self):
        if self.phonopy:
            self.qha.run_interface()
            write_fe_v(self.res_name_phonopy,self.qha.V_pho, self.qha.T, self.qha.F_pho)
        self.qha.run_calculation(self.phonopy)
        self.qha.write_results(self.directory + "/" + self.res_name)

