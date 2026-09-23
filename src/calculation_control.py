
from src.mag_compute import mag_compute
from src.result_writer import result_writer
from src.qha_compute import qha_compute
from src.tp_compute import tp_compute
from src.interfaces.interface_phonopy import load_fe_v, write_fe_v
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
        self.tc = tp_compute(arguments.epsilon, arguments.max_it, arguments.tol)
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
            self.tc.print_results()

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
        self.compute_magnon_dos = arguments.compute_magnon_dos
        self.compute_magnon_thermo = arguments.compute_magnon_thermo
        self.compute_magnon_projected_dos = arguments.compute_magnon_projected_dos
        self.compute_tc = arguments.compute_tc
        self.directory = arguments.directory_name_qha
        self.comp = mag_compute(arguments)
        data_type = arguments.type
        self.special = data_type in specials
        self.res_name = self.directory + "/fe-v-OstravaMag.dat"

    def make_directories_list(self):
        self.subdirs = []
        for item in os.listdir(self.directory):
            if item == self.res_name:
                continue
            if os.path.isdir(os.path.join(self.directory, item)):
                self.subdirs.append(os.path.join(self.directory, item))
        return len(self.subdirs)
    
    def read_phonopy_data(self):
        V, T, F = load_fe_v(self.directory + "/fe-v.dat")
        self.qha = qha_compute(len(self.subdirs), V, T, F)
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
       self.qha.add_energy(self.comp.thermo_quant[:, 1])
       V = self.comp.compute_volume()
       self.qha.add_volume(V)

    def run_qha(self):
        self.qha.run_calculation()
        if not os.path.exists(self.directory + self.res_name):
            os.makedirs(self.directory + self.res_name)
        write_fe_v(self.res_name,self.qha.V_pho, self.qha.T, self.qha.F_pho)

