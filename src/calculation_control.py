
from src.mag_compute import mag_compute
from src.result_writer import result_writer
from src.qha_compute import qha_compute
import os

class calculation_control:
    def __init__(self, arguments):
        self.compute_magnon_bands = arguments.compute_magnon_bands
        self.compute_magnon_dos = arguments.compute_magnon_dos
        self.compute_magnon_thermo = arguments.compute_magnon_thermo
        self.compute_magnon_projected_dos = arguments.compute_magnon_projected_dos
        self.comp = mag_compute(arguments)
        self.data_type = arguments.type
        self.comp.set_up()
        self.writer = result_writer(self.comp)

    def read_data(self):
        self.comp.read()
        if self.compute_magnon_bands:
            self.comp.read_path()


    def run_calculations(self):
            if self.data_type == "GPAW":
                self.comp.make_inverse_fourier_transform()
                self.writer.write_J()
            if self.compute_magnon_dos:
                if self.data_type != "GPAW":
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
            if self.compute_magnon_bands:
                self.comp.make_fourier_transform("path")
                self.comp.compute_magnon_energies("path")
                self.writer.write_bands()
 
class calculation_control_qha(calculation_control):
    def __init__(self, arguments):
        self.compute_magnon_bands = False
        self.compute_magnon_dos = True
        self.compute_magnon_thermo = True
        self.data_type = arguments.type
        self.directory = arguments.directory_name_qha
        self.comp = mag_compute(arguments)
        self.T_inf = (arguments.max_T, arguments.num_T)

 
    def make_directories_list(self):
        self.subdirs = []
        for item in os.listdir(self.directory):
            if os.path.isdir(os.path.join(self.directory, item)):
                self.subdirs.append(os.path.join(self.directory, item))
        self.qha = qha_compute(len(self.subdirs), self.T_inf)
        return len(self.subdirs)

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
