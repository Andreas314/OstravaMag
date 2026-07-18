
from mag_compute import mag_compute
from result_writer import result_writer
class calculation_control:
    def __init__(self, arguments):
        self.compute_magnon_bands = arguments.compute_magnon_bands
        self.compute_magnon_dos = arguments.compute_magnon_dos
        self.compute_magnon_thermo = arguments.compute_magnon_thermo
        self.data_type = arguments.type
        self.comp = mag_compute(arguments)
        self.comp.set_up()
        self.writer = result_writer(self.comp)

    def read_data(self):
        self.comp.read()
        if self.compute_magnon_bands:
            self.comp.read_path()


    def run_calculations(self):
            if self.data_type == "GPAW":
                self.comp.make_inverse_fourier_transform()
                self.write.write_J()
            if self.compute_magnon_bands:
                self.comp.make_fourier_transform("path")
                self.comp.compute_magnon_energies("path")
                self.writer.write_bands()
            if self.compute_magnon_dos:
                if self.data_type != "GPAW":
                    self.comp.make_grid()
                    self.comp.make_fourier_transform("grid")
                self.comp.compute_magnon_energies("grid")
                self.comp.compute_dos()
                self.writer.write_tdos()
            if self.compute_magnon_thermo:
                self.comp.compute_thermo()
                self.writer.write_thermo()
    
