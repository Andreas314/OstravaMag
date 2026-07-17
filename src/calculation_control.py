
from mag_compute import mag_compute

class calculation_control:
    def __init__(self, arguments):
        self.compute_magnon_bands = arguments.compute_magnon_bands
        self.compute_magnon_dos = arguments.compute_magnon_dos
        self.compute_magnon_thermo = arguments.compute_magnon_thermo
        self.data_type = arguments.type
        self.comp = mag_compute(arguments)
    

    def read_data(self):
        self.comp.read()
        if self.compute_magnon_bands:
            self.comp.read_path()


    def run_calculations(self):
        if self.data_type != "GPAW":
             if self.compute_magnon_bands:
                 self.comp.make_fourier_transform("path")
                 self.comp.compute_magnon_energies("path")
                 self.comp.write("path")
             if self.compute_magnon_dos:
                 self.comp.make_grid()
                 self.comp.make_fourier_transform("grid")
                 self.comp.compute_magnon_energies("grid")
                 self.comp.compute_dos()
                 self.comp.write("grid")
        else:
            self.comp.make_inverse_fourier_transform()
            self.comp.write("Jr")
            if self.compute_magnon_dos:
                self.comp.make_magnon_energies("grid")
                self.comp.compute_dos()
                self.comp.write("grid")
            if self.compute_magnon_bands:
                self.comp.make_fourier_transform("path")
                self.comp.compute_magnon_energies("path")
                self.comp.write("path")
        if self.compute_magnon_thermo:
            self.comp.compute_thermo()
            self.comp.write("thermo")
    
