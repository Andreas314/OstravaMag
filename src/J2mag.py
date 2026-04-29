#!/bin/env python3

import numpy as np
from arg_parse import arg_parse
from mag_compute import mag_compute
if __name__ == "__main__":
    arguments = arg_parse()
    handle_comp = mag_compute(arguments) 
    handle_comp.read()
    handle_comp.run_computations()

