#!/bin/env python3

import numpy as np
from arg_parse import arg_parse
from calculation_control import calculation_control
if __name__ == "__main__":
    arguments = arg_parse()
    handle_comp = calculation_control(arguments) 
    handle_comp.read_data()
    handle_comp.run_calculations()

