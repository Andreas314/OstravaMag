#!/bin/env python3

import numpy as np
from src.arg_parse import arg_parse
from src.calculation_control import calculation_control, calculation_control_qha
def OstravaMag():
    arguments = arg_parse()
    if  not arguments.compute_qha:
        handle_comp = calculation_control(arguments) 
        handle_comp.read_data()
        handle_comp.run_calculations()
    else:
        handle_comp = calculation_control_qha(arguments)
        num = handle_comp.make_directories_list()
        for i in range(num):
            current = handle_comp.advance()
            handle_comp.read_data()
            handle_comp.run_calculations()
            handle_comp.add_entries()
        print(handle_comp.qha.F)

if __name__ == "__main__":
    OstravaMag()
