import numpy as np

from lid_driven_cavity_python.TDMA_solver import solve_TDMA
from lid_driven_cavity_python.SIMPLE_solver import solve_SIMPLE

# lid driven cavity solver

# TODO: Create an example that runs the lid driven cavity problem

# example - 1
nx = 5  # number of CVs in x direction
ny = 5  # number of CVs in y direction
re_number_top = 1000  # Reynolds number at the top wall
re_number_bottom = 1000  # Reynolds number at the bottom wall

if __name__ == "__main__":

    # solve example - 1
    u_vel, v_vel, pressure, u_res, v_res, p_res = solve_SIMPLE(
        nx, ny, re_number_top, re_number_bottom
    )
