import numpy as np

from lid_driven_cavity_python.TDMA_solver import solve_TDMA
from lid_driven_cavity_python.SIMPLE_solver import solve_SIMPLE
# lid driven cavity solver

# part - A

# given data
nx = 5
ny = 5
re_number_top = 1000
re_number_bottom = 1000

if __name__ == "__main__":
     # testing TDMA solve
     # test_a = [[-2.6, 1, 0, 0],
     #           [1, -2.6, 1, 0],
     #           [0, 1, -2.6, 1],
     #           [0, 0, 1, -2.6]]
     # test_b = [-240, 0, 0, -150]

     # a = [[0, 1, 1, 1],
     #      [-2.6, -2.6, -2.6, -2.6],
     #      [1, 1, 1, 0]]

     # x = solve_TDMA(4, a, test_b)
     # print(x)

     # solve part - A
     # given data
     nx = 5
     ny = 5
     re_number_top = 1000
     re_number_bottom = 1000

     # solution
     u_vel, v_vel, pressure, u_res, v_res, p_res = solve_SIMPLE(nx, ny, re_number_top, re_number_bottom)
