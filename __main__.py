import numpy as np

from lid_driven_cavity_python.TDMA_solver import solve_TDMA
from lid_driven_cavity_python.SIMPLE_solver import solve_SIMPLE
from lid_driven_cavity_python.vector_plot import plot_velocities
from lid_driven_cavity_python.contour_plot import plot_pressure
from lid_driven_cavity_python.line_plot import plot_convergence

# lid driven cavity solver
def main():
    # example - 1
    nx = 10  # number of CVs in x direction
    ny = 10  # number of CVs in y direction
    re_number_top = 1000  # Reynolds number at the top wall
    re_number_bottom = 1000  # Reynolds number at the bottom wall

    # solve example - 1
    u_vel, v_vel, pressure, u_res, v_res, p_res = solve_SIMPLE(
        nx, ny, re_number_top, re_number_bottom
    )

    # plot results

    # plot 1: velocity vector plot with animation
    plot_velocities(u_vel, v_vel)

    # plot 2: pressure contour plot
    plot_pressure(pressure)

    # plot 3: convergence history plot
    plot_convergence(u_res, v_res, p_res)


if __name__ == "__main__":
    main()