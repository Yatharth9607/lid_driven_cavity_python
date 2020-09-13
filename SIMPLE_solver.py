# SIMPLE solver for 2D
import numpy as np

import lid_driven_cavity_python.properties as properties
from lid_driven_cavity_python.geometry_function import geometry_function
from lid_driven_cavity_python.ics_function import ics_function
from lid_driven_cavity_python.u_vel_coeff import u_vel_coeff
from lid_driven_cavity_python.v_vel_coeff import v_vel_coeff
from lid_driven_cavity_python.u_vel_solve import u_vel_solve
from lid_driven_cavity_python.v_vel_solve import v_vel_solve
from lid_driven_cavity_python.p_corr_coeff import p_corr_coeff
from lid_driven_cavity_python.p_correction import p_correction
from lid_driven_cavity_python.u_vel_correction import u_vel_correction
from lid_driven_cavity_python.v_vel_correction import v_vel_correction
from lid_driven_cavity_python.convergence import convergence


def solve_SIMPLE(nx, ny, re_number_top, re_number_bottom):
    # call geometry function
    ipcv, jpcv, dx, dy, l, _ = geometry_function(nx, ny)
    ipmax = ipcv + 2
    jpmax = jpcv + 2
    pressure = np.zeros((jpmax, ipmax))
    p_corr = np.zeros((jpmax, ipmax))
    iumax = ipcv + 1
    jumax = jpcv + 2
    u_vel = np.zeros((jumax, iumax))
    ivmax = ipcv + 2
    jvmax = jpcv + 1
    v_vel = np.zeros((jvmax, ivmax))

    # define properties
    rho = properties.DENSITY
    # k = properties.THERMAL_CONDUCTIVITY
    mu = properties.VISCOSITY
    # cp = properties.SPECIFIC_HEAT_CAPACITY

    # set initial conditions
    u_top_ref = re_number_top * mu / (rho * l)
    u_bottom_ref = re_number_bottom * mu / (rho * l)
    u_vel = ics_function(u_vel, u_top_ref, u_bottom_ref)

    # define omega values
    w_u = 0.5
    w_v = 0.5
    w_p = 1.0

    iteration = 0
    u_res = np.zeros(800)  # 800 as the limiting iterations
    v_res = np.zeros(800)
    p_res = np.zeros(800)

    while True:
        # step-1: compute u-velocity coefficients
        u_coeff = u_vel_coeff(u_vel, v_vel, dx, dy, rho, mu, iumax, jumax)

        # step-2: compute u-velocity
        u_vel = u_vel_solve(u_vel, u_coeff, iumax, jumax, w_u, dy, pressure)

        # step-3: compute v-velocity coefficients
        v_coeff = v_vel_coeff(u_vel, v_vel, dx, dy, rho, mu, ivmax, jvmax)

        # step-4: compute v-velocity
        v_vel = v_vel_solve(v_vel, v_coeff, ivmax, jvmax, w_v, dx, pressure)

        # step-5: compute pressure correction coefficients
        p_coeff = p_corr_coeff(dx, dy, rho, u_coeff[0], v_coeff[0], ipmax, jpmax)

        # step-6: solve pressure correction equation
        pressure, p_corr = p_correction(
            pressure, p_corr, p_coeff, ipmax, jpmax, dx, dy, w_p, rho, u_vel, v_vel
        )

        # step-7: solve u-velocity correction equation
        u_vel = u_vel_correction(u_vel, u_coeff[0], p_corr, dy)

        # step-8: solve v-velocity correction equation
        v_vel = v_vel_correction(v_vel, v_coeff[0], p_corr, dx)

        # reset pressure correction matrix
        p_corr = np.zeros_like(p_corr)

        # step-9: check for convergence
        u_res[iteration], v_res[iteration], p_res[iteration] = convergence(
            u_vel, u_coeff, v_vel, v_coeff, pressure, u_top_ref, dx, dy, rho, l
        )
        if not (u_res[iteration] > 1e-6 and v_res[iteration] > 1e-6 and p_res[iteration] > 1e-5):
            break

        # increase the iteration count
        iteration += 1


    return u_vel, v_vel, pressure, u_res, v_res, p_res
