# pressure solve function
import numpy as np
from lid_driven_cavity_python.TDMA_solver import solve_TDMA


def p_correction(
    pressure, p_corr, p_coeff, ipmax, jpmax, dy, dx, w_p, rho, u_vel, v_vel
):
    Ap_p = p_coeff[0]
    Aw_p = p_coeff[1]
    Ae_p = p_coeff[2]
    As_p = p_coeff[3]
    An_p = p_coeff[4]

    # x-sweep
    A = np.zeros((jpmax - 2, 3))

    for i in range(1, ipmax - 1):
        A[:, 0] = -As_p[:, i - 1]
        A[:, 1] = Ap_p[:, i - 1] / w_p
        A[:, 2] = -An_p[:, i - 1]
        b = Aw_p[:, i - 1] * p_corr[1:-1, i - 1]
        b = b + Ae_p[:, i - 1] * p_corr[1:-1, i + 1]
        b = b + Ap_p[:, i - 1] * p_corr[1:-1, i] * ((1 - w_p) / w_p)
        b = b + rho * (
            dy * (u_vel[1:-1, i - 1] - u_vel[1:-1, i])
            + dx * (v_vel[:-1, i] - v_vel[1:, i])
        )
        b[0] = b[0] + As_p[0, i - 1] * p_corr[0, i]
        b[-1] = b[-1] + An_p[-1, i - 1] * p_corr[-1, i]

        p_corr[1:-1, i] = solve_TDMA(jpmax - 2, A, b)

    # y-sweep
    A = np.zeros((ipmax - 2, 3))

    for j in range(1, jpmax - 1):
        A[:, 0] = -Aw_p[j - 1, :]
        A[:, 1] = Ap_p[j - 1, :] / w_p
        A[:, 2] = -Ae_p[j - 1, :]
        b = As_p[j - 1, :] * p_corr[j - 1, 1:-1]
        b = b + An_p[j - 1, :] * p_corr[j + 1, 1:-1]
        b = b + Ap_p[j - 1, :] * p_corr[j, 1:-1] * ((1 - w_p) / w_p)
        b = b + rho * (
            dy * (u_vel[j, :-1] - u_vel[j, 1:])
            + dx * (v_vel[j - 1, 1:-1] - v_vel[j, 1:-1])
        )
        b[0] = b[0] + Aw_p[j - 1, 0] * p_corr[j, 0]
        b[-1] = b[-1] + Ae_p[j - 1, -1] * p_corr[j, -1]

        p_corr[j, 1:-1] = solve_TDMA(ipmax - 2, A, b)

    # calculating pressure values
    pressure = pressure + 0.8 * p_corr

    # setting up boundary values for pressure field
    pressure[1:-1, 0] = pressure[1:-1, 1]
    pressure[1:-1, -1] = pressure[1:-1, -2]
    pressure[0, 1:-1] = pressure[1, 1:-1]
    pressure[-1, 1:-1] = pressure[-2, 1:-1]

    return pressure, p_corr
