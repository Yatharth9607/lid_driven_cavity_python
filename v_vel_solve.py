# v velocity solve function
import numpy as np

from lid_driven_cavity_python.TDMA_solver import solve_TDMA


def v_vel_solve(v_vel, v_coeff, ivmax, jvmax, w_v, dx, pressure):

    Ap_v = v_coeff[0]
    Aw_v = v_coeff[1]
    Ae_v = v_coeff[2]
    As_v = v_coeff[3]
    An_v = v_coeff[4]

    # x-sweep
    A = np.zeros((jvmax - 2, 3))

    for i in range(1, ivmax - 1):
        A[:, 0] = -As_v[:, i - 1]
        A[:, 1] = Ap_v[:, i - 1] / w_v
        A[:, 2] = -An_v[:, i - 1]
        b = Aw_v[:, i - 1] * v_vel[1:-1, i - 1]
        b = b + Ae_v[:, i - 1] * v_vel[1:-1, i + 1]
        b = b + Ap_v[:, i - 1] * v_vel[1:-1, i] * ((1 - w_v) / w_v)
        b = b + dx * (pressure[1:-2, i] - pressure[2:-1, i])
        b[0] = b[0] + As_v[0, i - 1] * v_vel[0, i]
        b[-1] = b[-1] + An_v[-1, i - 1] * v_vel[-1, i]

        v_vel[1:-1, i] = solve_TDMA(jvmax - 2, A, b)

    # y-sweep
    A = np.zeros((ivmax - 2, 3))

    for j in range(1, jvmax - 1):
        A[:, 0] = -Aw_v[j - 1, :]
        A[:, 1] = Ap_v[j - 1, :] / w_v
        A[:, 1] = -Ae_v[j - 1, :]
        b = As_v[j - 1, :] * v_vel[j - 1, 1:-1]
        b = b + An_v[j - 1, :] * v_vel[j + 1, 1:-1]
        b = b + Ap_v[j - 1, :] * v_vel[j, 1:-1] * ((1 - w_v) / w_v)
        b = b + dx * (pressure[j, 1:-1] - pressure[j + 1, 1:-1])
        b[0] = b[0] + Aw_v[j - 1, 0] * v_vel[j, 0]
        b[-1] = b[-1] + Ae_v[j - 1, -1] * v_vel[j, -1]

        v_vel[j][1:-1] = solve_TDMA(ivmax - 2, A, b)

    return v_vel
