# u velocity solve function
import numpy as np

from lid_driven_cavity_python.TDMA_solver import solve_TDMA


def u_vel_solve(u_vel, u_coeff, iumax, jumax, w_u, dy, pressure):

    Ap_u = u_coeff[0]
    Aw_u = u_coeff[1]
    Ae_u = u_coeff[2]
    As_u = u_coeff[3]
    An_u = u_coeff[4]

    # x-sweep
    A = np.zeros((jumax - 2, 3))
    for i in range(1, iumax - 1):
        A[:, 0] = -As_u[:, i - 1]
        A[:, 1] = Ap_u[:, i - 1] / w_u
        A[:, 2] = -An_u[:, i - 1]
        b = Aw_u[:, i - 1] * u_vel[1:-1, i - 1]
        b = b + Ae_u[:, i - 1] * u_vel[1:-1, i + 1]
        b = b + Ap_u[:, i - 1] * u_vel[1:-1, i] * ((1 - w_u) / w_u)
        b = b + dy * (pressure[1:-1, i] - pressure[1:-1, i + 1])
        b[0] = b[0] + As_u[0, i - 1] * u_vel[0, i]
        b[-1] = b[-1] + An_u[-1, i - 1] * u_vel[-1, i]

        u_vel[1:-1, i] = solve_TDMA(jumax - 2, A, b)

    # y-sweep
    A = np.zeros((iumax - 2, 3))

    for j in range(1, jumax - 1):
        A[:, 0] = -Aw_u[j - 1, :]
        A[:, 1] = Ap_u[j - 1, :] / w_u
        A[:, 2] = -Ae_u[j - 1, :]
        b = As_u[j - 1, :] * u_vel[j - 1, 1:-1]
        b = b + An_u[j - 1, :] * u_vel[j + 1, 1:-1]
        b = b + Ap_u[j - 1, :] * u_vel[j, 1:-1] * ((1 - w_u) / w_u)
        b = b + dy * (pressure[j, 1:-2] - pressure[j, 2:-1])
        b[0] = b[0] + Aw_u[j - 1, 0] * u_vel[j, 0]
        b[-1] = b[-1] + Ae_u[j - 1, -1] * u_vel[j, -1]

        u_vel[j, 1:-1] = solve_TDMA(iumax - 2, A, b)

    return u_vel
