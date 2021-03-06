# residuals calculation function
import numpy as np


def convergence(u_vel, u_coeff, v_vel, v_coeff, pressure, u_ref, dx, dy, rho, L):
    Ap_u = u_coeff[0]
    Aw_u = u_coeff[1]
    Ae_u = u_coeff[2]
    As_u = u_coeff[3]
    An_u = u_coeff[4]

    Ap_v = v_coeff[0]
    Aw_v = v_coeff[1]
    Ae_v = v_coeff[2]
    As_v = v_coeff[3]
    An_v = v_coeff[4]

    # u velocity residual calculation
    numerator = Ap_u * u_vel[1:-1, 1:-1]
    numerator = numerator - (
        Aw_u * u_vel[1:-1, :-2]
        + Ae_u * u_vel[1:-1, 2:]
        + As_u * u_vel[:-2, 1:-1]
        + An_u * u_vel[2:, 1:-1]
    )
    numerator = numerator - dx * (pressure[1:-1, 1:-2] - pressure[1:-1, 2:-1])
    numerator = np.sum(np.sum(np.absolute(numerator)))
    denom = np.sum(np.sum(np.absolute(Ap_u * u_vel[1:-1, 1:-1])))
    u_res = numerator / denom * 100

    # v velocity residual calculation
    numerator = Ap_v * v_vel[1:-1, 1:-1]
    numerator = numerator - (
        Aw_v * v_vel[1:-1, :-2]
        + Ae_v * v_vel[1:-1, 2:]
        + As_v * v_vel[:-2, 1:-1]
        + An_v * v_vel[2:, 1:-1]
    )
    numerator = numerator - dy * (pressure[1:-2, 1:-1] - pressure[2:-1, 1:-1])
    numerator = np.sum(np.sum(np.absolute(numerator)))
    denom = np.sum(np.sum(np.absolute(Ap_v * v_vel[1:-1, 1:-1])))
    v_res = numerator / denom * 100

    # pressure residual calculation
    numerator = rho * (u_vel[1:-1, :-1] - u_vel[1:-1, 1:]) * dy
    numerator = numerator + rho * (v_vel[:-1, 1:-1] - v_vel[1:, 1:-1]) * dx
    numerator = np.sum(np.sum(np.absolute(numerator)))
    denom = rho * u_ref * L
    p_res = numerator / denom * 1000

    return u_res, v_res, p_res
