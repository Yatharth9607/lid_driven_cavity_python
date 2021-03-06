# Setup pressure correction coefficients function
import numpy as np


def p_corr_coeff(dx, dy, rho, Ap_u, Ap_v, ipmax, jpmax):

    # calculate d_u and d_v
    d_u = dy * np.reciprocal(Ap_u)
    d_v = dx * np.reciprocal(Ap_v)

    # calculate a coefficients
    Aw_p = np.zeros((jpmax - 2, ipmax - 2))
    Ae_p = np.zeros((jpmax - 2, ipmax - 2))
    As_p = np.zeros((jpmax - 2, ipmax - 2))
    An_p = np.zeros((jpmax - 2, ipmax - 2))
    Aw_p[:, :-1] = rho * d_u * dy
    Ae_p[:, 1:] = rho * d_u * dy
    As_p[1:, :] = rho * d_v * dx
    An_p[:-1, :] = rho * d_v * dx

    Ap_p = Aw_p + Ae_p + As_p + An_p

    p_coeff = [Ap_p, Aw_p, Ae_p, As_p, An_p]

    return p_coeff
