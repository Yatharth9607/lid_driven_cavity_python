# Setup v_vel coefficients function
import numpy as np
from lid_driven_cavity_python.peclet_function import peclet_function


def v_vel_coeff(u_vel, v_vel, dx, dy, rho, gamma, ivmax, jvmax):
    ivcv = ivmax - 2
    jvcv = jvmax - 2

    # calculate physical distance delta_x (size = (jvcv, ivcv)) and delta_y (size = (jvcv, ivcv))
    delta_x = dx * np.ones((jvcv, ivcv))
    delta_y = dy * np.ones((jvcv, ivcv))
    delta_y[0][:] = 3 * dy / 2
    delta_y[-1][:] = 3 * dy / 2

    # calculate diffusion length del_x (size = (jvcv, ivcv)) and del_y (size = (jvcv, ivcv))
    del_x_w = dx * np.ones((jvcv, ivcv))
    del_x_w[:][0] = dx / 2
    del_x_e = dx * np.ones((jvcv, ivcv))
    del_x_e[:][-1] = dx / 2
    del_y_s = dy * np.ones((jvcv, ivcv))
    del_y_n = dy * np.ones((jvcv, ivcv))

    # calculate diffusion length array (size = (jvcv, ivcv))
    Dw = delta_y * gamma / del_x_w
    De = delta_y * gamma / del_x_e
    Ds = delta_x * gamma / del_y_s
    Dn = delta_x * gamma / del_y_n

    # calculate flow strength (size = (jvcv, ivcv))
    Fw = rho * ((u_vel[1:-2, :-1] + u_vel[2:-1, :-1]) / 2) * delta_y
    Fw[0, :] = rho * (
        ((u_vel[0, :-1] + u_vel[1, :-1]) / 2) * dy / 2
        + ((u_vel[1, :-1] + u_vel[2, :-1]) / 2) * dy
    )
    Fw[-1, :] = rho * (
        ((u_vel[-3, :-1] + u_vel[-2, :-1]) / 2) * dy
        + ((u_vel[-2, :-1] + u_vel[-1, :-1]) / 2) * dy / 2
    )
    Fe = rho * ((u_vel[1:-2, 1:] + u_vel[2:-1, 1:]) / 2) * delta_y
    Fe[0, :] = rho * (
        ((u_vel[0, 1:] + u_vel[1, 1:]) / 2) * dy / 2
        + ((u_vel[1, 1:] + u_vel[2, 1:]) / 2) * dy
    )
    Fe[-1, :] = rho * (
        ((u_vel[-3, 1:] + u_vel[-2, 1:]) / 2) * dy
        + ((u_vel[-2, 1:] + u_vel[-1, 1:]) / 2) * dy / 2
    )
    Fs = rho * ((v_vel[:-2, 1:-1] + v_vel[1:-1, 1:-1]) / 2) * delta_x
    Fs[0, :] = rho * v_vel[0, 1:-1] * delta_x[0, :]
    Fn = rho * ((v_vel[1:-1, 1:-1] + v_vel[2:, 1:-1]) / 2) * delta_x
    Fn[-1, :] = rho * v_vel[-1, 1:-1] * delta_x[-1, :]

    # calculate Peclet number
    Pw = Fw * np.reciprocal(Dw)
    Pe = Fe * np.reciprocal(De)
    Ps = Fs * np.reciprocal(Ds)
    Pn = Fn * np.reciprocal(Dn)

    # compute coefficient arrays
    Aw_v = Dw * peclet_function(np.absolute(Pw), "power_law") + np.maximum(
        Fw, np.zeros_like(Fw)
    )
    Ae_v = De * peclet_function(np.absolute(Pe), "power_law") + np.maximum(
        -Fe, np.zeros_like(Fe)
    )
    As_v = Ds * peclet_function(np.absolute(Ps), "power_law") + np.maximum(
        Fs, np.zeros_like(Fs)
    )
    An_v = Dn * peclet_function(np.absolute(Pn), "power_law") + np.maximum(
        -Fn, np.zeros_like(Fn)
    )
    Ap_v = Aw_v + Ae_v + As_v + An_v

    v_coeff = [Ap_v, Aw_v, Ae_v, As_v, An_v]

    return v_coeff
