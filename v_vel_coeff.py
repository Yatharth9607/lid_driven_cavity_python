# Setup v_vel coefficients function
import numpy as np

from lid_driven_cavity_python.peclet_function import peclet_function


def v_vel_coeff(u_vel, v_vel, dx, dy, rho, gamma, ivmax, jvmax):

    ivcv = ivmax - 2
    jvcv = jvmax - 2

    # calculate physical distance delta_x and delta_y
    delta_x = dx * np.ones(ivcv, jvcv)
    delta_y = dy * np.ones(ivcv, jvcv)
    delta_y[0][:] = 3 * dy / 2
    delta_y[-1][:] = 3 * dy / 2

    # calculate diffusion length del_x and del_y
    del_x_w = dx * np.ones(ivcv, jvcv)
    del_x_w[:][0] = dx / 2
    del_x_e = dx * np.ones(ivcv, jvcv)
    del_x_e[:][-1] = dx / 2
    del_y_s = dy * np.ones(ivcv, jvcv)
    del_y_n = dy * np.ones(ivcv, jvcv)

    # calculate diffusion length array
    Dw = delta_y * gamma / del_x_w
    De = delta_y * gamma / del_x_e
    Ds = delta_x * gamma / del_y_s
    Dn = delta_x * gamma / del_y_n

    # calculate flow strength
    Fw = rho * ((u_vel[1:-3][:-2] + u_vel[2:-2][:-2]) / 2) * delta_y
    Fw[0][:] = rho * (
        ((u_vel[0][:-2] + u_vel[1][:-2]) / 2) * dy / 2
        + ((u_vel[1][:-2] + u_vel[2][:-2]) / 2) * dy
    )
    Fw[-1][:] = rho * (
        ((u_vel[-3][:-2] + u_vel[-2][:-2]) / 2) * dy
        + ((u_vel[-2][:-2] + u_vel[-1][:-2]) / 2) * dy / 2
    )
    Fe = rho * ((u_vel[1:-3][1:] + u_vel[2:-2][1:]) / 2) * delta_y
    Fe[0][:] = rho * (
        ((u_vel[0][1:] + u_vel[1][1:]) / 2) * dy / 2
        + ((u_vel[1][1:] + u_vel[2][1:]) / 2) * dy
    )
    Fe[-1][:] = rho * (
        ((u_vel[-3][1:] + u_vel[-2][1:]) / 2) * dy
        + ((u_vel[-2][1:] + u_vel[-1][1:]) / 2) * dy / 2
    )
    Fs = rho * ((v_vel[:-3][1:-2] + v_vel[1:-2][1:-2]) / 2) * delta_x
    Fs[0][:] = rho * v_vel[0][1:-2] * delta_x[0][:]
    Fn = rho * ((v_vel[1:-2][1:-2] + v_vel[2:][1:-2]) / 2) * delta_x
    Fn[-1][:] = rho * v_vel[-1][1:-2] * delta_x[-1][:]

    # calculate Peclet number
    Pw = Fw * Dw ^ -1
    Pe = Fe * De ^ -1
    Ps = Fs * Ds ^ -1
    Pn = Fn * Dn ^ -1

    # compute coefficient arrays
    Aw_v = Dw * peclet_function(np.absolute(Pw), "power_law") + np.max(
        Fw, np.zeros_like(Fw)
    )
    Ae_v = De * peclet_function(np.absolute(Pe), "power_law") + np.max(
        -Fe, np.zeros_like(Fe)
    )
    As_v = Ds * peclet_function(np.absolute(Ps), "power_law") + np.max(
        Fs, np.zeros_like(Fs)
    )
    An_v = Dn * peclet_function(np.absolute(Pn), "power_law") + np.max(
        -Fn, np.zeros_like(Fn)
    )
    Ap_v = Aw_v + Ae_v + As_v + An_v

    v_coeff = [Ap_v, Aw_v, Ae_v, As_v, An_v]

    return v_coeff
