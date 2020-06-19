# compute u-velocity coefficients
import numpy as np

def u_vel_coeff(u_vel, v_vel, dx, dy, rho, omega, iumax, jumax):
    iucv = iumax - 2
    jucv = jumax - 2

    # calculate physical distance delta_x and delta_y
    delta_x = dx * np.ones(iucv, jucv)
    delta_x[:][0] = (3*dx/2)
    delta_x[:][-1] = (3*dx/2)
    delta_y = dy * np.ones(iucv, jucv)
    
    # calculate diffusion length del_x and del_y
    del_x_w = dx * np.ones(iucv, jucv)
    del_x_e = dx * np.ones(iucv, jucv)
    del_y_s = dy * np.ones(iucv, jucv)
    del_y_s[0][:] = dy/2
    del_y_n = dy * np.ones(iucv, jucv)
    del_y_n[-1][:] = dy/2
    
    # calculate diffusion length array
    Dw = (delta_y * gamma / del_x_w)
    De = (delta_y * gamma / del_x_e)
    Ds = (delta_x * gamma / del_y_s)
    Dn = (delta_x * gamma / del_y_n)
    
    # calculate flow strength
    Fw = rho * ((u_vel[1:-2][:-3] + u_vel[1:-2][1:-2]) / 2) * delta_y
    Fw[:][0] = rho * u_vel[1:-2][0] * delta_y[:][0]
    Fe = rho * ((u_vel[1:-2][1:-2] + u_vel[1:-2][2:-1]) / 2) * delta_y
    Fe[:][-1] = rho * u_vel[1:-2][-1] * delta_y[:][0]
    Fs = rho * ((v_vel[:-2][1:-3] + v_vel[:-2][2:-2]) / 2) * delta_x
    Fs[:][0] = rho * (((v_vel[:-2][0] + v_vel[:-2][1]) / 2) * dx/2 + ((v_vel[:-2][1] + v_vel[:-2][2]) / 2) * dx)
    Fs[:][-1] = rho * (((v_vel[:-2][-2] + v_vel[:-2][-1]) / 2) * dx/2 + ((v_vel[:-2][-3] + v_vel[:-2][-2]) / 2) * dx)
    Fn = rho * ((v_vel[1:-1][1:-3] + v_vel[1:-1][2:-2]) / 2) * delta_x
    Fn[:][0] = rho * (((v_vel[1:-1][0] + v_vel[1:-1][1]) / 2) * dx/2 + ((v_vel[1:-1][1] + v_vel[1:-1][2]) / 2) * dx)
    Fn[:][-1] = rho * (((v_vel[1:-1][-2] + v_vel[1:-1][-1]) / 2) * dx/2 + ((v_vel[1:-1][-3] + v_vel[1:-1][-2]) / 2) * dx)
    
    # calculate Peclet number
    Pw = Fw * np.inverse(Dw)
    Pe = Fe * np.inverse(De)
    Ps = Fs * np.inverse(Ds)
    Pn = Fn * np.inverse(Dn)

    # compute coefficient arrays
    Aw_u = Dw*peclet_function(np.absolute(Pw), "power_law") + np.maximum(Fw, np.zeros_like(Fw))
    Ae_u = De*peclet_function(np.absolute(Pe), "power_law") + np.maximum(-Fe, np.zeros_like(Fe))
    As_u = Ds*peclet_function(np.absolute(Ps), "power_law") + np.maximum(Fs, np.zeros_like(Fs))
    An_u = Dn*peclet_function(np.absolute(Pn), "power_law") + np.maximum(-Fn, np.zeros_like(Fn))
    Ap_u = Aw_u + Ae_u + As_u + An_u
    
    u_coeff = [Ap_u, Aw_u, Ae_u, As_u, An_u]

    return u_coeff