# u velocity correction function
import numpy as np


def u_vel_correction(u_vel, Ap_u, p_corr, dy):

    u_vel_corr = (dy * np.reciprocal(Ap_u)) * (p_corr[1:-1, 1:-2] - p_corr[1:-1, 2:-1])
    u_vel[1:-1, 1:-1] = u_vel[1:-1, 1:-1] + u_vel_corr
    return u_vel
