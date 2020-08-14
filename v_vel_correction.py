# v velocity correction function
import numpy as np


def v_vel_correction(v_vel, Ap_v, p_corr, dx):

    v_vel_corr = (dx * np.reciprocal(Ap_v)) * (p_corr[1:-2, 1:-1] - p_corr[2:-1, 1:-1])
    v_vel[1:-1, 1:-1] = v_vel[1:-1, 1:-1] + v_vel_corr
    return v_vel
