# v velocity correction function

def v_vel_correction(v_vel, Ap_v, p_corr, dx):

    v_vel_corr = (dx / Ap_v) * (p_corr[1:-3][1:-2] - p_corr[2:-2][1:-2])
    v_vel[1:-2][1:-2] = v_vel[1:-2][1:-2] + v_vel_corr
    return v_vel