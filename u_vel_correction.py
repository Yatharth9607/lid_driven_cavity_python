# u velocity correction function


def u_vel_correction(u_vel, Ap_u, p_corr, dy):

    u_vel_corr = (dy / Ap_u) * (p_corr[1:-2][1:-3] - p_corr[1:-2][2:-2])
    u_vel[1:-2][1:-2] = u_vel[1:-2][1:-2] + u_vel_corr
    return u_vel
