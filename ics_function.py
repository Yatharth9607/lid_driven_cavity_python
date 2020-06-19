# initial boundary conditions

def ics_function(u_vel, u_top_ref, u_bottom_ref):
    u_vel[0] = u_top_ref
    u_vel[-1] = u_bottom_ref
    u_vel[1:-2][1:-2] = 1e-4    # to have some momentum

    return u_vel