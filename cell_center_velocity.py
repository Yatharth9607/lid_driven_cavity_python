# compute velocities at cell center from face center values
import numpy as np


def cell_velocity(u_vel_face, v_vel_face):
    nx = np.shape(v_vel_face)[0] - 1
    ny = np.shape(u_vel_face)[1] - 1

    u_vel_cell = np.zeros((nx, ny))
    v_vel_cell = np.zeros((nx, ny))

    u_vel_cell[:, :] = (u_vel_face[1:-1, :-1] + u_vel_face[1:-1, 1:]) / 2
    v_vel_cell[:, :] = (v_vel_face[:-1, 1:-1] + v_vel_face[1:, 1:-1]) / 2

    # u_vel_cell = (u_vel_face[1:-1, :-1] + u_vel_face[1:-1, 1:]) / 2
    # v_vel_cell = (v_vel_face[:-1, 1:-1] + v_vel_face[1:, 1:-1]) / 2
    return u_vel_cell, v_vel_cell
