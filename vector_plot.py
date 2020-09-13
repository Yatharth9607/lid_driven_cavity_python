# vector plot with animation
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from lid_driven_cavity_python.cell_center_velocity import cell_velocity

def plot_velocities(u_vel, v_vel):
    """receives face centered velocity values"""
    # step 1: create data to plot
    nx = np.shape(u_vel)[1] - 1
    ny = np.shape(u_vel)[0] - 2
    x, y = np.meshgrid(np.arange(0, nx), np.arange(0, ny))
    U, V = cell_velocity(u_vel, v_vel)

    u_animation = [0] * 10
    v_animation = [0] * 10

    for i in range(10):
        u_animation[i] = U * 0.1 * i
        v_animation[i] = V * 0.1 * i

    # step 2: setup plotting area
    fig, ax = plt.subplots()
    plt.xlim((-0.5, 9.5))
    plt.ylim((-0.5, 9.5))

    # step 3: plot first set of vectors
    vector_plot = ax.quiver(x, y, u_animation[-1], v_animation[-1], units="width", scale=0.005)

    # step 4: create a function to update the vectors
    def animate(i):
        vector_plot.set_UVC(u_animation[i], v_animation[i])
        return vector_plot

    # step 5: call FuncAnimation and show
    anim = FuncAnimation(
        fig, animate, interval=100, frames=10
    )

    # step 6: garnish the plots
    # generate gridlines for cells
    for i in range(nx):
        ax.axvline(i-0.5, color='gray', linestyle='--', linewidth=0.5)
    for j in range(ny):
        ax.axhline(j-0.5, color='gray', linestyle='--', linewidth=0.5)

    # set x and y axis
    plt.xticks(ticks=np.arange(nx), labels=np.arange(1, nx+1))
    plt.yticks(ticks=np.arange(nx), labels=np.arange(1, ny+1))

    # set titles
    ax.set_title("Velocity Vector Plot")
    plt.xlabel("x-direction")
    plt.ylabel("y-direction")

    # tada...
    plt.show()
