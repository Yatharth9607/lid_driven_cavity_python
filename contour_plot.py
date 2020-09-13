# contour plot with matplotlib
import numpy as np
import matplotlib.pyplot as plt

def plot_pressure(pressure):
    """receives cell centered pressure values"""
    # step 1: create data to plot
    nx = np.shape(pressure)[1] - 2
    ny = np.shape(pressure)[0] - 2
    x_axis = np.arange(-1.0, nx+1)
    x_axis[0] += 0.5
    x_axis[-1] -= 0.5
    y_axis = np.arange(-1.0, ny+1)
    y_axis[0] += 0.5
    y_axis[-1] -= 0.5
    x, y = np.meshgrid(x_axis, y_axis)

    # step 2: setup plotting area
    fig1, ax1 = plt.subplots()
    plt.xlim((-0.5, 9.5))
    plt.ylim((-0.5, 9.5))

    # step 3: plot first set of contours
    plt.contourf(x, y, pressure, levels=100)

    # step 6: garnish the plots
    # generate gridlines for cells
    for i in range(nx):
        ax1.axvline(i-0.5, color='gray', linestyle='--', linewidth=0.5)
    for j in range(ny):
        ax1.axhline(j-0.5, color='gray', linestyle='--', linewidth=0.5)
    plt.colorbar(format="%.5f")
    
    # set x and y axis
    plt.xticks(ticks=np.arange(nx), labels=np.arange(1, nx+1))
    plt.yticks(ticks=np.arange(nx), labels=np.arange(1, ny+1))

    # set titles
    ax1.set_title("Pressure Contour Plot")
    plt.xlabel("x-direction")
    plt.ylabel("y-direction")

    # tada...
    plt.show()
