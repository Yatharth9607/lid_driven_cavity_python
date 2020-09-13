# line plot with matplotlib
import numpy as np
import matplotlib.pyplot as plt

def plot_convergence(u_res, v_res, p_res):
    """receives total residual values"""
    # step 1: create data to plot
    p = np.count_nonzero(u_res)
    q = np.count_nonzero(v_res)
    r = np.count_nonzero(p_res)
    iterations = max(p, q, r)
    x_axis = np.arange(iterations)

    # step 2: setup plotting area
    fig1, ax1 = plt.subplots()

    # step 3: plot first set of contours
    plt.plot(x_axis, u_res[:iterations], 'r--', x_axis, v_res[:iterations], 'b--', x_axis, p_res[:iterations], 'g--')

    # step 6: garnish the plots
    labels = ["x-velocity", "y-velocity", "pressure"]
    plt.legend(labels)
    
    # set titles
    ax1.set_title("Residual Line Plot")
    plt.xlabel("Iterations")
    plt.ylabel("Residuals")

    # tada...
    plt.show()
