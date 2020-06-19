# peclet function (A(|P|))
import math
import numpy as np

# parameters required:
# p: Peclet number (array)
# solver_type: type of solver ('central_diff', 'exponential', 'power_law')

# return output:
# a_p: Peclet function (A(|P|)), depending upon the solver type

def peclet_function(p, solver_type):

    # central difference scheme
    if solver_type == 'central_diff':
        a_p = 1 - 0.5 * np.absolute(p)
    
    # exponential scheme
    elif solver_type == 'exponential':
        a_p = np.exp(np.absolute(p)) / (np.exp(np.absolute(p)) - np.ones_like(p))
        
    # power law
    elif solver_type == 'power_law':
        mult = (np.ones_like(p) - 0.1 * np.absolute(p))
        a_p = np.max(np.zeros_like(p), np.power(mult, 5))

    else:
        raise ValueError

    return a_p