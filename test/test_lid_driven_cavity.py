import pytest
import re
import numpy as np

from lid_driven_cavity_python.TDMA_solver import solve_TDMA


def test_TDMA_solver():
    """
    Test the matrix solver function solve()
    """
    # define test matrix

    test_a = [[0, -2.6, 1], [1, -2.6, 1], [1, -2.6, 1], [1, -2.6, 0]]
    test_b = [-240.0, 0.0, 0.0, -150.0]

    # define the solution to the test matrix
    correct_x = [118.1122, 67.0916, 56.3261, 79.3562]

    # run the TDMA solver
    x = solve_TDMA(len(test_b), np.asarray(test_a), np.asarray(test_b))

    # test the solution with correct answer
    assert x == pytest.approx(correct_x, abs=1e-4)
