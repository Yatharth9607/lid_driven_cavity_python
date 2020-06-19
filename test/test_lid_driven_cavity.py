import pytest
import re
import numpy as np

from lid_driven_cavity_python.TDMA_solver import solve_TDMA


def test_TDMA_solver():
    """
    Test the matrix solver function solve()
    """
    # define test matrix
    test_a = [[0, 1, 1, 1],
              [-2.6, -2.6, -2.6, -2.6],
              [1, 1, 1, 0]]
    test_b = [-240, 0, 0, -150]

    # define the solution to the test matrix
    correct_x = [118.1122, 67.0916, 56.3261, 79.3562]

    # run the TDMA solver
    x = solve_TDMA(len(test_b), test_a, test_b)

    # test the solution with correct answer
    assert((x == pytest.approx(correct_x, abs=1e-4)))
    