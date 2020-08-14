# Matrix solver (TDMA)

# parameters required:
# n: number of unknowns (length of x)
# a: coefficient matrix
# b: RHS/constant array

# return output:
# b: solution array


def solve_TDMA(n, a, b):
    # forward substitution
    a[0][2] = -a[0][2] / a[0][1]
    b[0] = b[0] / a[0][1]
    for i in range(1, n):
        a[i][2] = -a[i][2] / (a[i][1] + a[i][0] * a[i - 1][2])
        b[i] = (b[i] - a[i][0] * b[i - 1]) / (a[i][1] + a[i][0] * a[i - 1][2])
    # backward elimination
    for i in range(n - 2, -1, -1):
        b[i] = a[i][2] * b[i + 1] + b[i]
    return b
