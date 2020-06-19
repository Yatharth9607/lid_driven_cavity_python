# Matrix solver (TDMA)

# parameters required:
# n: number of unknowns (length of x)
# a: coefficient matrix
# b: RHS/constant array

# return output:
# b: solution array

def solve_TDMA(n, a, b):
    # forward substitution
    a[2][0] = -a[2][0] / a[1][0]
    b[0] = b[0] / a[1][0]
    for i in range(1,n):
        a[2][i] = -a[2][i] / (a[1][i] + a[0][i]*a[2][i-1])
        b[i] = (b[i] - a[0][i]*b[i-1]) / (a[1][i] + a[0][i]*a[2][i-1])
    
    # backward elimination
    for i in range(n-2,-1,-1):
        b[i] = a[2][i]*b[i+1] + b[i]
    return b