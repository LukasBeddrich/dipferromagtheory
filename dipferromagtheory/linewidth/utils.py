"""
Mathematical utility functions.
"""

import numpy as np

def kdelta(i, j):
    """
    Kronecker delta function
    """

    assert i in (0,1)
    assert j in (0,1)

    return i*j

def rhom(rho, eta):
    """

    Note
    ----
    This becomes 0 if eta = 1 and rho = 1, which leads to trouble
    when calculating rhom(rho, eta)**-2
    """
    return (1 + rho**2 - 2*rho*eta)**0.5

# def rhom_old(y, g, rho, eta):
#     """

#     """
#     return g**0.5 / y * (1 - 2*eta*rho + rho**2)

def sum_of_squared_difference(arr):
    res = np.zeros(len(arr) - 1)
    for idx, subarr in enumerate(arr[:-1]):
        res[idx] = np.power(subarr - arr[idx+1], 2).sum()
    return res