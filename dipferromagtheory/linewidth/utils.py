"""
Mathematical utility functions.
"""

def kdelta(i, j):
    """
    Kronecker delta function
    """

    assert i in (0,1)
    assert j in (0,1)

    return i*j

def rhom(rho, eta):
    """

    """
    return (1 + rho**2 - 2*rho*eta)**0.5

# def rhom_old(y, g, rho, eta):
#     """

#     """
#     return g**0.5 / y * (1 - 2*eta*rho + rho**2)