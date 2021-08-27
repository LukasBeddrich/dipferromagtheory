"""
Implements all the vertex functions.
"""

import numpy as np

from .utils import kdelta, rhom

def vLTT(k, q, g, eta):
    """

    """
    return 2.0 * np.power(eta, 2) * np.power(q, 4) * np.power(k * eta / q - 0.5, 2)

#------------------------------------------------------------------------------

def vLLT(k, q, g, eta):
    """

    """
    return 2.0 * (1 - np.power(eta, 2)) * np.power(q, 4) * np.power(k * eta / q - 0.5 + g / 2 / q**2, 2)

#------------------------------------------------------------------------------

def vTTT(k, q, g, eta):
    """

    """
    return (1 - np.power(eta, 2)) * np.power(q, 4) * (1 + 0.5 * np.power(q / (q**2 + k**2 - q * k * eta), 2)) * np.power(k * eta / q - 0.5, 2)

#------------------------------------------------------------------------------

def vTLT(k, q, g, eta):
    """

    """
    return np.power(q, 4) * (2.0 - (1 - np.power(eta, 2)) * (1 + np.power(q / (q**2 + k**2 - q * k * eta), 2))) * np.power(k * eta / q - 0.5 + g * 0.5 / q**2, 2)

#------------------------------------------------------------------------------

def vTLL(k, q, g, eta):
    """

    """
    return (1 - np.power(eta, 2)) \
        * np.power(q, 4) \
        * 0.5 * np.power(q / (q**2 + k**2 - q * k * eta), 2) \
        * np.power(k * eta / q - 0.5, 2)

###############################################################################

def vLTTscaled(y, g, rho, eta):
    """

    """
    return np.power(rho * eta - 0.5, 2) * 2.0 * np.power(eta, 2)

#------------------------------------------------------------------------------

def vLLTscaled(y, g, rho, eta):
    """

    """
    return np.power(rho * eta - 0.5 * (1 + np.power(y, 2)), 2) * 2 * (1 - np.power(eta, 2))

#------------------------------------------------------------------------------

def vTTTscaled(y, g, rho, eta):
    """

    """
    return np.power(rho * eta - 0.5, 2) * (1 - np.power(eta, 2)) * (1 + 0.5 * np.power(rhom(y, g, rho, eta), -2))

#------------------------------------------------------------------------------

def vTLTscaled(y, g, rho, eta):
    """

    """
    return (2 - (1 - np.power(eta, 2)) * (1 + np.power(rhom(y, g, rho, eta), -2))) * np.power(rho * eta - 0.5 * (1 + np.power(y, 2)), 2)

#------------------------------------------------------------------------------

def vTLLscaled(y, g, rho, eta):
    """

    """
    return (1 - np.power(eta, 2)) * 0.5 * np.power(rhom(y, g, rho, eta), -2) * np.power(rho * eta - 0.5, 2)

###############################################################################