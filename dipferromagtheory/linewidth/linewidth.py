"""
Implementation of the integral equations determining the longitudinal
and transverse linewdith.
"""
import numpy as np
#------------------------------------------------------------------------------
from .vertices import *
from .susceptibility import Chi

class DynScalingFunc():
    """

    """
    def __init__(self, vertices, susceptibilities):
        """

        """
        self.vertices = vertices
        self. susceptibilities = susceptibilities


    def __call__(self, x, y):
        """

        """
        pass

    def calc(self, x, y):
        """

        """
        pass

    def verify(self, gammaxy):
        """

        """
        pass

###############################################################################
###############################################################################
###############################################################################

def kappaNi(T, Tc):
    """

    """
    return 0.81 * np.abs((T - Tc)/Tc)**0.701

j=0.1
g = 1.0
xi = 1 / kappaNi(T=629, Tc=628)

q = np.linspace(0.0, 0.03, 101)
gammaq = np.array([5.1326 * ] * len(q))

chi = Chi(j=j, g=g, xi=xi)

lwL = 